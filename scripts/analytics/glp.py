#!/usr/bin/env python3
import glob
import os
import re
import time
from datetime import datetime, timedelta
from rich.live import Live
from rich.table import Table

# --- Global Simulation Constants ---
TOTAL_NS = 100.0  # Target simulation time in ns
DT_PS = 0.002     # 2 fs timestep = 0.002 ps
TOTAL_STEPS = int((TOTAL_NS * 1000) / DT_PS)  # 50,000,000 steps
STALE_THRESHOLD_MINUTES = 30

# --- Global State Tracker for Delta-based Velocity ---
class SimulationStateTracker:
    def __init__(self):
        # Key: sim_name -> Dict holding previous measurements and smoothed EMA speed
        self.history = {}

    def update_velocity(self, name: str, current_step: int | None) -> tuple[float | None, float | None]:
        """
        Calculates instantaneous ns/day and remaining seconds based on the step delta 
        between live refresh ticks. Applies Exponential Moving Average (EMA) smoothing.
        """
        now_ts = time.time()

        if current_step is None:
            return None, None

        if name not in self.history:
            self.history[name] = {
                "last_step": current_step,
                "last_ts": now_ts,
                "ema_ns_per_sec": None
            }
            return None, None

        prev = self.history[name]
        step_delta = current_step - prev["last_step"]
        time_delta = now_ts - prev["last_ts"]

        # Only update if time has elapsed
        if time_delta > 1.0:
            if step_delta > 0:
                inst_ns_sec = ((step_delta * DT_PS) / 1000.0) / time_delta
                
                # Apply Exponential Moving Average (alpha = 0.3) for smooth reading
                if prev["ema_ns_per_sec"] is None:
                    prev["ema_ns_per_sec"] = inst_ns_sec
                else:
                    prev["ema_ns_per_sec"] = (0.3 * inst_ns_sec) + (0.7 * prev["ema_ns_per_sec"])
            elif step_delta == 0 and prev["ema_ns_per_sec"] is not None:
                # If step hasn't updated on this tick (due to output buffering), decay slightly
                pass

            prev["last_step"] = current_step
            prev["last_ts"] = now_ts

        ema_ns_sec = prev.get("ema_ns_per_sec")
        if ema_ns_sec and ema_ns_sec > 0:
            ns_per_day = ema_ns_sec * 86400.0
            current_ns = (current_step * DT_PS) / 1000.0
            remaining_ns = max(0.0, TOTAL_NS - current_ns)
            remaining_sec = remaining_ns / ema_ns_sec
            return ns_per_day, remaining_sec

        return None, None


# Global tracker instance persisting across Live refresh loops
STATE_TRACKER = SimulationStateTracker()


def read_tail(file_path: str, max_bytes: int = 65536) -> list[str]:
    """Reads the tail end of a text/log file cleanly without high memory overhead."""
    try:
        file_size = os.path.getsize(file_path)
        with open(file_path, "rb") as f:
            if file_size > max_bytes:
                f.seek(file_size - max_bytes)
            raw_bytes = f.read()
        return raw_bytes.decode("utf-8", errors="ignore").splitlines()
    except Exception:
        return []


def format_eta(remaining_seconds: float) -> str:
    """Formats calculated time remaining into clean clock or calendar targets."""
    if remaining_seconds <= 0:
        return "[bold green]Done[/bold green]"
    
    now = datetime.now()
    eta_dt = now + timedelta(seconds=remaining_seconds)
    
    if eta_dt.date() == now.date():
        return eta_dt.strftime("Today %H:%M")
    elif eta_dt.date() == (now + timedelta(days=1)).date():
        return eta_dt.strftime("Tomorrow %H:%M")
    else:
        return eta_dt.strftime("%b %d %H:%M")


def format_gmx_will_finish(will_finish_str: str) -> str:
    """Standardizes raw GROMACS 'will finish' dates into clean display strings."""
    if not will_finish_str:
        return "--"
    try:
        dt = datetime.strptime(will_finish_str.strip(), "%a %b %d %H:%M:%S %Y")
        now = datetime.now()
        if dt.date() == now.date():
            return dt.strftime("Today %H:%M:%S")
        elif dt.date() == (now + timedelta(days=1)).date():
            return dt.strftime("Tomorrow %H:%M:%S")
        else:
            return dt.strftime("%b %d %H:%M:%S")
    except ValueError:
        return will_finish_str.strip()


def parse_local_log(log_path: str) -> tuple[int | None, str | None]:
    """
    Parses local_production_mdrun.log for the latest frame/step index 
    and raw GROMACS 'will finish' estimates.
    """
    lines = read_tail(log_path, max_bytes=65536)
    if not lines:
        return None, None

    last_step = None
    will_finish_raw = None

    for line in reversed(lines):
        # 1. Look for GROMACS 'step X, will finish...' printouts
        step_finish_match = re.search(r"step\s+(\d+),\s+will finish\s+(.*)", line, re.IGNORECASE)
        if step_finish_match:
            if last_step is None:
                last_step = int(step_finish_match.group(1))
            if will_finish_raw is None:
                will_finish_raw = step_finish_match.group(2).strip()

        # 2. Alternative search for step & ns progress
        if last_step is None:
            step_m = re.search(r"Step\s+(\d+)\s+Time\s+([\d\.]+)", line, re.IGNORECASE)
            if step_m:
                last_step = int(step_m.group(1))

        if last_step is not None and will_finish_raw is not None:
            break

    return last_step, will_finish_raw


def parse_md_logs() -> Table:
    table = Table(title="Live GROMACS MD Simulation Tracker (100 ns Runs)", expand=True)
    table.add_column("System / Replicate", style="cyan", no_wrap=True)
    table.add_column("Target", style="bold white", no_wrap=True)
    table.add_column("First Started", style="blue")
    table.add_column("Last Active", style="dim white")
    table.add_column("Sim Time", style="green")
    table.add_column("Progress (Steps / Total)", style="yellow")
    table.add_column("Performance", style="magenta")
    table.add_column("ETA (Est.)", style="bold yellow")
    table.add_column("Will Finish (GMX)", style="bold cyan")

    step_logs = sorted(glob.glob("results/gromacs/*/*/*/standard_100ns/JOB/simulation_steps.log"))

    if not step_logs:
        table.add_row("No simulations found", "-", "-", "-", "-", "-", "-", "-", "-")
        return table

    now = datetime.now()

    for log_path in step_logs:
        try:
            lines = read_tail(log_path, max_bytes=32768)
            if not lines:
                continue

            # 1. Standardize System & Replicate Naming
            parts = log_path.split(os.sep)
            if "gromacs" in parts:
                idx = parts.index("gromacs")
                pdb, source, model_id = parts[idx+1], parts[idx+2], parts[idx+3]
                name = f"{pdb}/{source}/{model_id}"
            else:
                name = log_path
                pdb, source, model_id = None, None, None

            # 2. Extract Timestamps
            first_match = re.search(r"^\[(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2})\]", lines[0])
            first_start_dt = (
                datetime.strptime(first_match.group(1), "%Y-%m-%d %H:%M:%S")
                if first_match
                else datetime.fromtimestamp(os.path.getctime(log_path))
            )
            first_start_str = first_start_dt.strftime("%b %d %H:%M")

            latest_line = lines[-1]
            last_match = re.search(r"^\[(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2})\]", latest_line)
            last_active_dt = (
                datetime.strptime(last_match.group(1), "%Y-%m-%d %H:%M:%S")
                if last_match
                else datetime.fromtimestamp(os.path.getmtime(log_path))
            )
            last_active_str = last_active_dt.strftime("%b %d %H:%M")
            is_stale = (now - last_active_dt).total_seconds() > (STALE_THRESHOLD_MINUTES * 60)

            # 3. Detect Execution Target (HPC vs Local)
            full_log_text = "\n".join(lines)
            is_hpc = bool(
                re.search(r"EXEC_MODE.*HPC", full_log_text, re.IGNORECASE) or
                re.search(r"SLURM", full_log_text, re.IGNORECASE) or
                "komondor" in full_log_text.lower()
            )
            target_str = "[cyan]HPC[/cyan]" if is_hpc else "[green]Local[/green]"

            # 4. Step Code & Queue Context Search (Bottom-Up)
            step_match = re.search(r"\[STEP:\s*([^\]]+)\]\s*(.*)", latest_line)
            step_code = step_match.group(1) if step_match else "UNKNOWN"
            step_details = step_match.group(2) if step_match else ""

            hb_step_val = None
            hb_ns = None
            hb_dt = None
            gmx_will_finish_raw = None
            is_running = False

            for line in reversed(lines):
                if "running -" in line:
                    is_running = True

                    step_finish_match = re.search(
                        r"step\s+(\d+),\s+will finish\s+(.*)", line, re.IGNORECASE
                    )
                    if step_finish_match:
                        if hb_step_val is None:
                            hb_step_val = int(step_finish_match.group(1))
                            hb_ns = (hb_step_val * DT_PS) / 1000.0
                        if gmx_will_finish_raw is None:
                            gmx_will_finish_raw = step_finish_match.group(2).strip()

                    ns_match = re.search(r"running\s*-\s*([\d\.]+)\s*ns", line, re.IGNORECASE)
                    if ns_match and hb_ns is None:
                        hb_ns = float(ns_match.group(1))
                        hb_step_val = int((hb_ns * 1000.0) / DT_PS)

                    if hb_dt is None:
                        dt_m = re.search(r"^\[(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2})\]", line)
                        if dt_m:
                            hb_dt = datetime.strptime(dt_m.group(1), "%Y-%m-%d %H:%M:%S")

                    if hb_step_val is not None and gmx_will_finish_raw is not None:
                        break

            is_queued = (not is_running) and ("State: PD" in step_details or "PD" in step_code)

            # 5. Local Production Log Search
            local_step = None
            if not is_hpc and pdb and source and model_id:
                local_log_path = f"logs/{pdb}/{source}/{model_id}/local_production_mdrun.log"
                if not os.path.exists(local_log_path):
                    fallback_logs = glob.glob(f"logs/{pdb}/{source}/{model_id}/*local*.log")
                    local_log_path = fallback_logs[0] if fallback_logs else None

                if local_log_path and os.path.exists(local_log_path):
                    l_step, l_wf_raw = parse_local_log(local_log_path)
                    local_step = l_step
                    if l_wf_raw and gmx_will_finish_raw is None:
                        gmx_will_finish_raw = l_wf_raw

            # 6. GROMACS Output Directory Log Search
            sim_dir = os.path.dirname(log_path)
            gmx_logs = [
                f for f in glob.glob(os.path.join(sim_dir, "*_md.log")) + glob.glob(os.path.join(sim_dir, "*.log"))
                if os.path.basename(f) != "simulation_steps.log" and not os.path.basename(f).startswith(("cg", "st"))
            ]

            gmx_ns = None
            gmx_step_val = None

            if gmx_logs and not is_queued:
                for glog in gmx_logs:
                    glines = read_tail(glog, max_bytes=65536)
                    for idx, line in enumerate(reversed(glines)):
                        if gmx_will_finish_raw is None:
                            wf_match = re.search(r"will finish\s+(.*)", line, re.IGNORECASE)
                            if wf_match:
                                gmx_will_finish_raw = wf_match.group(1).strip()

                        if gmx_ns is None and "Step" in line and "Time" in line:
                            actual_idx = len(glines) - 1 - idx
                            if actual_idx + 1 < len(glines):
                                val_line = glines[actual_idx + 1].strip()
                                parts_val = val_line.split()
                                if len(parts_val) >= 2 and parts_val[0].isdigit():
                                    try:
                                        gmx_step_val = int(parts_val[0])
                                        gmx_ns = float(parts_val[1]) / 1000.0
                                    except ValueError:
                                        pass
                        if gmx_step_val is not None:
                            break

            # 7. Priority Metrics Consolidation & Delta-Velocity Update
            current_step = local_step if local_step is not None else (gmx_step_val if gmx_step_val is not None else hb_step_val)
            current_ns = (current_step * DT_PS / 1000.0) if current_step is not None else (gmx_ns if gmx_ns is not None else hb_ns)

            # Compute real-time velocity between script refresh cycles
            delta_ns_day, delta_rem_sec = STATE_TRACKER.update_velocity(name, current_step)

            # Format Progress string as Steps / Total Steps (%)
            if current_step is not None:
                pct = min(100.0, (current_step / TOTAL_STEPS) * 100.0)
                progress_str = f"{current_step:,} / {TOTAL_STEPS:,} ({pct:.1f}%)"
            else:
                progress_str = f"0 / {TOTAL_STEPS:,} (0.0%)"

            will_finish_str = format_gmx_will_finish(gmx_will_finish_raw)

            # 8. State & Row Rendering
            if is_queued:
                job_match = re.search(r"Job\s*(\d+)", step_details, re.IGNORECASE)
                job_str = f"Job {job_match.group(1)}" if job_match else "Queued"
                sim_time_str = f"0.0 / {TOTAL_NS:.0f} ns"
                progress_str = f"0 / {TOTAL_STEPS:,} (0.0%)"
                perf_str = f"Queued ({job_str})"
                eta_str = "[yellow]Pending (PD)[/yellow]"
                will_finish_str = "--"

            elif step_code in ["DONE", "SIMULATION_SUMMARY"] or (current_ns and current_ns >= TOTAL_NS):
                sim_time_str = f"{TOTAL_NS:.1f} / {TOTAL_NS:.0f} ns"
                progress_str = f"{TOTAL_STEPS:,} / {TOTAL_STEPS:,} (100.0%)"
                perf_str = "[green]Finished[/green]"
                eta_str = "[bold green]Done[/bold green]"
                will_finish_str = "[bold green]Done[/bold green]"

            elif delta_ns_day is not None and not is_stale:
                sim_ns = current_ns if current_ns is not None else 0.0
                sim_time_str = f"{sim_ns:.2f} / {TOTAL_NS:.0f} ns"
                perf_str = f"{delta_ns_day:.1f} ns/day"
                eta_str = format_eta(delta_rem_sec)

            elif current_ns is not None:
                sim_time_str = f"{current_ns:.2f} / {TOTAL_NS:.0f} ns"
                ref_dt = hb_dt if hb_dt else last_active_dt
                elapsed_sec = (ref_dt - first_start_dt).total_seconds()

                if elapsed_sec > 60 and current_ns > 0 and not is_stale:
                    cumulative_ns_day = (current_ns / elapsed_sec) * 86400.0
                    remaining_sec = ((TOTAL_NS - current_ns) / current_ns) * elapsed_sec
                    perf_str = f"{cumulative_ns_day:.1f} ns/day (avg)"
                    eta_str = format_eta(remaining_sec)
                elif is_stale:
                    perf_str = "[bold yellow]STALLED[/bold yellow]"
                    eta_str = "Inactive"
                    will_finish_str = "Inactive"
                else:
                    perf_str = "Measuring..."
                    eta_str = "Calculating..."

            else:
                sim_time_str = f"0.0 / {TOTAL_NS:.0f} ns"
                progress_str = f"0 / {TOTAL_STEPS:,} (0.0%)"

                if step_code == "RETRIEVE_START":
                    perf_str = "Syncing"
                    eta_str = "Downloading..."
                elif step_code == "ERROR":
                    perf_str = "[bold red]FAILED[/bold red]"
                    eta_str = "[red]Check Log[/red]"
                else:
                    perf_str = "Initializing"
                    eta_str = "Starting..."

            table.add_row(
                name,
                target_str,
                first_start_str,
                last_active_str,
                sim_time_str,
                progress_str,
                perf_str,
                eta_str,
                will_finish_str
            )

        except Exception:
            continue

    return table


if __name__ == "__main__":
    try:
        with Live(parse_md_logs(), refresh_per_second=0.2, transient=False) as live:
            while True:
                time.sleep(5)
                live.update(parse_md_logs())
    except KeyboardInterrupt:
        pass