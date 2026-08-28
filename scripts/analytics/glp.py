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


def read_tail(file_path: str, max_bytes: int = 32768) -> list[str]:
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
    """Formats estimated time remaining into clean clock or calendar targets."""
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


def parse_md_logs() -> Table:
    table = Table(title="Live GROMACS MD Simulation Tracker (100 ns Runs)", expand=True)
    table.add_column("System / Replicate", style="cyan", no_wrap=True)
    table.add_column("Target", style="bold white", no_wrap=True)
    table.add_column("First Started", style="blue")
    table.add_column("Last Active", style="dim white")
    table.add_column("Sim Time", style="green")
    table.add_column("Progress", style="yellow")
    table.add_column("Performance", style="magenta")
    table.add_column("ETA", style="bold yellow")

    # Target path matching: results/gromacs/<PDB>/<MODEL>/<REP>/standard_100ns/JOB/
    step_logs = sorted(glob.glob("results/gromacs/*/*/*/standard_100ns/JOB/simulation_steps.log"))

    if not step_logs:
        table.add_row("No simulations found", "-", "-", "-", "-", "-", "-", "-")
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
                name = "/".join(parts[idx+1:idx+4])
            else:
                name = log_path

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

            # 4. Step Code & Status Context
            step_match = re.search(r"\[STEP:\s*([^\]]+)\]\s*(.*)", latest_line)
            step_code = step_match.group(1) if step_match else "UNKNOWN"
            step_details = step_match.group(2) if step_match else ""

            # Check if current state explicitly indicates queued on HPC
            is_queued = "State: PD" in step_details or (step_code == "HEARTBEAT" and "PD" in step_details)

            # 5. GROMACS Log Direct Extraction (Production Only)
            sim_dir = os.path.dirname(log_path)
            
            # Target ONLY production MD logs (*_md.log or logs within JOB/)
            gmx_logs = [
                f for f in glob.glob(os.path.join(sim_dir, "*_md.log")) + glob.glob(os.path.join(sim_dir, "*.log"))
                if os.path.basename(f) != "simulation_steps.log" and not os.path.basename(f).startswith(("cg", "st"))
            ]

            gmx_ns = None
            checkpoints = []

            # Skip production log parsing if the job is explicitly pending in Slurm
            if gmx_logs and not is_queued:
                for glog in gmx_logs:
                    glines = read_tail(glog, max_bytes=65536)
                    for idx, line in enumerate(reversed(glines)):
                        # Match checkpoint write
                        cp_match = re.search(r"Writing checkpoint, step (\d+) at (.*)", line)
                        if cp_match:
                            step_val = int(cp_match.group(1))
                            time_str = cp_match.group(2).strip()
                            try:
                                parsed_time = datetime.strptime(time_str, "%a %b %d %H:%M:%S %Y")
                                checkpoints.append((step_val, parsed_time))
                            except ValueError:
                                pass

                        # Match step/time block
                        if gmx_ns is None and "Step" in line and "Time" in line:
                            actual_idx = len(glines) - 1 - idx
                            if actual_idx + 1 < len(glines):
                                val_line = glines[actual_idx + 1].strip()
                                parts_val = val_line.split()
                                if len(parts_val) >= 2 and parts_val[0].isdigit():
                                    try:
                                        time_ps = float(parts_val[1])
                                        gmx_ns = time_ps / 1000.0
                                    except ValueError:
                                        pass

                        if len(checkpoints) >= 2 and gmx_ns is not None:
                            break
                    if checkpoints or gmx_ns is not None:
                        break

            # 6. Extract Wrapper Heartbeats (Fallback)
            hb_ns = None
            hb_dt = None
            if not is_queued:
                for line in reversed(lines):
                    hb_match = re.search(
                        r"^\[(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2})\].*?running\s*-\s*([\d\.]+)\s*ns",
                        line,
                        re.IGNORECASE
                    )
                    if hb_match:
                        hb_dt = datetime.strptime(hb_match.group(1), "%Y-%m-%d %H:%M:%S")
                        hb_ns = float(hb_match.group(2))
                        break

            # 7. Priority Metrics Evaluation
            current_ns = gmx_ns if gmx_ns is not None else hb_ns

            if is_queued:
                job_match = re.search(r"Job\s*(\d+)", step_details, re.IGNORECASE)
                job_str = f"Job {job_match.group(1)}" if job_match else "Queued"
                sim_time_str = f"0.0 / {TOTAL_NS:.0f} ns"
                progress_str = "0.0%"
                perf_str = f"Queued ({job_str})"
                eta_str = "[yellow]Pending (PD)[/yellow]"

            elif step_code in ["DONE", "SIMULATION_SUMMARY"] or (current_ns and current_ns >= TOTAL_NS):
                sim_time_str = f"{TOTAL_NS:.1f} / {TOTAL_NS:.0f} ns"
                progress_str = "100.0%"
                perf_str = "[green]Finished[/green]"
                eta_str = "[bold green]Done[/bold green]"

            elif len(checkpoints) >= 2:
                latest_step, latest_time = checkpoints[0]
                prev_step, prev_time = checkpoints[1]
                
                sim_ns = current_ns if current_ns is not None else (latest_step * DT_PS / 1000.0)
                pct = min(100.0, (sim_ns / TOTAL_NS) * 100.0)

                step_delta = latest_step - prev_step
                delta_t = (latest_time - prev_time).total_seconds()

                if delta_t > 0 and step_delta > 0 and not is_stale:
                    ns_per_sec = (step_delta * DT_PS / 1000.0) / delta_t
                    ns_per_day = ns_per_sec * 86400.0
                    remaining_sec = (TOTAL_NS - sim_ns) / ns_per_sec if ns_per_sec > 0 else 0
                    
                    perf_str = f"{ns_per_day:.1f} ns/day"
                    eta_str = format_eta(remaining_sec)
                elif is_stale:
                    perf_str = "[bold yellow]STALLED[/bold yellow]"
                    eta_str = "Inactive"
                else:
                    perf_str = "Calculating..."
                    eta_str = "N/A"

                sim_time_str = f"{sim_ns:.2f} / {TOTAL_NS:.0f} ns"
                progress_str = f"{pct:.1f}%"

            elif current_ns is not None:
                pct = min(100.0, (current_ns / TOTAL_NS) * 100.0)
                sim_time_str = f"{current_ns:.2f} / {TOTAL_NS:.0f} ns"
                progress_str = f"{pct:.1f}%"

                ref_dt = hb_dt if hb_dt else last_active_dt
                elapsed_sec = (ref_dt - first_start_dt).total_seconds()

                if elapsed_sec > 60 and current_ns > 0 and not is_stale:
                    ns_per_day = (current_ns / elapsed_sec) * 86400.0
                    remaining_sec = ((TOTAL_NS - current_ns) / current_ns) * elapsed_sec
                    perf_str = f"{ns_per_day:.1f} ns/day"
                    eta_str = format_eta(remaining_sec)
                elif is_stale:
                    perf_str = "[bold yellow]STALLED[/bold yellow]"
                    eta_str = "Inactive"
                else:
                    perf_str = "Calculating..."
                    eta_str = "Running"

            else:
                sim_time_str = f"0.0 / {TOTAL_NS:.0f} ns"
                progress_str = "0.0%"

                if "SLURM" in step_code or "EXEC_MODE" in step_code or step_code == "HEARTBEAT":
                    job_match = re.search(r"Job\s*(\d+)", step_details, re.IGNORECASE)
                    job_str = f"Job {job_match.group(1)}" if job_match else "Queued"
                    perf_str = f"HPC ({job_str})"
                    eta_str = "Queued"
                elif step_code == "RETRIEVE_START":
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
                eta_str
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