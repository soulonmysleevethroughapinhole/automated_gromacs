import glob
import os
import re
import time
from datetime import datetime, timedelta
from rich.live import Live
from rich.table import Table

TOTAL_NS = 100.0  # Target length in ns
DT_PS = 0.002     # 2 fs timestep = 0.002 ps
TOTAL_STEPS = int((TOTAL_NS * 1000) / DT_PS)  # 50,000,000 steps

def parse_md_logs():
    table = Table(title="Live GROMACS MD Simulation Tracker (100 ns Runs)")
    table.add_column("System / Replicate", style="cyan")
    table.add_column("First Started", style="bold blue")
    table.add_column("Last Active", style="dim white")
    table.add_column("Sim Time", style="green")
    table.add_column("Progress", style="yellow")
    table.add_column("Performance", style="magenta")
    table.add_column("ETA", style="bold yellow")

    step_logs = sorted(glob.glob("results/gromacs/*/*/*/standard_100ns/simulation_steps.log"))

    if not step_logs:
        table.add_row("No simulations found", "-", "-", "-", "-", "-", "-")
        return table

    for log_path in step_logs:
        try:
            with open(log_path, "r", encoding="utf-8", errors="ignore") as f:
                lines = [line.strip() for line in f if line.strip()]

            if not lines:
                continue

            parts = log_path.split("/")
            if "results" in parts:
                idx = parts.index("results")
                name = f"{parts[idx+2]}/{parts[idx+3]}/{parts[idx+4]}"
            else:
                name = log_path

            # 1. Start & Last Active Timestamps
            first_line = lines[0]
            first_match = re.search(r"^\[(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2})\]", first_line)
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

            # 2. Latest Step & Details
            step_match = re.search(r"\[STEP:\s*([^\]]+)\]\s*(.*)", latest_line)
            step_code = step_match.group(1) if step_match else "UNKNOWN"
            step_details = step_match.group(2) if step_match else ""

            # 3. Direct Progress Extraction from HEARTBEAT lines
            hb_ns = None
            hb_dt = None
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

            # 4. Checkpoints Fallback
            sim_dir = os.path.dirname(log_path)
            gmx_logs = [
                f for f in glob.glob(os.path.join(sim_dir, "*.log")) + glob.glob(os.path.join(sim_dir, "JOB", "*.log"))
                if os.path.basename(f) != "simulation_steps.log"
            ]

            checkpoints = []
            if gmx_logs:
                for glog in gmx_logs:
                    try:
                        with open(glog, "r", encoding="utf-8", errors="ignore") as gf:
                            glines = gf.readlines()
                        for line in reversed(glines[-200:]):
                            match = re.search(r"Writing checkpoint, step (\d+) at (.*)", line)
                            if match:
                                step_val = int(match.group(1))
                                time_str = match.group(2).strip()
                                try:
                                    parsed_time = datetime.strptime(time_str, "%a %b %d %H:%M:%S %Y")
                                    checkpoints.append((step_val, parsed_time))
                                except ValueError:
                                    continue
                            if len(checkpoints) >= 2:
                                break
                    except Exception:
                        pass
                    if checkpoints:
                        break

            # 5. Determine State & Progress Formatting
            if step_code in ["DONE", "SIMULATION_SUMMARY"]:
                sim_time_str = "100.0 / 100 ns"
                progress_str = "100.0%"
                perf_str = "Finished"
                eta_str = "Done"

            elif hb_ns is not None:
                pct = min(100.0, (hb_ns / TOTAL_NS) * 100.0)
                sim_time_str = f"{hb_ns:.2f} / {TOTAL_NS:.0f} ns"
                progress_str = f"{pct:.1f}%"

                # Calculate speed using time elapsed since start or heartbeat interval
                elapsed_sec = (last_active_dt - first_start_dt).total_seconds()
                if elapsed_sec > 60 and hb_ns > 0:
                    ns_per_day = (hb_ns / elapsed_sec) * 86400.0
                    remaining_sec = ((TOTAL_NS - hb_ns) / hb_ns) * elapsed_sec
                    eta_dt = datetime.now() + timedelta(seconds=remaining_sec)
                    
                    if eta_dt.date() == datetime.now().date():
                        eta_str = eta_dt.strftime("Today %H:%M")
                    else:
                        eta_str = eta_dt.strftime("%b %d %H:%M")
                    perf_str = f"{ns_per_day:.1f} ns/day"
                else:
                    perf_str = "Calculating..."
                    eta_str = "Running"

            elif checkpoints:
                latest_step, latest_time = checkpoints[0]
                sim_ps = latest_step * DT_PS
                sim_ns = sim_ps / 1000.0
                pct = min(100.0, (latest_step / TOTAL_STEPS) * 100.0)

                if len(checkpoints) >= 2:
                    prev_step, prev_time = checkpoints[1]
                    step_delta = latest_step - prev_step
                    time_delta_sec = (latest_time - prev_time).total_seconds()
                else:
                    step_delta = latest_step
                    time_delta_sec = (latest_time - last_active_dt).total_seconds()

                if time_delta_sec > 0 and step_delta > 0:
                    ns_per_sec = (step_delta * DT_PS / 1000.0) / time_delta_sec
                    ns_per_day = ns_per_sec * 86400.0
                    remaining_ns = max(0.0, TOTAL_NS - sim_ns)
                    remaining_sec = remaining_ns / ns_per_sec if ns_per_sec > 0 else 0
                    
                    eta_dt = datetime.now() + timedelta(seconds=remaining_sec)
                    if eta_dt.date() == datetime.now().date():
                        eta_str = eta_dt.strftime("Today %H:%M")
                    else:
                        eta_str = eta_dt.strftime("%b %d %H:%M")
                    perf_str = f"{ns_per_day:.1f} ns/day"
                else:
                    perf_str = "Calculating..."
                    eta_str = "N/A"

                sim_time_str = f"{sim_ns:.2f} / {TOTAL_NS:.0f} ns"
                progress_str = f"{pct:.1f}%"

            else:
                sim_time_str = "0.0 / 100 ns"
                progress_str = "0.0%"

                if "SLURM" in step_code or "EXEC_MODE" in step_code or step_code == "HEARTBEAT":
                    job_id_match = re.search(r"Job\s*(\d+)", step_details, re.IGNORECASE)
                    job_str = f"Job {job_id_match.group(1)}" if job_id_match else "Queued"
                    perf_str = f"HPC ({job_str})"
                    eta_str = "Running on HPC"
                elif step_code == "RETRIEVE_START":
                    perf_str = "Syncing"
                    eta_str = "Downloading..."
                elif step_code == "ERROR":
                    perf_str = "FAILED"
                    eta_str = "Check Log"
                else:
                    perf_str = "Initializing"
                    eta_str = "Starting..."

            table.add_row(
                name,
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
    with Live(parse_md_logs(), refresh_per_second=0.1) as live:
        while True:
            time.sleep(10)
            live.update(parse_md_logs())