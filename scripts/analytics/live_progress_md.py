import glob
import re
import os
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
    table.add_column("Start Time", style="dim white")
    table.add_column("Sim Time", style="green")
    table.add_column("Progress", style="yellow")
    table.add_column("Performance", style="magenta")
    #table.add_column("ETA", style="bold bold_yellow" if True else "white")
    table.add_column("ETA", style="bold yellow")
    log_files = glob.glob("results/gromacs/*/*/*/standard_100ns/JOB/*_md.log") + \
                glob.glob("logs/*/*/*/*_md.log")
    
    for log in sorted(log_files):
        try:
            with open(log, "r", encoding="utf-8", errors="ignore") as f:
                lines = f.readlines()
            
            if not lines:
                continue

            # 1. Determine Start Time from GROMACS header log or file creation
            start_time = None
            for line in lines[:50]:  # GROMACS prints log start header near the top
                match = re.search(r"Log file opened on (.*)", line)
                if match:
                    try:
                        start_time = datetime.strptime(match.group(1).strip(), "%a %b %d %H:%M:%S %Y")
                    except ValueError:
                        pass
                    break
            
            # Fallback to file creation time if string parsing missed
            if not start_time:
                start_time = datetime.fromtimestamp(os.path.getctime(log))

            start_str = start_time.strftime("%b %d %H:%M")

            # 2. Collect Checkpoint Timestamps for Step Rate Calculation
            checkpoints = []
            for line in reversed(lines[-200:]): # Look through recent buffer
                match = re.search(r"Writing checkpoint, step (\d+) at (.*)", line)
                if match:
                    step_val = int(match.group(1))
                    time_str = match.group(2).strip()
                    try:
                        parsed_time = datetime.strptime(time_str, "%a %b %d %H:%M:%S %Y")
                        checkpoints.append((step_val, parsed_time))
                    except ValueError:
                        continue
                if len(checkpoints) >= 2:  # Need at least two points to compute performance
                    break

            if not checkpoints:
                # No checkpoints written yet
                parts = log.split("/")
                idx = parts.index("results") if "results" in parts else 0
                name = f"{parts[idx+2]}/{parts[idx+3]}/{parts[idx+4]}" if "results" in parts else log
                table.add_row(name, start_str, "0.0 / 100 ns", "0.0%", "Starting...", "N/A")
                continue

            # Latest step data
            latest_step, latest_time = checkpoints[0]
            sim_ps = latest_step * DT_PS
            sim_ns = sim_ps / 1000.0
            pct = (latest_step / TOTAL_STEPS) * 100.0

            # 3. Calculate Performance (ns/day) and ETA
            if len(checkpoints) >= 2:
                prev_step, prev_time = checkpoints[1]
                step_delta = latest_step - prev_step
                time_delta_sec = (latest_time - prev_time).total_seconds()
            else:
                # Fallback to total run time if only 1 checkpoint exists
                step_delta = latest_step
                time_delta_sec = (latest_time - start_time).total_seconds()

            if time_delta_sec > 0 and step_delta > 0:
                ns_per_sec = (step_delta * DT_PS / 1000.0) / time_delta_sec
                ns_per_day = ns_per_sec * 86400.0
                
                remaining_ns = TOTAL_NS - sim_ns
                remaining_sec = remaining_ns / ns_per_sec if ns_per_sec > 0 else 0
                
                eta_dt = datetime.now() + timedelta(seconds=remaining_sec)
                
                # Format ETA: Show time if today, else include date
                if eta_dt.date() == datetime.now().date():
                    eta_str = eta_dt.strftime("Today at %H:%M")
                else:
                    eta_str = eta_dt.strftime("%b %d %H:%M")
                
                perf_str = f"{ns_per_day:.1f} ns/day"
            else:
                perf_str = "Calculating..."
                eta_str = "N/A"

            # Parse path identifier
            parts = log.split("/")
            if "results" in parts:
                idx = parts.index("results")
                name = f"{parts[idx+2]}/{parts[idx+3]}/{parts[idx+4]}"
            else:
                name = log

            table.add_row(
                name,
                start_str,
                f"{sim_ns:.2f} / {TOTAL_NS:.0f} ns",
                f"{pct:.1f}%",
                perf_str,
                eta_str
            )
        except Exception:
            continue
            
    return table

if __name__ == "__main__":
    with Live(parse_md_logs(), refresh_per_second=1) as live:
        while True:
            time.sleep(3)
            live.update(parse_md_logs())