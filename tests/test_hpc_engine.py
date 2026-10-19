"""Characterization tests for scripts/smk/sims/hpc_engine.py.

These tests pin the *current* observable behaviour of the HPC orchestration
layer so the upcoming refactor of sims.smk / hpc_engine is provably safe:
if a test breaks after refactoring, either the refactor changed semantics
(bad) or we consciously accepted a behaviour change (document it here).

Covered surface:
  - command construction for ssh / squeue / sbatch / rsync
  - job-id parsing from sbatch and squeue output
  - Slurm submit-limit retry loop (backoff + eventual failure)
  - Lmod/Cray stderr noise filtering
  - monitor_job state machine: complete / checkpoint-resubmit / dead job
  - execute_pipeline idempotency (local sentinel, remote sentinel, resume)
"""

import os
import subprocess

import pytest

import hpc_engine
from conftest import CompletedProcess


def flat(cmd):
    """Flatten a recorded subprocess call (argv list or shell string) to text."""
    if isinstance(cmd, list):
        return " ".join(cmd)
    return " ".join(str(cmd).split())


# ---------------------------------------------------------------------------
# Construction / plumbing
# ---------------------------------------------------------------------------

class TestConstruction:
    def test_paths_are_absolute(self, tmp_path, logger):
        mgr = hpc_engine.HPCJobManager(
            ssh_target="host", remote_dir="/remote", local_dir="relative/local",
            job_name="j", job_script="job.sh", log_path="rel.log", logger=logger,
        )
        assert os.path.isabs(mgr.local_dir)
        assert os.path.isabs(mgr.log_path)

    def test_hpc_user_defaults_to_env(self, tmp_path, logger, monkeypatch):
        monkeypatch.setenv("USER", "alice")
        mgr = hpc_engine.HPCJobManager(
            ssh_target="host", remote_dir="/r", local_dir=str(tmp_path),
            job_name="j", job_script="job.sh", log_path=str(tmp_path / "l.log"),
            logger=logger,
        )
        assert mgr.hpc_user == "alice"

    def test_status_file_lives_in_local_dir(self, manager):
        assert manager.local_status_file == os.path.join(
            manager.local_dir, "job_execution_status.log")


# ---------------------------------------------------------------------------
# _run_ssh: exact argv contract
# ---------------------------------------------------------------------------

class TestRunSsh:
    def test_uses_batchmode_and_list_argv(self, manager, fake_runner):
        manager._run_ssh("echo hi")
        assert fake_runner.calls == [
            ["ssh", "-o", "BatchMode=yes", manager.ssh_target, "echo hi"]
        ]


# ---------------------------------------------------------------------------
# Remote state probes
# ---------------------------------------------------------------------------

class TestIsRemoteFinished:
    def test_yes_when_sentinel_present(self, manager, fake_runner):
        fake_runner.responses["[ -f"] = CompletedProcess(0, "YES\n", "")
        assert manager.is_remote_finished("md.gro") is True
        # sentinel path must be joined onto remote_dir
        cmd = flat(fake_runner.calls[-1])
        assert f"{manager.remote_dir}/md.gro" in cmd

    def test_no_when_absent(self, manager, fake_runner):
        fake_runner.responses["[ -f"] = CompletedProcess(0, "NO\n", "")
        assert manager.is_remote_finished("md.gro") is False

    def test_no_on_ssh_failure(self, manager, fake_runner):
        fake_runner.responses["[ -f"] = CompletedProcess(255, "", "ssh: connect failed")
        assert manager.is_remote_finished("md.gro") is False


class TestGetActiveJobId:
    def test_parses_id_and_ignores_state(self, manager, fake_runner):
        fake_runner.responses["squeue"] = CompletedProcess(0, "123456 R\n", "")
        assert manager.get_active_job_id() == "123456"
        cmd = flat(fake_runner.calls[-1])
        assert "-u testuser" in cmd and "-n md_1abc_A_1_standard_100ns" in cmd

    def test_none_when_queue_empty(self, manager, fake_runner):
        fake_runner.responses["squeue"] = CompletedProcess(0, "\n", "")
        assert manager.get_active_job_id() is None


class TestHasCheckpoint:
    @pytest.mark.parametrize("stdout,expected", [
        ("0\n", False),
        ("2\n", True),
        ("garbage", False),   # ValueError -> False (fail-safe)
    ])
    def test_counts_cpt_files(self, manager, fake_runner, stdout, expected):
        fake_runner.responses[".cpt"] = CompletedProcess(0, stdout, "")
        assert manager.has_checkpoint() is expected


# ---------------------------------------------------------------------------
# submit_job
# ---------------------------------------------------------------------------

class TestSubmitJob:
    def test_returns_job_id_from_sbatch_output(self, manager, fake_runner):
        fake_runner.responses["sbatch"] = CompletedProcess(
            0, "Submitted batch job 987654\n", "")
        assert manager.submit_job() == "987654"

    def test_unpack_flag_controls_tar_extraction(self, manager, fake_runner):
        fake_runner.responses["sbatch"] = CompletedProcess(0, "Submitted batch job 1\n", "")
        manager.submit_job(unpack_job_archive=True)
        assert "tar -xzf JOB.tar.gz" in flat(fake_runner.calls[-1])
        manager.submit_job(unpack_job_archive=False)
        assert "tar -xzf" not in flat(fake_runner.calls[-1])

    def test_missing_script_raises(self, manager, fake_runner):
        fake_runner.responses["sbatch"] = CompletedProcess(
            1, "", "Error: Explicit script 'md_job.job' not found in /scratch/x")
        with pytest.raises(RuntimeError, match="sbatch submission failed"):
            manager.submit_job()

    def test_queue_limit_retries_then_succeeds(self, manager, fake_runner, fast_sleep):
        calls = {"n": 0}

        def flaky(text):
            calls["n"] += 1
            if calls["n"] < 3:
                return CompletedProcess(1, "", "slurm_load_jobs error: "
                                        "AssocGrpSubmitJobsLimit reached")
            return CompletedProcess(0, "Submitted batch job 555\n", "")

        fake_runner.responses["sbatch"] = flaky
        assert manager.submit_job() == "555"
        assert fast_sleep == [300, 300]          # two backoffs before success

    def test_queue_limit_exhausts_retries(self, manager, fake_runner, fast_sleep):
        fake_runner.responses["sbatch"] = CompletedProcess(
            1, "", "AssocGrpSubmitJobsLimit reached")
        with pytest.raises(RuntimeError, match="after 6 attempts"):
            manager.submit_job()
        assert len(fast_sleep) == 6              # slept on every attempt...
        assert len(fake_runner.calls) == 6       # ...and tried exactly 6 times

    def test_lmod_noise_is_filtered(self, manager, fake_runner):
        noisy = (
            "Lmod is automatically replaced...\n"
            "Due to MODULEPATH changes, the following have been reloaded:\n"
            "  1) gromacs/2024     2) cuda/12.2\n"
            "real problem: sbatch: error: Batch job submission failed\n"
        )
        fake_runner.responses["sbatch"] = CompletedProcess(1, "", noisy)
        with pytest.raises(RuntimeError) as exc:
            manager.submit_job()
        msg = str(exc.value)
        assert "real problem" in msg
        assert "Lmod" not in msg and "MODULEPATH" not in msg


# ---------------------------------------------------------------------------
# pull_results
# ---------------------------------------------------------------------------

class TestPullResults:
    def test_rsync_contract(self, manager, fake_runner):
        manager.pull_results(target_subdir="md_results_HPC")
        cmd = flat(fake_runner.calls[-1])
        assert "rsync" in cmd
        # Current implementation wires paths through bash variables rather than
        # interpolating them directly into the rsync argv. Pin both halves.
        assert f'SSH_TARGET="{manager.ssh_target}"' in cmd
        assert f'REMOTE_DIR="{manager.remote_dir}"' in cmd
        assert '"$SSH_TARGET:$REMOTE_DIR/"' in cmd
        local_target = os.path.join(manager.local_dir, "md_results_HPC")
        assert f'LOCAL_TARGET_DIR="{local_target}"' in cmd
        assert '--exclude="JOB"' in cmd and '--exclude="JOB.tar.gz"' in cmd

    def test_creates_local_target_dir(self, manager, fake_runner):
        target = os.path.join(manager.local_dir, "out")
        manager.pull_results(target_subdir="out")
        assert os.path.isdir(target)

    def test_relative_remote_path_gets_rooted(self, tmp_path, logger, fake_runner):
        mgr = hpc_engine.HPCJobManager(
            ssh_target="host", remote_dir="scratch/user/run", local_dir=str(tmp_path),
            job_name="j", job_script="job.sh", log_path=str(tmp_path / "l.log"),
            logger=logger,
        )
        mgr.pull_results()
        assert 'REMOTE_DIR="/scratch/user/run"' in flat(fake_runner.calls[-1])


# ---------------------------------------------------------------------------
# monitor_job state machine
# ---------------------------------------------------------------------------

def _squeue_state(responses_module, state):
    """Helper: canned squeue output for a given Slurm state letter."""
    return CompletedProcess(0, f"{state}\n", "")


class TestMonitorJob:
    def test_running_then_completed(self, manager, fake_runner, fast_sleep):
        fake_runner.add_sequence("-j", [
            _squeue_state(None, "R"),   # poll 1: running
            _squeue_state(None, ""),    # poll 2: left queue -> verify completion
        ])
        fake_runner.responses["[ -f"] = CompletedProcess(0, "YES", "")

        manager.monitor_job("123", "md.gro", poll_interval_sec=1)

        status = open(manager.local_status_file).read()
        assert "RUNNING" in status
        assert os.path.exists(manager.local_status_file)

    def test_resubmits_on_walltime_with_checkpoint(self, manager, fake_runner, fast_sleep):
        """Job vanishes without sentinel but .cpt exists -> resubmit, then finish."""
        fake_runner.add_sequence("squeue -j", [_squeue_state(None, "")])  # both jobs leave queue
        fake_runner.responses[".cpt"] = CompletedProcess(0, "1\n", "")  # checkpoint exists
        fake_runner.add_sequence("[ -f", [
            CompletedProcess(0, "NO\n", ""),    # not finished at first look
            CompletedProcess(0, "YES\n", ""),   # finished after resubmit
        ])
        fake_runner.responses["sbatch '"] = CompletedProcess(0, "Submitted batch job 222\n", "")

        manager.monitor_job("111", "md.gro", poll_interval_sec=1)

        submitted = [c for c in fake_runner.calls if "sbatch" in flat(c)]
        assert len(submitted) == 1             # exactly one auto-resubmission

    def test_gives_up_without_sentinel_or_checkpoint(self, manager, fake_runner, fast_sleep):
        fake_runner.responses["-j"] = CompletedProcess(0, "", "")
        fake_runner.responses[".cpt"] = CompletedProcess(0, "0\n", "")
        fake_runner.responses["[ -f"] = CompletedProcess(0, "NO\n", "")

        manager.monitor_job("999", "md.gro", poll_interval_sec=1)  # must terminate

    def test_queued_state_heartbeat(self, manager, fake_runner, fast_sleep):
        fake_runner.add_sequence("-j", [
            _squeue_state(None, "PD"),
            _squeue_state(None, ""),
        ])
        fake_runner.responses["[ -f"] = CompletedProcess(0, "YES", "")

        manager.monitor_job("5", "md.gro", poll_interval_sec=1)
        assert "QUEUED" in open(manager.local_status_file).read()


# ---------------------------------------------------------------------------
# execute_pipeline (top-level orchestration)
# ---------------------------------------------------------------------------

class TestExecutePipeline:
    def test_local_sentinel_short_circuits_everything(self, manager, fake_runner, tmp_path):
        sentinel = os.path.join(manager.local_dir, "md_results", "md_completed.txt")
        os.makedirs(os.path.dirname(sentinel))
        open(sentinel, "w").close()

        manager.execute_pipeline("md.gro")
        assert fake_runner.calls == []          # zero cluster traffic

    def test_already_finished_remotely_only_pulls(self, manager, fake_runner):
        fake_runner.responses["[ -f"] = CompletedProcess(0, "YES\n", "")
        manager.execute_pipeline("md.gro")
        text = " ".join(flat(c) for c in fake_runner.calls)
        assert "rsync" in text
        assert "sbatch" not in text

    def test_attaches_to_existing_active_job(self, manager, fake_runner, fast_sleep):
        fake_runner.add_sequence("[ -f", [
            CompletedProcess(0, "NO\n", ""),    # pre-check
            CompletedProcess(0, "YES\n", ""),   # post-monitor verification
        ])
        # squeue -u (active-job lookup) finds the running job; note the sbatch
        # command also contains "-n <jobname>", so keys must disambiguate.
        fake_runner.responses["-u testuser"] = CompletedProcess(0, "4242 R\n", "")
        fake_runner.add_sequence("squeue -j", [
            CompletedProcess(0, "R\n", ""),
            CompletedProcess(0, "", ""),
        ])

        manager.execute_pipeline("md.gro", poll_interval_sec=1)

        assert not any("sbatch" in flat(c) for c in fake_runner.calls)

    def test_full_submit_monitor_pull_cycle(self, manager, fake_runner, fast_sleep):
        fake_runner.add_sequence("[ -f", [
            CompletedProcess(0, "NO\n", ""),    # pre-check: not done
            CompletedProcess(0, "YES\n", ""),   # post-monitor verification
        ])
        fake_runner.responses["-u testuser"] = CompletedProcess(0, "", "")  # nothing queued
        fake_runner.responses["sbatch '"] = CompletedProcess(0, "Submitted batch job 77\n", "")
        fake_runner.add_sequence("squeue -j", [
            CompletedProcess(0, "R\n", ""),
            CompletedProcess(0, "", ""),
        ])

        manager.execute_pipeline("md.gro", poll_interval_sec=1)

        text = " ".join(flat(c) for c in fake_runner.calls)
        assert "mkdir -p" in text and "sbatch" in text and "rsync" in text

    def test_raises_if_never_completes(self, manager, fake_runner, fast_sleep):
        fake_runner.responses["[ -f"] = CompletedProcess(0, "NO\n", "")
        fake_runner.responses["-u testuser"] = CompletedProcess(0, "", "")  # nothing queued
        fake_runner.responses["sbatch '"] = CompletedProcess(0, "Submitted batch job 77\n", "")
        fake_runner.responses["squeue -j"] = CompletedProcess(0, "", "")    # job already gone
        fake_runner.responses[".cpt"] = CompletedProcess(0, "0\n", "")     # no checkpoint

        with pytest.raises(RuntimeError, match="without producing"):
            manager.execute_pipeline("md.gro", poll_interval_sec=1)


# ---------------------------------------------------------------------------
# get_simulation_progress
# ---------------------------------------------------------------------------

class TestGetSimulationProgress:
    def test_by_job_id_greps_err_file(self, fake_runner):
        fake_runner.responses[".err"] = CompletedProcess(0, "Step 100, Time 0.2\n", "")
        out = hpc_engine.get_simulation_progress("host", "/remote", job_id="123")
        assert out == "Step 100, Time 0.2"
        cmd = flat(fake_runner.calls[-1])
        assert "_123.err" in cmd

    def test_fallback_uses_log(self, fake_runner):
        fake_runner.responses[".log"] = CompletedProcess(0, "ns/day: 18.5\n", "")
        out = hpc_engine.get_simulation_progress("host", "/remote")
        assert out == "ns/day: 18.5"

    def test_none_on_empty_output(self, fake_runner):
        assert hpc_engine.get_simulation_progress("host", "/remote") is None
