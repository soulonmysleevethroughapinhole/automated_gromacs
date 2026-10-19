"""Shared pytest fixtures for HPC orchestration tests.

Strategy: hpc_engine.py only talks to the cluster through subprocess.run /
subprocess.check_output. We monkeypatch those with a FakeRunner that records
every invocation and replays canned responses keyed by substring of the
command line. This lets us assert on the exact remote commands issued
(sbatch/squeue/rsync paths, quoting, retry behaviour) without any network.
"""

import logging
import sys
import types
from pathlib import Path

import pytest

# Make `import yaml` safe even if PyYAML is absent in a bare checkout.
if "yaml" not in sys.modules:
    try:
        import yaml  # noqa: F401
    except ImportError:  # pragma: no cover
        stub = types.ModuleType("yaml")
        stub.safe_load = lambda *a, **k: {}
        stub.safe_dump = lambda *a, **k: None
        sys.modules["yaml"] = stub

import hpc_engine


class CompletedProcess:
    """Minimal stand-in for subprocess.CompletedProcess."""

    def __init__(self, returncode=0, stdout="", stderr=""):
        self.returncode = returncode
        self.stdout = stdout
        self.stderr = stderr


class FakeRunner:
    """Records subprocess calls and replays canned responses.

    Responses are {substring: response} pairs; the LONGEST matching substring
    wins (deterministic regardless of dict insertion order, so specific keys
    like "squeue -j" beat generic ones like "-j"). A response is either a
    CompletedProcess, a callable(command_text) -> CompletedProcess, or an
    Exception to raise. Unmatched commands return an empty success (so
    incidental mkdir etc. don't fail the test).
    """

    def __init__(self, responses=None):
        self.calls = []           # list of command (list or str)
        self.responses = responses or {}
        self._seq_keys = set()    # keys whose value is a sequence (consumed in order)

    def add_sequence(self, key, items):
        """Register an ordered list of responses for `key`.

        Each matching call consumes the next item; once exhausted, every
        further call reuses the LAST item. This mirrors how the real cluster
        behaves: transient states (queue polls) advance, terminal states stick.
        """
        self.responses[key] = list(items)
        self._seq_keys.add(key)

    def _pick(self, key):
        resp = self.responses[key]
        if key in self._seq_keys:
            if len(resp) > 1:
                return resp.pop(0)
            return resp[0]
        return resp

    def _normalize(self, cmd):
        if isinstance(cmd, list):
            return " ".join(cmd)
        # collapse whitespace in shell strings so keys like "sbatch" or
        # "_123.err" match regardless of indentation/line continuations
        return " ".join(str(cmd).split())

    def _respond(self, cmd):
        text = self._normalize(cmd)
        for key in self.responses:
            if key in text:
                resp = self._pick(key)
                if isinstance(resp, Exception):
                    raise resp
                if callable(resp):
                    return resp(text)
                return resp
        return CompletedProcess(0, "", "")

    def _match_key(self, text):
        """Longest matching response key wins (deterministic priority)."""
        matches = [k for k in self.responses if k in text]
        if not matches:
            return None
        return max(matches, key=len)

    def run(self, cmd, *args, **kwargs):
        self.calls.append(cmd)
        key = self._match_key(self._normalize(cmd))
        if key is None:
            return CompletedProcess(0, "", "")
        resp = self._pick(key)
        if isinstance(resp, Exception):
            raise resp
        if callable(resp):
            return resp(self._normalize(cmd))
        return resp

    def check_output(self, cmd, *args, **kwargs):
        res = self.run(cmd, *args, **kwargs)
        if res.returncode != 0:
            raise RuntimeError(f"check_output failed: {cmd}")
        return res.stdout.encode()


@pytest.fixture
def fake_runner(monkeypatch):
    """Patch subprocess inside hpc_engine; returns the FakeRunner."""
    runner = FakeRunner()
    monkeypatch.setattr(hpc_engine.subprocess, "run", runner.run)
    monkeypatch.setattr(hpc_engine.subprocess, "check_output", runner.check_output)
    return runner


@pytest.fixture
def fast_sleep(monkeypatch):
    """Make time.sleep instant and record how long each call *would* wait."""
    slept = []
    monkeypatch.setattr(hpc_engine.time, "sleep", lambda s: slept.append(s))
    return slept


@pytest.fixture
def logger():
    log = logging.getLogger("test-hpc")
    log.setLevel(logging.CRITICAL)
    log.handlers = [logging.NullHandler()]
    return log


@pytest.fixture
def manager(tmp_path, logger):
    """A ready-to-use HPCJobManager rooted in a temp directory.

    NOTE: hpc_engine.monitor_job appends to <local_dir>/job_execution_status.log
    without creating local_dir first. In the real workflow Snakemake guarantees
    the output dir exists; here we mirror that precondition with mkdir. (The
    refactor should make the engine create it itself -- see tech-debt note.)
    """
    (tmp_path / "local").mkdir()
    return hpc_engine.HPCJobManager(
        ssh_target="testuser@test.cluster.edu",
        remote_dir="/scratch/testuser/md_runs/1abc/A/1/standard_100ns",
        local_dir=str(tmp_path / "local"),
        job_name="md_1abc_A_1_standard_100ns",
        job_script="md_job.job",
        log_path=str(tmp_path / "logs" / "md.log"),
        logger=logger,
        hpc_user="testuser",
    )
