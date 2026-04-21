"""Smoke tests ensuring the new CLIs expose a working --help surface."""

import subprocess
import sys


def _run_help(module: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, "-m", f"ExposoGraph.{module}", "--help"],
        capture_output=True,
        text=True,
        timeout=30,
    )


def test_flux_cli_help_exits_clean():
    result = _run_help("flux_cli")
    assert result.returncode == 0
    assert "usage" in result.stdout.lower() or "usage" in result.stderr.lower()


def test_exposure_cli_help_exits_clean():
    result = _run_help("exposure_cli")
    assert result.returncode == 0


def test_interaction_cli_help_exits_clean():
    result = _run_help("interaction_cli")
    assert result.returncode == 0
