"""
Tests for maws.tools module.

This module tests the executable utility functions:
- find_exe(): Locate executables on PATH
- run(): Execute commands via subprocess
- ExecError: Custom exception for missing executables
"""

import subprocess

import pytest

from maws.tools import ExecError, find_exe, run


class TestFindExe:
    """Tests for find_exe() function."""

    def test_find_exe_finds_common_command(self):
        """find_exe() finds common commands like 'echo' or 'python'."""
        # 'echo' exists on all Unix systems
        path = find_exe("echo")
        assert path is not None
        assert "echo" in path

    def test_find_exe_raises_for_missing(self):
        """find_exe() raises ExecError for non-existent executables."""
        with pytest.raises(ExecError) as exc_info:
            find_exe("this_command_does_not_exist_xyz123")
        assert "not found on PATH" in str(exc_info.value)

    def test_exec_error_is_runtime_error(self):
        """ExecError is a RuntimeError subclass."""
        assert issubclass(ExecError, RuntimeError)


class TestRun:
    """Tests for run() function."""

    def test_run_executes_simple_command(self):
        """run() executes a simple command successfully."""
        # This should not raise
        run(["echo", "hello"])

    def test_run_raises_on_failure(self):
        """run() raises CalledProcessError on non-zero exit."""
        with pytest.raises(subprocess.CalledProcessError):
            run(["false"])  # 'false' always exits with 1

    def test_run_with_cwd(self, tmp_path):
        """run() respects the cwd argument."""
        # Create a file in temp directory
        test_file = tmp_path / "test.txt"
        test_file.write_text("hello")

        # Run 'ls' in the temp directory
        run(["ls"], cwd=tmp_path)
        # If it didn't raise, the cwd was valid
