# tools.py
import shutil, subprocess, tempfile
from pathlib import Path

class ExecError(RuntimeError): pass

def find_exe(name: str) -> str:
    exe = shutil.which(name)
    if not exe:
        raise ExecError(f"{name} not found on PATH. Install AmberTools and ensure it’s on PATH.")
    return exe

def run(cmd: list[str], cwd: str | Path | None = None) -> None:
    subprocess.run(cmd, cwd=cwd, check=True)
