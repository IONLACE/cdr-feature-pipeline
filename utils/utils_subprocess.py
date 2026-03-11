from subprocess import run, CalledProcessError
from pathlib import Path
from typing import List


def run_command(cmd: List[str], cwd: Path | str | None = None):
    """
    Safely run external command and stop pipeline if it fails.

    Parameters
    ----------
    cmd : list[str]
        Command and arguments (no shell=True)
    cwd : Path | str | None
        Working directory to run command in.

    Returns
    -------
    subprocess.CompletedProcess
    """

    print("\n▶ Running command:")
    print(" ".join(cmd))

    try:
        result = run(
            cmd,
            check=True,
            capture_output=True,
            text=True,
            cwd=str(cwd) if cwd is not None else None,
        )

        print("✔ Success")
        return result

    except CalledProcessError as e:
        print("\nCOMMAND FAILED")
        print("Command:", " ".join(cmd))
        print("\nSTDOUT:\n", e.stdout)
        print("\nSTDERR:\n", e.stderr)

        # Stop entire pipeline immediately
        raise SystemExit("Pipeline stopped due to external tool failure.")
