"""PLINK2 subprocess wrappers and pfile counting utilities."""

import os
import subprocess
from pathlib import Path
from typing import List, Optional, Union

from config import PLINK


def run_plink2_command(
    args: List[str],
    check: bool = True,
    capture_output: bool = False,
    verbose: bool = True,
    memory_mb: Optional[int] = None,
) -> subprocess.CompletedProcess:
    """Execute a PLINK2 command."""
    cmd = [PLINK.EXECUTABLE, "--threads", str(PLINK.THREADS)]
    if memory_mb is not None:
        cmd += ["--memory", str(memory_mb)]
    cmd += args

    if verbose:
        print(f"Running: {' '.join(cmd)}")

    return subprocess.run(
        cmd,
        check=check,
        capture_output=capture_output,
        text=True if capture_output else None,
    )


def count_variants(pfile_prefix: Union[str, Path]) -> int:
    """Count variants in a PLINK2 pfile (.pvar)."""
    pvar_path = f"{pfile_prefix}.pvar"
    if os.path.exists(pvar_path):
        with open(pvar_path) as f:
            return sum(1 for line in f if not line.startswith("#"))
    return -1


def count_samples(pfile_prefix: Union[str, Path]) -> int:
    """Count samples in a PLINK2 pfile (.psam), excluding the header."""
    psam_path = f"{pfile_prefix}.psam"
    if os.path.exists(psam_path):
        with open(psam_path) as f:
            return sum(1 for line in f if not line.startswith("#")) - 1
    return -1
