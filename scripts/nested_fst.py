"""Per-fold, leak-free FST pool selection via plink2 (Option B).

The candidate-set leakage fix requires the FST pool to be nominated using ONLY
the training samples of each CV fold. Rather than reimplement Hudson FST, we call
the same `plink2 --fst` the rest of the pipeline uses (04b), restricted to the
fold's training samples with `--keep`. plink2 then estimates allele frequencies —
and therefore per-variant Hudson FST — from those samples alone.

Alignment guards (Option B is only safe if these hold):
  * the numpy genotype-matrix sample order equals the .psam #IID order
    (verified once; re-asserted here per call via the keep set), and
  * every FST variant ID maps to a numpy column (asserted; 100% on full data).

Returns FST-pool membership as *column indices into the numpy genotype matrix*,
so the caller can slice the same array the CV splits operate on.
"""

from __future__ import annotations

import glob
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence

import numpy as np
import pandas as pd


def _run_plink_fst(
    pfile: str,
    keep_file: str,
    out_prefix: str,
    pop_column: str,
    plink: str,
    threads: int,
) -> None:
    cmd = [
        plink,
        "--pfile", pfile,
        "--keep", keep_file,
        "--fst", pop_column, "report-variants",
        "--out", out_prefix,
    ]
    if threads and threads > 0:
        cmd += ["--threads", str(threads)]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            f"plink2 --fst failed (rc={proc.returncode})\n"
            f"CMD: {' '.join(cmd)}\nSTDERR:\n{proc.stderr[-2000:]}"
        )


def _assert_kept_sample_count(log_path: str, expected: int) -> None:
    """Defense-in-depth: confirm plink2 actually kept `expected` samples.

    Parses the .log for the sample count reported after --keep. Phrasing varies
    slightly across plink2 builds, so we look for the standard patterns and, if
    none match, raise rather than silently trust the keep file.
    """
    if not os.path.exists(log_path):
        raise RuntimeError(f"plink2 log not found for sample-count check: {log_path}")
    text = Path(log_path).read_text()
    patterns = [
        r"--keep:\s*([\d,]+)\s+samples? remaining",
        r"([\d,]+)\s+samples? remaining after",
    ]
    for pat in patterns:
        m = re.search(pat, text)
        if m:
            got = int(m.group(1).replace(",", ""))
            if got != expected:
                raise AssertionError(
                    f"plink2 kept {got} samples but fold expected {expected} "
                    f"(sample-alignment guard tripped)"
                )
            return
    # If --keep removed nothing, plink2 may not print a "remaining after --keep"
    # line. Fall back to the "loaded" count and require it to equal `expected`.
    m = re.search(r"([\d,]+)\s+samples?\b.*loaded", text)
    if m and int(m.group(1).replace(",", "")) == expected:
        return
    raise AssertionError(
        f"could not confirm plink2 kept exactly {expected} samples from log "
        f"{log_path} (Option B sample-alignment guard cannot be verified)"
    )


def fst_pool_for_samples(
    sample_ids: Sequence[str],
    *,
    pfile: str,
    snp_col_index: Dict[str, int],
    psam_ids: Optional[Sequence[str]] = None,
    top_n: int = 1000,
    pop_column: str = "pop",
    plink: str = "plink2",
    threads: int = 0,
    workdir: Optional[str] = None,
    cleanup: bool = True,
    fst_col: str = "HUDSON_FST",
    id_col: str = "ID",
    return_details: bool = False,
):
    """Nominate the FST pool from `sample_ids` only.

    Args:
        sample_ids: training-fold sample #IIDs (e.g. numpy_samples[train_idx]).
        pfile: LD-pruned pfile prefix; its .psam must carry `pop_column`.
        snp_col_index: {variant_id -> numpy column index}.
        psam_ids: optional full psam #IID list, for a membership assertion.
        top_n: variants per pairwise comparison (union across pairs).
    Returns:
        np.ndarray of sorted unique numpy column indices (the FST pool),
        or (indices, details_dict) if return_details.
    """
    sample_ids = [str(s) for s in sample_ids]
    if len(set(sample_ids)) != len(sample_ids):
        raise ValueError("sample_ids contains duplicates")
    if psam_ids is not None:
        missing = set(sample_ids) - set(map(str, psam_ids))
        if missing:
            raise ValueError(
                f"{len(missing)} fold sample(s) absent from .psam, e.g. "
                f"{list(missing)[:3]} — sample-alignment guard failed"
            )

    own_workdir = workdir is None
    workdir = workdir or tempfile.mkdtemp(prefix="foldfst_")
    os.makedirs(workdir, exist_ok=True)
    try:
        keep_file = os.path.join(workdir, "keep.txt")
        with open(keep_file, "w") as fh:
            fh.write("#IID\n")
            fh.write("\n".join(sample_ids) + "\n")

        out_prefix = os.path.join(workdir, "fold_fst")
        _run_plink_fst(pfile, keep_file, out_prefix, pop_column, plink, threads)
        _assert_kept_sample_count(out_prefix + ".log", len(sample_ids))

        fst_files = sorted(glob.glob(out_prefix + ".*.fst.var"))
        if not fst_files:
            raise RuntimeError(f"no .fst.var files produced at {out_prefix}.*")

        pool_ids: set = set()
        per_pair = {}
        for f in fst_files:
            df = pd.read_csv(f, sep=r"\s+", usecols=[id_col, fst_col])
            df = df[df[fst_col].notna()]
            top_ids = df.nlargest(top_n, fst_col)[id_col].astype(str).tolist()
            pair = ".".join(Path(f).stem.replace(".fst", "").split(".")[-2:])
            per_pair[pair] = top_ids
            pool_ids.update(top_ids)

        mapped = [snp_col_index[i] for i in pool_ids if i in snp_col_index]
        n_unmapped = len(pool_ids) - len(mapped)
        if n_unmapped:
            raise KeyError(
                f"{n_unmapped}/{len(pool_ids)} FST variant IDs not in numpy columns "
                f"— ID-namespace guard failed"
            )
        indices = np.array(sorted(mapped), dtype=np.int64)

        if return_details:
            return indices, {
                "pool_ids": pool_ids,
                "per_pair": per_pair,
                "n_pairs": len(fst_files),
            }
        return indices
    finally:
        if cleanup and own_workdir:
            shutil.rmtree(workdir, ignore_errors=True)
