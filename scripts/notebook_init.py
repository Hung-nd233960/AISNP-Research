"""
One-call setup for Jupyter notebooks.

Usage (first code cell of every notebook):
    import sys, os
    sys.path.insert(0, str(__import__('pathlib').Path.cwd().parent.parent / 'scripts')
                       if (__import__('pathlib').Path.cwd() / 'scripts').exists() is False
                       else str(__import__('pathlib').Path.cwd() / 'scripts'))

Simpler — just call setup():
    from notebook_init import setup
    cfg, PATHS, POPULATIONS, HARD_FILTERS, SITUATIONAL_FILTERS, ML = setup()
"""

import os
import sys
import warnings
from pathlib import Path


def setup(suppress_warnings: bool = True):
    """
    Locate project root, fix working directory, add scripts/ to sys.path,
    and return the config module plus its main constants.

    Returns:
        cfg, PATHS, POPULATIONS, HARD_FILTERS, SITUATIONAL_FILTERS, ML
    """
    # Walk up from CWD until we find the project root (contains scripts/config.py)
    root = Path.cwd()
    while not (root / "scripts" / "config.py").exists():
        if root == root.parent:
            raise RuntimeError(
                "Cannot find project root. Launch Jupyter from the project directory."
            )
        root = root.parent

    os.chdir(root)
    scripts_dir = str(root / "scripts")
    if scripts_dir not in sys.path:
        sys.path.insert(0, scripts_dir)

    if suppress_warnings:
        warnings.filterwarnings("ignore")

    import importlib
    import config as cfg
    importlib.reload(cfg)  # re-read paths.local.yaml and config.py on every setup() call

    print(f"Project root : {root}")
    print(f"genomes_data : {cfg.PATHS._genomes_data}")
    print(f"output root  : {cfg.PATHS._output_root}")

    return cfg, cfg.PATHS, cfg.POPULATIONS, cfg.HARD_FILTERS, cfg.SITUATIONAL_FILTERS, cfg.ML
