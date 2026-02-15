#!/usr/bin/env python3
"""Main driver for example workflows in Analytical-mantle-flow-modeling.

Modes:
- quick: validate precomputed outputs and repo structure only (default).
- full: execute pressure + age + dip-comparison scripts.
"""

from __future__ import annotations

import argparse
import importlib
import os
import shutil
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent
FLOW_DIR = ROOT / "flow_computations"


def _check_python_deps() -> tuple[bool, list[str]]:
    required = [
        "numpy",
        "scipy",
        "matplotlib",
        "mpl_toolkits.basemap",
        "shapely",
        "xarray",
        "netCDF4",
        "geographiclib",
    ]
    missing = []
    for mod in required:
        try:
            importlib.import_module(mod)
        except Exception:
            missing.append(mod)
    return len(missing) == 0, missing


def _check_external_tools() -> tuple[bool, list[str]]:
    required = ["convert"]
    missing = [tool for tool in required if shutil.which(tool) is None]
    return len(missing) == 0, missing


def quick_validate() -> int:
    checks = [
        ROOT / "dip_observations" / "dip_catalogues" / "Slab2_const-depth" / "AllDips.txt",
        FLOW_DIR / "inputs" / "Subbon_Slab2.0Final_NoJapTail_nnr_FS.inp",
        FLOW_DIR / "inputs" / "Subfil_Slab2.0Final_NoJapTail_nnr_FS.inp",
        FLOW_DIR / "inputs" / "Subgrd.inp",
        FLOW_DIR / "inputs" / "Subgrd_Fast.inp",
        FLOW_DIR / "text_files" / "Slab2.0Final_NoJapTail_nnr_FS.3e+20_VcSlabFlux_width500000.0_alpha0.0.NoTailFlux" / "DP.txt",
        FLOW_DIR / "text_files" / "Slab2.0Final_NoJapTail_nnr_FS.3e+20_VcSlabFlux_width500000.0_alpha0.0.NoTailFlux" / "Pcoeff.txt",
        FLOW_DIR / "plots" / "dip_comparisons" / "Slab2.0Final_NoJapTail_nnr_FS.3e+20_VcSlabFlux_width500000.0_alpha0.0.NoTailFlux.fact1.327.png",
    ]
    missing = [p for p in checks if not p.exists()]
    if missing:
        print("Quick validation failed. Missing files:")
        for path in missing:
            print("-", path)
        return 1

    print("Quick validation passed.")
    print("Key precomputed artifacts are present.")
    return 0


def full_run(python_cmd: str) -> int:
    py_ok, py_missing = _check_python_deps()
    tool_ok, tool_missing = _check_external_tools()

    if not py_ok or not tool_ok:
        print("Full workflow prerequisites are missing.")
        if py_missing:
            print("Missing Python modules:")
            for m in py_missing:
                print("-", m)
        if tool_missing:
            print("Missing external tools:")
            for t in tool_missing:
                print("-", t)
        print("Run with '--mode quick' to validate shipped outputs without rerunning models.")
        return 2

    env = os.environ.copy()
    env.setdefault(
        "DIPS_OBS_TXT",
        str(ROOT / "dip_observations" / "dip_catalogues" / "Slab2_const-depth" / "AllDips.txt"),
    )
    env.setdefault("MPLCONFIGDIR", "/tmp/mplcfg")
    env.setdefault("XDG_CACHE_HOME", "/tmp")

    commands = [
        [
            python_cmd,
            "global_pressure_withPressurePlot.py",
            "Slab2.0Final_NoJapTail_nnr_FS",
            "3.0e20",
            "2",
            "500000",
            "0",
            "1",
            "Subgrd_Fast.inp",
            "4.0e20",
        ],
        [python_cmd, "get_SPages.py", "Slab2.0Final_NoJapTail_nnr_FS"],
        [
            python_cmd,
            "plot_DipComparison_varyDPfactor.py",
            "Slab2.0Final_NoJapTail_nnr_FS",
            "3.0e20",
            "2",
            "500000",
            "0",
            "1",
            "Subgrd_Fast.inp",
            "2",
            "12",
            "4",
            "5",
        ],
    ]

    try:
        for cmd in commands:
            print()
            print("Running:", " ".join(cmd))
            subprocess.run(cmd, cwd=str(FLOW_DIR), check=True, env=env)
    except subprocess.CalledProcessError as exc:
        print("Workflow failed with exit code", exc.returncode)
        return exc.returncode or 1

    print()
    print("Full workflow completed.")
    print("Check outputs in:")
    print("-", FLOW_DIR / "text_files")
    print("-", FLOW_DIR / "plots" / "pressure_fields")
    print("-", FLOW_DIR / "plots" / "dip_comparisons")
    return 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run example workflows for the mantle-flow modeling repo.")
    parser.add_argument(
        "--mode",
        choices=["quick", "full"],
        default="quick",
        help="quick=validate shipped outputs; full=rerun pressure/age/dip workflow",
    )
    parser.add_argument(
        "--python",
        default=sys.executable,
        help="Python executable used for full mode (default: current interpreter)",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.mode == "quick":
        return quick_validate()
    return full_run(args.python)


if __name__ == "__main__":
    raise SystemExit(main())
