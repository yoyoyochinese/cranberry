"""Phase 6 axis 2: configuration diagrams (signal flow + optical layout) of every scheme / fusion on F_liu20, through
the CLI `python -m cranberry diagram` (D84-3, D87-2), into docs/diagrams/phase6/axis2/.

  uv run python reviews/scripts/p6_axis2_diagrams.py
"""

from __future__ import annotations

import importlib.util
import subprocess
import sys
from pathlib import Path

_sp = importlib.util.spec_from_file_location("p6_axis2_lib", Path(__file__).with_name("p6_axis2_lib.py"))
lib = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(lib)

OUT = lib.ROOT / "docs" / "diagrams" / "phase6" / "axis2"
LABELS = ("base", "L-LOdeliv", "pair_dk16", "P1_calframe", "P2_borealis_si", "P3_wdm_chop6k", "P4_wdm_R0.3", "P5_dbs",
          "P6_steps", "P8_two_tone", "P8_intermittent", "F_P3+P5", "F_P2+P5", "F_P1slow+Liu", "F_P8int+P3")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    base = "F_liu20"
    cfg = lib.configs(base)
    for label in LABELS:
        ovs = cfg[label][0]
        cmd = [sys.executable, "-m", "cranberry", "diagram", str(lib.ROOT / lib.fs.nominal(base)), "--out", str(OUT),
               "--name", f"{base}_{label.replace('+', 'plus')}", "-f", "svg"]
        for o in ovs:
            cmd += ["--overlay", o]
        r = subprocess.run(cmd, capture_output=True, text=True, cwd=lib.ROOT)
        print(label, "ok" if r.returncode == 0 else f"ERROR {r.stderr[-600:]}", flush=True)
    # P7: the identification wiring
    cmd = [sys.executable, "-m", "cranberry", "diagram", str(lib.ROOT / lib.fs.nominal(base)), "--out", str(OUT),
           "--name", f"{base}_P7_ident", "-f", "svg", "--overlay", "ident/wire_F", "--overlay",
           "ident/in_lock_stepped"]
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=lib.ROOT)
    print("P7_ident", "ok" if r.returncode == 0 else f"ERROR {r.stderr[-600:]}", flush=True)


if __name__ == "__main__":
    main()
