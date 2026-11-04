"""Build check of every axis-2 configuration on every base (short run): prints open axes / errors / CHK codes.

  uv run python reviews/scripts/p6_axis2_probe.py [label ...]
"""

from __future__ import annotations

import importlib.util
import sys
import warnings
from pathlib import Path

_sp = importlib.util.spec_from_file_location("p6_axis2_lib", Path(__file__).with_name("p6_axis2_lib.py"))
lib = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(lib)

from cranberry.runner.run import run  # noqa: E402

warnings.simplefilter("ignore")
labels = sys.argv[1:]
for base in lib.BASES:
    for label in labels or lib.configs(base):
        try:
            s = lib.build(base, label, duration_s=0.02, seed=1)
            r = run(s)
            q = r.summary()
            print(f"OK  {base:8s} {label:22s} valid={r.valid} ii={1e3 * q.get('ii.port1.sigma_total_rad', float('nan')):7.1f} mrad "
                  f"L_k={s.offline['params']['L_k_m']:.0f} chk={sorted({e['code'] for e in r.events})}", flush=True)
        except Exception as e:  # noqa: BLE001
            print(f"ERR {base:8s} {label:22s} {type(e).__name__}: {str(e)[:400]}", flush=True)
