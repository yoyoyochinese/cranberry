"""Phase 6 axis 2, P7: in-lock identification (ident/wire_F + ident/in_lock_stepped, port 0 full stack) of the future
nominal L-Liu loop on F_liu20 / F_wg6, with the fpga_generic, Moku:Delta and DSC1 controllers and with the 915B
multimode lab overlay. The identified gain margin sets the upper BoTorch bound of the loop gain:
f_u,max = f_u · 10^((GM − 6 dB)/20) (6 dB kept as the robustness margin, the same |S| ≤ 2 rule as the servo-bump limit).

  uv run python reviews/scripts/p6_axis2_ident.py      → runs/p6_axis2/ident.json
"""

from __future__ import annotations

import importlib.util
import json
import math
import warnings
from pathlib import Path

import numpy as np

_sp = importlib.util.spec_from_file_location("p6_axis2_lib", Path(__file__).with_name("p6_axis2_lib.py"))
lib = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(lib)

from cranberry.runner.run import run  # noqa: E402

MARG = ("f_u_hz", "phase_margin_deg", "f_180_hz", "gain_margin_db")
CASES = {"fpga": ((), {}), "moku": (("controller/moku_delta",), {}), "dsc1": (("controller/dsc1",), lib.DSC1),
         "fpga_mm915b": (("lab/multimode915b_future",), lib.MM915B),
         # the 25-point 200 Hz – 50 kHz grid steps over the Q = 25 main mode: refined sweep 15–35 kHz (review I-3)
         "fpga_mm915b_refined_15-35kHz": (("lab/multimode915b_future",),
                                          {**lib.MM915B, "ident_f_lo_hz": 1.5e4, "ident_f_hi_hz": 3.5e4})}


def main():
    warnings.simplefilter("ignore")
    out = {}
    for base in ("F_liu20", "F_wg6"):
        for name, (ovs, pins) in CASES.items():
            spec = lib.fs.load(lib.fs.nominal(base), overlays=(*ovs, "ident/wire_F", "ident/in_lock_stepped"),
                               run={"seed": 1, "seed_mode": "common"})
            from cranberry.runner.sweep import expand
            spec = expand(spec, {k: [v] for k, v in lib._p({**pins, "port_k": 31}).items()}, limit=2)[0].spec
            r = run(spec)
            m = r.measurements["loop_identification"]["lock0"]
            row = {k: float(m[k]) for k in MARG} | {f"model_{k}": float(m.get(f"model_{k}", math.nan)) for k in MARG}
            row["coh_min"] = float(np.nanmin(m["coherence"]))
            f_u_nom = float(spec.offline["params"]["f_u_hz"])
            gm = row["gain_margin_db"]
            row["f_u_nominal_hz"] = f_u_nom
            row["f_u_max_hz"] = f_u_nom * 10 ** ((gm - 6.0) / 20.0) if math.isfinite(gm) else math.nan
            row["f_u_stab_hz"] = f_u_nom * 10 ** (gm / 20.0) if math.isfinite(gm) else math.nan   # BoTorch bound
            row["chk"] = sorted({e["code"] for e in r.events})
            out[f"{base}/{name}"] = row
            print(base, name, {k: round(v, 2) for k, v in row.items()}, flush=True)
    lib.OUT.mkdir(parents=True, exist_ok=True)
    (lib.OUT / "ident.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
