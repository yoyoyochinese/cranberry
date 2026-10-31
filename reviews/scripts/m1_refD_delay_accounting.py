"""Phase 5 M1 review (Ref D): where the 1.5 µs of D61-2 actually sit, and the −180° frequency.

1. Phase of the ZOH-discretised 914 (refD) minus the continuous two-pole phase: the extra lag is the
   half-sample ZOH interface delay h/2 = 0.5 µs at 1 MHz (it is part of the 1 µs "500 kHz ZOH").
2. The same with the 914 ``poles`` toggle off (or ``discretisation: tustin``): the half sample is gone,
   so the loop is 0.5 µs short of τ_L.
3. Continuous check of T15 §2.7(b): f_180 of 120 kHz + 100 kHz poles + τ_L = 5.67 / 5.7 µs, and of the
   8.3 µs pure delay.

Run: uv run python reviews/scripts/m1_refD_delay_accounting.py
"""

import importlib.util
import math
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[2]
_sp = importlib.util.spec_from_file_location("refD_support", ROOT / "tests" / "integration" / "refD_support.py")
rs = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(rs)


def main():
    f = np.array([2.0e4, 4.0e4, 6.0e4])
    s = 2j * np.pi * f
    cont = 1 / ((1 + s / (2 * np.pi * 2.95e4)) * (1 + s / (2 * np.pi * 1.12e5)))
    for label, tog, blocks in (("zoh (refD)", {}, None),
                               ("tustin", {}, [{"name": "amp914", "params": {"discretisation": "tustin"}}]),
                               ("poles off", {"amp914.poles": False}, None)):
        spec = rs.load(toggles=tog or None, blocks=blocks)
        amp = rs.built_leaves(spec)[0]["amp914"].freq_response(f)
        ref = cont if label != "poles off" else np.ones_like(cont)
        extra = -np.angle(amp / ref) / (2 * np.pi * f)
        print(f"914 {label:10s}: extra delay vs continuous at 20/40/60 kHz = {np.round(extra * 1e6, 3)} µs")

    def f180(tau, poles=True):
        def ph(x):
            p = math.atan(x / 1.2e5) + math.atan(x / 1.0e5) if poles else 0.0
            return p + 2 * math.pi * x * tau - math.pi
        return brentq(ph, 1e4, 2e5)

    print(f"T15 fit: τ = 5.67 µs → f_180 = {f180(5.67e-6) / 1e3:.2f} kHz; 5.7 µs → {f180(5.7e-6) / 1e3:.2f} kHz; "
          f"8.3 µs pure → {f180(8.3e-6, False) / 1e3:.2f} kHz; 5.7 + 1.5 µs (double count) → {f180(7.2e-6) / 1e3:.2f} kHz")


if __name__ == "__main__":
    main()
