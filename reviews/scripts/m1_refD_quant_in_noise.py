"""Phase 5 M1 review (Ref D): does the DAC limit cycle survive the nominal noise? (D61-3 adjudication)

refD_nominal with every noise source on (laser white FM, mechanical, thermal drift, RIN, PDB, ADC noise),
DAC quantisation on vs off, for ΔL = 9 m (laser-dominated) and ΔL = 0 (laser term vanishes), rigid and
nominal-resonance 915B. Same seed, so the difference is the quantiser. Reports the true-phase residual
rms (last 30 ms) and the excess variance attributable to the DAC.

Run: uv run python reviews/scripts/m1_refD_quant_in_noise.py
"""

import importlib.util
import math
from pathlib import Path

import numpy as np

from cranberry.runner.run import run

ROOT = Path(__file__).resolve().parents[2]
_sp = importlib.util.spec_from_file_location("refD_support", ROOT / "tests" / "integration" / "refD_support.py")
rs = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(rs)

REC = {"record": {"slots": ["lock.lockD.err"], "decim": 1}}


def residual(dl, resonance, dac_q, h_ref=None):
    params = {"dL_m": dl}
    if h_ref is not None:
        params["h_ref_rad2_hz"] = h_ref
    spec = rs.load(params=params, run={"duration_s": 0.08}, measure=REC,
                   toggles={"dsc1.dac.quantization": dac_q, "stretcher915b.resonance": resonance})
    e = run(spec).series["lock.lockD.err"][-30000:]
    return float(np.std(e))


def main():
    lsb = 2 * math.pi * 16e-6 / (1.5488e-6 * 10) * 8 / 4096
    print(f"1 LSB = {lsb * 1e3:.2f} mrad")
    for dl in (9.0, 0.0):
        for res in (False, True):
            on, off = residual(dl, res, True), residual(dl, res, False)
            ex = math.sqrt(max(on**2 - off**2, 0.0))
            print(f"dL = {dl:4.1f} m, {'20 kHz Q10' if res else 'rigid    '}: rms DAC-q on {on * 1e3:7.3f} mrad, "
                  f"off {off * 1e3:7.3f} mrad, excess {ex * 1e3:6.3f} mrad ({ex / lsb:.2f} LSB)")
    # quiet environment (lowest h_ref of the sweep), ΔL = 0: closest to the noise-free limit cycle
    for res in (False, True):
        on, off = residual(0.0, res, True, h_ref=1e-3), residual(0.0, res, False, h_ref=1e-3)
        ex = math.sqrt(max(on**2 - off**2, 0.0))
        print(f"dL = 0, h_ref = 1e-3, {'20 kHz Q10' if res else 'rigid    '}: on {on * 1e3:.3f}, off {off * 1e3:.3f}, "
              f"excess {ex * 1e3:.3f} mrad ({ex / lsb:.2f} LSB)")


if __name__ == "__main__":
    main()
