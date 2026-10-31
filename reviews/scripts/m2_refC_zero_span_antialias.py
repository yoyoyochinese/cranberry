"""Phase 5 M2 Ref C review: the zero-span hook builds V′(t) from the RECORDED δ series (runner._zero_span),
so an anti-alias filter on the record (D68-5 `antialias: boxcar | fir`, which the RecordAliasWarning asks
for) averages the fast phase jitter out of δ BEFORE cos 2δ is taken. In the time-series dither overlay the
1 MHz dither (m/√2 = 0.0555 rad) lives only in δ(t), so the zero-span squeezed level then loses the dither.
Compares measure.esa.sqz_30m for antialias none / boxcar with the full-rate metric moments (v_prime of the
swept spectrum bin at 30 MHz uses the full-rate ⟨cos 2δ⟩ and is the reference).

Run: uv run python reviews/scripts/m2_refC_zero_span_antialias.py   (≈ 2 × 0.625 s at 4 MHz, RF-rate QNL)
"""
import importlib.util
import math
from pathlib import Path

import numpy as np

from cranberry.runner.run import run

ROOT = Path(__file__).resolve().parents[2]
_sp = importlib.util.spec_from_file_location("refC_support", ROOT / "tests" / "integration" / "refC_support.py")
rs = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(rs)

for aa in ("none", "boxcar"):
    spec = rs.load(params={"x": rs.X_T04}, toggles=rs.NOISE_OFF, overlays=(rs.TIMESERIES,), esa=True,
                   measure={"record": {"slots": ["lock.lockC.err"], "decim": 25, "antialias": aa}})
    res = run(spec)
    zs = res.measurements["esa"]["sqz_30m"]
    nl = res.metrics.noise_levels["i.port0"]

    sp = res.measurements["esa"]["spectrum_sqz"]
    i30 = int(np.argmin(np.abs(np.asarray(sp["f_hz"]) - 3e7)))   # esa_spec bins: 10 MHz … 1.3 GHz, 130 points
    print(f"antialias {aa:<7s}: zero span 30 MHz V′ = {10*math.log10(zs['r_lin']):+.3f} dB; swept-spectrum bin "
          f"(full-rate moments) = {10*math.log10(np.atleast_1d(sp['r_lin'])[i30]):+.3f} dB; "
          f"σ_total(30 MHz) {rs.at(nl, 'sigma_total_rad', 3e7)*1e3:.1f} mrad")
