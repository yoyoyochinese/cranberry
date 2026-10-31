"""M5 review: the ``lo_temporal_shape`` block's fitted overlaps on the Borealis fixture (x_ss = 2.4,
δt = 1/(40κ_c)) vs D60-5 (Gaussian 0.955, T17 0.951) and T17 §6.3 (rect 0.864, one-sided exp 0.841–0.845).

Run: uv run python reviews/scripts/m5_lo_shape_block.py
"""
import numpy as np

from cranberry.blocks.los.lo_shape import make_lo_temporal_shape
from cranberry.blocks.source.pulsed_testing import BOREALIS

if __name__ == "__main__":
    kc = np.pi * 66e6
    src = {"block": "pulsed_opo_source", "params": {**BOREALIS, "x_ss": 2.4}}
    eta = {}
    for sh in ("gauss", "rect", "one_sided_exp"):
        p = {"source": src, "shape": sh, "fit": True, "observer": "bhd", "source_tag": "t", "toggles": {"shape_mismatch": True}}
        if sh == "one_sided_exp":
            p["decay_rate_rad_s"] = kc
        b = make_lo_temporal_shape(p)
        eta[sh] = b.eta_t
        print(f"{sh:14s} eta_t = {b.eta_t:.4f}  params·κ_c = {[round(x * kc, 3) for x in b.shape_params]}")
    assert abs(eta["gauss"] - 0.955) < 1.5e-3 and abs(eta["gauss"] - 0.951) > 3e-3  # D60-5 confirmed, T17 refuted
    assert abs(eta["rect"] - 0.866) < 2e-3 and 0.839 < eta["one_sided_exp"] < 0.846
