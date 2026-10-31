"""VT-41 (loss budgets; 00_summary §4 VT-41 ← T12; PLAN §6.2). Ref C part (WP-I3); the Ref B budget (WP-I4)
and the 32-port / GBS-vs-BHD budgets (WP-I5) are added by their owners.

Ref C [Ref C §3, §4; T02 §6.2, V-C1]: η_tot = η_esc · (LossRegistry η at the BHD) with
* before the fibre 15 % "applied to both stations", incl. table-top propagation, mode matching and the BHD
  QE 97 % (escape efficiency 0.99 [Tohermes 2024 via T02] inside it): η_pre = 0.85;
* fibre in/out coupling + 1 km at 0.2 dB/km: 8 % "(0.4 dB)" [Ref C §4] vs the fitted total 21 % [Ref C §3,
  Fig. 5]. The preset keeps both (choice axis eta_fiber_link): 0.79/0.85 (total = the fit, 21 %) and 0.92
  (total 21.8 %); they differ by 0.8 percentage points.
Values are dark-noise corrected [Ref C Fig. 4 caption]: the BHD's electronic element is registered but
disabled (subtract_dark), and the 0.83 GHz interference dip does not enter η (it is a gain dip, T13 §1.4).
The LO-path actuators (EOM, stretcher, mirror) register no loss on the squeezed light (T16 §3.2). The η of
metric (i) equals the registry's (no second bookkeeping).
"""

import importlib.util
from pathlib import Path

import pytest

from cranberry.runner.run import run

_sp = importlib.util.spec_from_file_location("refC_support", Path(__file__).parents[1] / "integration" / "refC_support.py")
rs = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(rs)

pytestmark = [pytest.mark.numba, pytest.mark.integration]


def offline(params=None):
    res = run(rs.load(params=params, run={"duration_s": 0.025}))
    return res, res.offline.loss_registry, res.offline.source_outputs["opo"]


def test_vt41_refC_budget():
    res, reg, src = offline()
    p = rs.params_of(res.spec)
    rho = p["escape_efficiency"]
    els = {e.id: e for e in reg.elements_on("sqz", "bhd", include_disabled=True)}
    assert set(els) == {"optics.stations", "optics.fiber_link", "bhd_det.qe", "bhd_det.electronic", "bhd_det.imbalance"}
    assert not els["bhd_det.electronic"].enabled and not els["bhd_det.imbalance"].enabled
    assert els["bhd_det.qe"].eta == pytest.approx(0.97, rel=1e-12)
    eta = reg.eta("sqz", "bhd")
    assert rho * eta == pytest.approx(0.79, rel=1e-6)                                    # fitted 21 %
    pre = rho * eta / els["optics.fiber_link"].eta
    assert pre == pytest.approx(0.85, rel=1e-12)                                         # 15 % before the fibre
    assert res.metrics.noise_levels["i.port0"]["eta"] == pytest.approx(eta, rel=1e-12)
    assert {e.id for e in reg.elements_on("lo", "bhd")} == {"bhd_det.qe"}                # only the detector; no LO-path element
    # the source carries the escape efficiency (ρ inside cov), not the registry
    assert src.basis == "cw_sideband"


def test_vt41_refC_loss_analysis_variant():
    res, reg, _ = offline({"eta_fiber_link": 0.92})
    rho = rs.params_of(res.spec)["escape_efficiency"]
    total = rho * reg.eta("sqz", "bhd")
    assert total == pytest.approx(0.85 * 0.92, rel=1e-12)                               # 21.8 %
    assert total == pytest.approx(0.79, abs=0.01)                                        # consistent with the fit (0.8 points)
