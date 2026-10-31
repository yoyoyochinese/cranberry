"""VT-16 (prompt §7-8, D51-4): the conventional-tap trade-off, both definitions, and that the
phase-detection-OPA topology removes the tap loss.

Sources: S §4 VT-16; T12 §6, V10; T14 §3.3 (table), §3.4. Presets: presets/reference/refA.yaml (L-A,
no tap) and presets/reference/refA_tap.yaml (L-tap, tap ratio r_t on SQ).

(a) T12 §6.1 model (V10): L_tot(r_t) from the runner's LossRegistry of `refA_tap` at each r_t, and
    θ² = θ0² + (π/2) f_BW S_φ,sens with S_φ,sens = κ_sn hν/(η_d r_t P_pr) + NEP²/(m r_t P_pr)²
    (θ0 = 9 mrad, f_BW = 10 kHz, η_d = 0.9, κ_sn = 1/2, m = 1, P_pr = 10 µW, NEP = 1 pW/√Hz):
    R′−(r_t) has an interior minimum at r_t ≈ 3×10⁻³ (within one grid point) of −9.94 ± 0.02 dB.
(a′) The same r_t sweep through the closed-loop simulation of `refA_tap`: the tap-PD sensor noise
    sets lock 1's residual (∝ 1/r_t) but reaches Δ only through lock 2, scaled by A−/A+ = 1/123
    (T14 §3.4, VT-15), so Δ's residual does not change with r_t (±2 %) and R′−(r_t) has no interior
    minimum — the (a) trade-off assumes that the tap SNR enters θ 1:1. This is the C23 question
    (Ref A credits the lock-1 SNR with 14 → 9 mrad); recorded for the Phase 5 review, not resolved here.
(b) θ fixed at the simulated residual (9.0 ± 0.2 mrad, VT-02; 9.11 mrad after the Phase 5 M1-A-01
    ledger correction, Δν_FC not re-tuned), T_tap = 0.1 / 1 / 3 %: R′− equals the single-mode reference
    at the simulated θ and L = 1 − 0.92(1 − r_t) within ±0.01 dB, and the T14 §3.3 table (θ = 9 mrad)
    −10.11 / −9.76 / −9.07 within ±0.02 dB (θ 9.0 → 9.11 mrad moves R′− by ≈ 0.012 dB). The PD-OPA
    topology (`refA`) has no tap element and R′− = the reference at L = 8 % (±0.01 dB; −10.15 ± 0.02 dB)
    whatever the probe power at its PD (the tap loss is gone).
"""

import math

import numpy as np
import pytest

from cranberry import config
from cranberry.config.yaml_io import Delete, Replace, load_yaml
from cranberry.core.units import V_VAC, lin_to_db
from cranberry.runner.run import run

pytestmark = [pytest.mark.numba, pytest.mark.integration]
H_NU = 6.62607015e-34 * 299792458.0 / 1545.32e-9


def tap_spec(r_t, duration_s, preset="refA_tap", toggles=None, params=None):
    t = load_yaml(f"presets/reference/{preset}.yaml")
    t.setdefault("run", {})["duration_s"] = duration_s  # refA_tap merges it onto refA's run (extends)
    if "extends" in t:
        t.setdefault("measure", {})["esa"] = Delete()   # no zero-span traces on short records
    else:
        t["measure"].pop("esa")
    if toggles:
        t["toggles"] = dict(toggles)
    for b in t.get("pipeline", []):
        b = b.value if isinstance(b, Replace) else b  # refA_tap's `!replace` element (pd1)
        if not isinstance(b, dict):
            continue
        if b.get("name") == "tap" and r_t is not None:
            b["params"]["tap_ratio"] = r_t
        for (name, key), v in (params or {}).items():
            if b.get("name") == name:
                b["params"][key] = v
    return config.load(t)


def r_minus_db(eta, theta):
    """§7-1 exact Gaussian average (T12 §1.5) for the refA source with total transmission η."""
    r = math.sqrt(9.06 * 0.64)
    rm, rp = (1 - eta) + eta * math.exp(-2 * r), (1 - eta) + eta * math.exp(2 * r)
    return float(lin_to_db((rm + rp) / 2 + (rm - rp) / 2 * math.exp(-2 * theta**2)))


def test_vt16a_t12_tradeoff_with_runner_losses():
    grid = [1e-4, 3e-4, 1e-3, 2e-3, 3e-3, 5e-3, 1e-2, 3e-2, 1e-1]
    out = []
    for r_t in grid:
        res = run(tap_spec(r_t, 2e-7), mode="reference", keep_state=False)  # offline stage + 4 samples
        eta = res.offline.loss_registry.eta("sqz", "bhd")
        assert eta == pytest.approx(0.92 * (1 - r_t), rel=1e-12)  # L → 1 − (1 − L)(1 − r_t) (T12 §6.1)
        p_tap = r_t * 10e-6
        s = 0.5 * H_NU / (0.9 * p_tap) + (1e-12 / (1.0 * p_tap)) ** 2
        theta = math.sqrt(0.009**2 + math.pi / 2 * 1e4 * s)
        out.append(r_minus_db(eta, theta))
    k = int(np.argmin(out))
    assert abs(k - grid.index(3e-3)) <= 1 and 0 < k < len(grid) - 1  # interior minimum near 3e-3
    assert out[grid.index(3e-3)] == pytest.approx(-9.94, abs=0.02)


def test_vt16a_closed_loop_has_no_snr_tradeoff():
    rows = []
    for r_t in (1e-4, 1e-3, 1e-2, 3e-2):
        # tap PD NEP 10 pW/√Hz (T12 §6.2 row; ℜ = 1 A/W): sensor-noise-limited lock 1 at small r_t
        res = run(tap_spec(r_t, 0.03, params={("pd1", "i_n_a_rthz"): 1e-11}), keep_state=False)
        nl = res.metrics.noise_levels["i.port0"]
        rows.append((r_t, res.series["lock.lock1.err"].std() / 2, nl["sigma_total_rad"], nl["r_minus_db"]))
    u = [r[1] for r in rows]
    th = [r[2] for r in rows]
    rm = [r[3] for r in rows]
    assert u[0] > 5 * u[-1]                               # lock 1 gets much noisier at small r_t
    np.testing.assert_allclose(th, th[-1], rtol=0.02)     # Δ does not see it (× A−/A+)
    assert all(a < b for a, b in zip(rm, rm[1:]))         # R′− monotonic: smallest tap is best


def _reference_r_minus(theta, loss):
    from cranberry.blocks.source.single_mode import make_single_mode_reference

    b = make_single_mode_reference(dict(loss=loss, alpha_per_w=9.06, pump_w=0.64, theta_rms_rad=theta,
                                        phi_opa0_rad=0.0, average="small_angle"))
    return b.r_prime_db()[0]


def test_vt16b_fixed_theta_and_pd_opa_topology():
    for r_t, target in ((1e-3, -10.11), (1e-2, -9.76), (3e-2, -9.07)):
        res = run(tap_spec(r_t, 0.05), keep_state=False)
        nl = res.metrics.noise_levels["i.port0"]
        assert nl["sigma_total_rad"] == pytest.approx(0.009, abs=2e-4)
        assert nl["eta"] == pytest.approx(0.92 * (1 - r_t), rel=1e-12)
        ref = _reference_r_minus(nl["sigma_total_rad"], 1 - 0.92 * (1 - r_t))
        assert nl["r_minus_db"] == pytest.approx(ref, abs=0.01)
        assert nl["r_minus_db"] == pytest.approx(target, abs=0.02)
    # phase-detection-OPA topology: no tap, the PD-OPA probe power does not touch the squeezed path
    for p_in in (1e-4, 1e-3, 1e-2):
        res = run(tap_spec(None, 0.05, preset="refA", params={("pd1", "p_in_w"): p_in}), keep_state=False)
        reg = res.offline.loss_registry
        assert not [e for e in reg.elements_on("sqz", "bhd") if e.id.startswith("tap")]
        nl = res.metrics.noise_levels["i.port0"]
        assert nl["eta"] == pytest.approx(0.92, rel=1e-12)
        assert nl["r_minus_db"] == pytest.approx(_reference_r_minus(nl["sigma_total_rad"], 0.08), abs=0.01)
        assert nl["r_minus_db"] == pytest.approx(-10.15, abs=0.02)
    assert V_VAC == 0.5
