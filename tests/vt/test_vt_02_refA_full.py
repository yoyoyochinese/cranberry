"""VT-02 (prompt §7-1, D48): the full `refA` simulation reproduces R′− = −10.15 ± 0.10 dB and
R′+ = +20.56 ± 0.10 dB with the lumped-loss reading (waveguide loss as a lumped element) and the
reference loss budget [Ref A §4.2], with the lock residual tuned to θ ≈ 9 mrad.

Sources: S §4 VT-02; D48; T04 §11; T12 V4; T01 §7.3. Pipeline: presets/reference/refA.yaml (phase-
detection-OPA topology L-A; how 9 mrad is reached is documented in its header).

Tolerance ±0.10 dB (D48), justified as in T04 §11: (1) statistical error of θ² from a finite record —
here the residual is dominated by a flat spectrum up to ~MHz, so a 0.2 s record at 20 MHz has
N_eff ~ 10⁵ and the seed spread of R′− is ~0.001 dB (checked in test_refA_convergence.py); (2) the loss
model (lumped vs distributed, T01: R′+ −0.27 dB for the distributed option, not used here; the 8 % is
reproduced to 1e−12 by the registry); (3) ESA averaging type (log vs power, ≤ 0.033 dB at 9 mrad,
T04 §5.2); (4) radiometer noise of the simulated trace mean (≈ 0.011 dB, T04 §5.3). The linear sum
≈ 0.11 dB is rounded to D48's 0.10 dB. The checks below are tighter where they test the
implementation rather than the measurement emulation: η = 0.92 to 1e−12, θ = 9.0 ± 0.2 mrad,
noise-level R′∓ within ±0.10 dB and the ESA zero-span means within ±0.10 dB.
"""

import math

import pytest

from cranberry import config
from cranberry.blocks.meas.esa import radiometer_std_db
from cranberry.runner.run import run

pytestmark = [pytest.mark.numba, pytest.mark.integration]
PRESET = "presets/reference/refA.yaml"


@pytest.fixture(scope="module")
def res():
    return run(config.load(PRESET))


def test_loss_budget_is_ref_a(res):
    reg = res.offline.loss_registry
    assert reg.eta("sqz", "bhd") == pytest.approx(0.92, rel=1e-12)  # L = 8 % [Ref A §4.1]
    ids = {e.id: e.eta for e in reg.elements_on("sqz", "bhd")}
    assert ids["sqz_opa.wg"] == pytest.approx(0.98)            # waveguide 2 % [Ref A §4.2]
    assert ids["optics.components"] == pytest.approx(0.99)     # mirrors, lenses 1 % [Ref A §4.2]
    assert ids["bhd_det.qe"] == pytest.approx(0.99)            # PD QE ≈ 99 % [Ref A §3.1]
    assert 1 - ids["bhd_det.electronic"] == pytest.approx(1.58e-3, rel=0.01)   # 28 dB clearance (T12 V8)
    assert 1 - ids["mm.overlap"] == pytest.approx(0.0407, abs=5e-4)  # mode mismatch ≈ 4 % by subtraction [Ref A §4.2]
    assert reg.eta("sqz", "gbs") == pytest.approx(0.98 * 0.99)  # BHD-only elements excluded (D46)


def test_r_prime_lumped_loss(res):
    nl = res.metrics.noise_levels["i.port0"]
    assert nl["sigma_total_rad"] == pytest.approx(9.0e-3, abs=0.2e-3)
    assert nl["r_minus_db"] == pytest.approx(-10.15, abs=0.10)
    assert nl["r_plus_db"] == pytest.approx(20.56, abs=0.10)
    # without jitter the source + registry give the §7-1 R∓ (T04 §9.1: −10.582 / +20.556 dB)
    assert nl["v_minus_db"] == pytest.approx(-10.582, abs=1e-3)
    assert nl["v_plus_db"] == pytest.approx(20.556, abs=1e-3)
    assert res.metrics.metric_ii_status == "not_applicable"


def test_esa_zero_span_traces(res):
    """Fig. 3 locked traces through measure.esa (zero span, 3 MHz, RBW 1 MHz, VBW 100 Hz)."""
    sq, anti = res.measurements["esa"]["squeezed"], res.measurements["esa"]["anti_squeezed"]
    assert sq["mean_db"] == pytest.approx(-10.15, abs=0.10)
    assert anti["mean_db"] == pytest.approx(20.56, abs=0.10)
    # trace scatter is radiometer noise (T04 §5.3; Ref A §3.2 reports 0.07 dB), not phase noise
    sd = radiometer_std_db(1e6, 100.0, 1.0, "gaussian", "power")
    assert sq["std_db"] == pytest.approx(sd, rel=0.25)
    assert sd == pytest.approx(0.07, abs=0.01)


def test_lock_points(res):
    """Static lock-point offset of Δ: the 21 mrad non-common probe phase (f_sh τ difference of pr_bhd and
    pr_det) reaches Δ scaled by A−/A+ = e^{−2r} (T14 §3.4): ε ≈ 0.021/123 = 0.17 mrad."""
    nl = res.metrics.noise_levels["i.port0"]
    u0 = 2 * math.pi * 1e6 * 1.0 * 1.0003 / 299792458.0  # SQ segment: pr_bhd − pr_det static f_off·τ
    eps = math.atan(math.exp(-2 * math.sqrt(9.06 * 0.64)) * math.tan(u0))
    assert nl["mean_loop_rad"] == pytest.approx(eps, rel=0.1)
