"""Phase 5 M1 review (finding M1-A-01): laser x filter-cavity out-of-band term double counted?

Claim checked: the in-loop white-FM random walk is an *exact sample* of the continuous Brownian phase
(laser.py: Var[phi(t) - phi(t - k h)] = 2 pi^2 h0 k h), and filter_cavity's discrete low pass acting on it
gives Var(H_FC phi - phi) = 2 pi^2 h0 h (1 + a) / (4 (1 - a)) = pi h0 / dnu_FC + O(h^2) -- the FULL
continuous variance (aliasing included). The ledger (oob_pairs {dnu_fc_hz}) then adds
int_{f_N}^{f_max} |1 - H_FC|^2 h0 / f^2 df ~= h0 / f_N on top: a double count.

Part 1: closed form of the discrete filter (no simulator).
Part 2 (--sim): the refA preset with only the white FM on (mech off, sensor noise off, quantisation off),
        f_sim = 10 / 20 / 40 MHz: loop variance vs continuous truth incl. the LO-squeezed delay, and the
        ledger's sigma^2_HF.
Run: uv run python reviews/scripts/m1_fc_ledger_double_count.py [--sim]
"""

import copy
import math
import sys

import numpy as np

H0 = 200.0 / math.pi
DNU = 3.3e6
C = 299792458.0


def discrete_var(f_sim: float) -> float:
    h = 1.0 / f_sim
    a = math.exp(-math.pi * DNU * h)
    return 2 * math.pi**2 * H0 * h * (1 + a) / (4 * (1 - a))


def continuous_fc_delay(tau: float, f_hi: float = 1e13) -> float:
    """int_0^f_hi |H e^{-i w tau} - 1|^2 h0 / f^2 df (LO through FC, delayed by tau vs the squeezed light)."""
    f = np.geomspace(1e-2, f_hi, 2_000_001)
    hfc = 1.0 / (1.0 + 2j * f / DNU)
    s = np.abs(hfc * np.exp(-2j * np.pi * f * tau) - 1.0) ** 2 * H0 / f**2
    return float(np.trapezoid(s * f, np.log(f)))


print("Part 1: discrete FC filter on the exactly sampled random walk")
print(f"  continuous pi h0/dnu = {math.pi * H0 / DNU:.4e} rad^2")
for fs in (1e7, 2e7, 4e7):
    print(f"  f_sim={fs:.0e}: discrete Var = {discrete_var(fs):.4e} rad^2; ledger adds ~h0/f_N = {H0 / (fs / 2):.3e}"
          f" (+{100 * H0 / (fs / 2) / (math.pi * H0 / DNU):.1f} % of the FC term)")
assert abs(discrete_var(2e7) / (math.pi * H0 / DNU) - 1) < 3e-2  # O(h^2): eps^2/12, eps = pi dnu h



def discrete_var_with_delays(f_sim: float, n: int = 20000) -> float:
    """Expected in-loop Var(Delta) of refA's discrete model (no lock): LO = FC low pass of phi(t - tau_lo)
    (linear interpolation, as filter_cavity / beam_assembler), squeezed = phi(t - tau_sq) (linear interp.)."""
    h = 1.0 / f_sim
    a = math.exp(-math.pi * DNU * h)

    def interp(tau):
        k = int(tau / h)
        mu = tau / h - k
        g = np.zeros(n)
        g[k] += 1 - mu
        g[k + 1] += mu
        return g

    xb = interp(TAU_LO)
    y = np.zeros(n)
    py = px = 0.0
    for i in range(n):
        y[i] = a * py + (1 - a) * 0.5 * (xb[i] + px)
        py, px = y[i], xb[i]
    gw = np.cumsum(y - interp(TAU_SQ))  # response to one random-walk increment
    return 2 * math.pi**2 * H0 * h * float(np.sum(gw**2))


TAU_LO = (5.0 * 1.468 + 0.3 * 1.0003) / C
TAU_SQ = (1.0 * 1.0003 + 0.3 * 1.0003 + 0.2 * 1.0003 + 0.045 * 2.2 + 1.0 * 1.0003) / C
TRUTH = continuous_fc_delay(TAU_LO - TAU_SQ)
print(f"\nPart 1b: expected in-loop variance with the LO-sqz delay (linear interpolation), no lock; truth {TRUTH:.4e}")
for fs in (1e7, 2e7, 4e7, 8e7, 1.6e8):
    v = discrete_var_with_delays(fs)
    print(f"  f_sim={fs:.1e}: in-loop {v:.4e}  + ledger h0/f_N = {v + H0 / (fs / 2):.4e}  (/truth {(v + H0 / (fs / 2)) / TRUTH:.4f})")

if "--sim" not in sys.argv:
    sys.exit(0)

from cranberry import config  # noqa: E402
from cranberry.config.yaml_io import load_yaml  # noqa: E402
from cranberry.runner.run import run  # noqa: E402

PRESET = "presets/reference/refA.yaml"
base = load_yaml(PRESET)

# path delays of the bookkeeping beams (preset segment table)
tau_lo = (5.0 * 1.468 + 0.3 * 1.0003) / C
tau_sq = (1.0 * 1.0003 + 0.3 * 1.0003 + 0.2 * 1.0003 + 0.045 * 2.2 + 1.0 * 1.0003) / C
truth = continuous_fc_delay(tau_lo - tau_sq)
print(f"\nPart 2: refA, white FM only. delay LO - sqz = {1e9 * (tau_lo - tau_sq):.2f} ns; "
      f"continuous truth (FC + delay, all f) = {truth:.4e} rad^2 -> {1e3 * math.sqrt(truth):.3f} mrad")


def blk(tree, name):
    for b in tree["pipeline"]:
        if b.get("name") == name:
            return b
    raise KeyError(name)


for fs in (1e7, 2e7, 4e7):
    t = copy.deepcopy(base)
    t["run"]["f_sim_hz"] = fs
    t["run"]["duration_s"] = 0.05
    t["run"]["chunk_samples"] = 262144
    blk(t, "mech")["params"]["toggles"]["mech"] = False
    for n in ("pd1", "bhdbeat"):
        blk(t, n)["params"]["toggles"].update(shot_noise=False, electronic_noise=False)
    for n in ("adc1", "adc2", "dac1", "dac2"):
        blk(t, n)["params"]["toggles"]["quantization"] = False
    t["measure"].pop("esa", None)
    res = run(config.load(t))
    nl = res.metrics.noise_levels["i.port0"]
    loop2 = nl["sigma_loop_rad"] ** 2
    hf2 = nl["sigma_hf_rad"] ** 2
    print(f"  f_sim={fs:.0e}: loop {loop2:.4e}  ledger HF {hf2:.4e}  total {loop2 + hf2:.4e}"
          f"  (total/truth = {(loop2 + hf2) / truth:.4f}, loop/truth = {loop2 / truth:.4f});"
          f"  sigma_total = {1e3 * nl['sigma_total_rad']:.3f} mrad, R'- = {nl['r_minus_db']:.4f} dB")
