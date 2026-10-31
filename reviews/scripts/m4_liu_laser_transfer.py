"""M4 review (D62-5): in-band laser-phase transfer of the Liu-type lock, port k vs port 0.

Independent derivation from the conventions' path model (conv §7, §9; T10 §2.3):
  squeezed (GBS, t):  φ_L(t − T_a − τ_k) + a_k(t − τ_k + τ_xa)
  probe (detector):   φ_L(t − T_pd − τ_k − T_px) + a_k(t − T_px − τ_xa)      (T_b = T_pd + T_px)
  reference:          φ_L(t − T_ref)
  ideal lock: probe − reference = const  →  a_k(s) = φ_L(s + T_px + τ_xa − T_ref) − φ_L(s + τ_xa − T_pd − τ_k)
  ⇒ θ_k/φ_L = e^{−iω(T_a+τ_k)} − e^{−iω(T_b − T_px + 2τ_k − 2τ_xa)} + e^{−iω(T_ref − T_px + τ_k − 2τ_xa)}
With τ_xa = T_px = 0 this is B8's H (D62-5). Parts:
  A  the general formula vs a sample-level replay through the project's F6 path (beam_assembler,
     emission times), with an extraction segment T_px and the actuator upstream or downstream;
  B  T10 eq. (13) vs B8's H; which delays null H; one pair only for port-independent delays;
  C  finite loop gain: residual = H_ideal + (P_k/(1+L_k) − P_0/(1+L_0))·(A_s/A_p) — compare magnitudes;
  D  linewidth requirement numbers (white FM).

Run: uv run python reviews/scripts/m4_liu_laser_transfer.py
"""
import math

import numpy as np

import cranberry.path.testing as pt
from cranberry.path.topology import C_M_PER_S

N_G = 1.468
FS = 1e9


def seg_len(samples):
    return samples * C_M_PER_S / (N_G * FS)


def dtft(h, f):
    n = np.arange(len(h))
    return np.exp(-2j * np.pi * np.outer(np.atleast_1d(f), n) / FS) @ h


def e(f, d):
    return np.exp(-2j * np.pi * f * d / FS)


def H_general(f, Ta, Tpd, Tpx, Tref, t0, tk, xa0, xak):
    """θ_k − θ_0 per unit φ_L (ideal locks), delays in samples; xa_j = τ_xa of port j."""
    def th(tj, xa):
        return e(f, Ta + tj) - e(f, Tpd + 2 * tj - 2 * xa) + e(f, Tref - Tpx + tj - 2 * xa)
    return th(tk, xak) - th(t0, xa0)


def H_b8(f, Ta, Tb, Tref, t0, tk):
    return (e(f, tk) - e(f, t0)) * (e(f, Ta) + e(f, Tref) - e(f, Tb) * (e(f, tk) + e(f, t0)))


def topo(Ta, Tpd, Tpx, Tref, t0, tk, actuator):
    segs = [pt.seg_spec("up", seg_len(Ta)), pt.seg_spec("pd", seg_len(Tpd)), pt.seg_spec("px", seg_len(Tpx)),
            pt.seg_spec("rf", seg_len(Tref)), pt.seg_spec("c0", seg_len(t0)), pt.seg_spec("ck", seg_len(tk)),
            pt.seg_spec("a0", 0.0, drivers={"stretcher": "uniform"}),
            pt.seg_spec("ak", 0.0, drivers={"stretcher": "uniform"})]
    s = pt.step_spec

    def sq(j):
        return [s("up"), s(f"a{j}"), s(f"c{j}")] if actuator == "upstream" else [s("up"), s(f"c{j}"), s(f"a{j}")]

    def pr(j):
        mid = [s(f"c{j}", -1), s(f"a{j}", -1)] if actuator == "upstream" else [s(f"a{j}", -1), s(f"c{j}", -1)]
        return [s("pd"), *mid, s("px")]

    beams = [pt.beam_spec("sq0", sq(0), "g0"), pt.beam_spec("sqk", sq("k"), "gk"),
             pt.beam_spec("pr0", pr(0), "d0"), pt.beam_spec("prk", pr("k"), "dk"), pt.beam_spec("ref", [s("rf")], "d0")]
    coef = pt.coef_spec("a0", "stretcher", "*", 1.0) + pt.coef_spec("ak", "stretcher", "*", 1.0)
    return pt.topo_spec(segs, beams, ["g0", "gk", "d0", "dk"], coef)


def replay(Ta, Tpd, Tpx, Tref, t0, tk, actuator, n=700, n0=40):
    tp = topo(Ta, Tpd, Tpx, Tref, t0, tk, actuator)
    beams = ["sq0", "sqk", "pr0", "prk", "ref"]
    pipe = [pt.pathtest("pathtest_table", "act0", "drv.a0.stretcher", (), "rad", n_rows=n),
            pt.pathtest("pathtest_table", "actk", "drv.ak.stretcher", (), "rad", n_rows=n),
            pt.pathtest("pathtest_table", "las", "laser.las.phi", (), "rad", n_rows=n),
            *pt.path_pipeline(tp, {"las": "laser.las.phi"}), pt.record([f"beam.{b}.phi" for b in beams], n)]
    las = np.zeros(n)
    las[n0] = 1.0
    free, _, _ = pt.run_and_record(pipe, FS, n, {"laser.las.phi_table": las}, compiled=True)
    tabs = {"laser.las.phi_table": las}
    for j, tj in (("0", t0), ("k", tk)):
        sh = Tpx + (0 if actuator == "upstream" else tj)  # the probe sees a(t − T_px − τ_xa): anticipate
        want = free["beam.ref.phi"] - free[f"beam.pr{j}.phi"]
        a = np.zeros(n)
        a[: n - sh] = want[sh:]
        tabs[f"drv.a{j}.stretcher_table"] = a
    locked, _, _ = pt.run_and_record(pipe, FS, n, tabs, compiled=True)
    assert np.allclose(locked["beam.prk.phi"][:n - 2 * tk - 10], locked["beam.ref.phi"][:n - 2 * tk - 10], atol=1e-12)
    return locked["beam.sqk.phi"] - locked["beam.sq0.phi"]


if __name__ == "__main__":
    f = np.array([1e5, 1e6, 3e6, 1e7, 2.3e7, 4.1e7])
    print("A. general formula vs F6 replay (delays in ns at f_s = 1 GHz)")
    for d in [(30, 20, 7, 47, 5, 36), (12, 9, 3, 40, 3, 31), (25, 25, 0, 25, 4, 60)]:
        Ta, Tpd, Tpx, Tref, t0, tk = d
        for act in ("upstream", "downstream"):
            h = dtft(replay(*d, act), f) / e(f, 40)  # remove the impulse time n0 = 40
            xa0, xak = (0, 0) if act == "upstream" else (t0, tk)
            want = H_general(f, Ta, Tpd, Tpx, Tref, t0, tk, xa0, xak)
            err = np.max(np.abs(h - want))
            print(f"  {d} {act:10s}: max |H_F6 − H_general| = {err:.1e}")
            assert err < 1e-9
            if Tpx == 0 and act == "upstream":
                assert np.max(np.abs(h - H_b8(f, Ta, Tpd, Tref, t0, tk))) < 1e-9  # = B8's H
            if Tpx > 0 or act == "downstream":  # B8's H with T_b = T_pd + T_px is NOT the transfer then
                assert np.max(np.abs(want - H_b8(f, Ta, Tpd + Tpx, Tref, t0, tk))) > 1e-3

    print("B. T10 eq. (13) vs the path model (upstream actuator, T_px = 0)")
    tk_list = np.array([0, 1, 5, 10, 15, 20, 25, 31]) * 10  # τ_k − τ_0 = k·10 ns (100 MHz, N_b = 1)
    f1 = 1e4
    Ta = Tb = 20
    for Tref in (20, 60):
        h = np.array([abs(H_b8(f1, Ta, Tb, Tref, 0, tk)) for tk in tk_list])
        t10 = 4 * np.abs(np.sin(np.pi * f1 * (Ta - Tb) / FS) * np.sin(np.pi * f1 * tk_list / FS))
        print(f"  T_a = T_b = 20 ns, T_ref = {Tref} ns: |H|(10 kHz) for k = 0,1,5,…,31: "
              + " ".join(f"{x:.2e}" for x in h) + f"   T10 eq.(13): max {t10.max():.1e}")
        assert t10.max() == 0.0 and h[1:].min() > 0  # T10 says 0; the path model does not
    # leading order |H| ≈ ω²|τ_k − τ_0|·|T_a + T_ref − 2T_b − τ_k − τ_0| (any port-independent delays)
    rng = np.random.default_rng(1)
    for _ in range(5):
        Ta, Tb, Tref, t0 = rng.integers(0, 80, 4)
        tk = t0 + rng.integers(1, 300)
        w = 2 * np.pi * 1e3 / FS
        lead = w**2 * (tk - t0) * abs(Ta + Tref - 2 * Tb - tk - t0)
        assert abs(abs(H_b8(1e3, Ta, Tb, Tref, t0, tk)) - lead) <= 1e-3 * lead + 1e-15
    # B8's H vanishes iff {T_a, T_ref} = {T_b + τ_k, T_b + τ_0}: one port k only
    Tb, t0 = 15, 0
    Ta, Tref = Tb + 150, Tb + t0  # nulls port k = 15 (τ_15 = 150 ns)
    vals = {k: float(np.max(np.abs(H_b8(f, Ta, Tb, Tref, t0, 10 * k)))) for k in (1, 15, 16, 31)}
    print("  T_a = T_b + τ_15, T_ref = T_b + τ_0: max|H| over f for k = 1, 15, 16, 31:", {k: f"{v:.1e}" for k, v in vals.items()})
    assert vals[15] < 1e-12 and min(vals[1], vals[16], vals[31]) > 1e-3
    # best port-independent choice: Y = T_a + T_ref − 2T_b − 2τ_0 minimises max_k |u(Y − u)|, u = τ_k − τ_0 ∈ [0, U]
    U = 310.0
    Ys = np.linspace(0, 2 * U, 20001)
    u = np.linspace(0, U, 3101)
    worst = np.array([np.max(np.abs(u * (Y - u))) for Y in Ys])
    Yopt = Ys[np.argmin(worst)]
    print(f"  min over Y of max_k |u(Y − u)|: Y = {Yopt / U:.4f}·U, worst = {worst.min() / U**2:.4f}·U² "
          f"(Y = 0: {worst[0] / U**2:.2f}·U²)")
    assert abs(Yopt / U - 2 * (math.sqrt(2) - 1)) < 1e-3 and abs(worst.min() / U**2 - (3 - 2 * math.sqrt(2))) < 1e-3

    print("C. finite loop gain vs the ideal-lock term (integrator loop, unity gain f_u, loop delay τ_loop)")
    U_s, Y_s = 310e-9, 2 * 310e-9  # worst |Y| ~ 2U for unmatched delays
    for f_u in (1e3, 1e4, 1e5):
        for fr in (0.01, 0.1, 0.3):
            fq = fr * f_u
            w = 2 * np.pi * fq
            L = f_u / (1j * fq) * np.exp(-1j * w * 1e-6)
            leak = abs(w * U_s / (1 + L))  # (P_k − P_0)/(1 + L) ≈ −iω(τ_k − τ_0)/(1 + L)
            ideal = w**2 * U_s * Y_s
            print(f"  f_u = {f_u:.0e} Hz, f = {fq:.0e} Hz: |leak| = {leak:.2e}, |H_ideal| = {ideal:.2e}, ratio {leak / ideal:.0f}")
            if fr <= 0.1:  # well inside the band: amplitude ratio → 1/(2π f_u Y)
                assert abs(leak / ideal * 2 * np.pi * f_u * Y_s - 1) < 0.1
    # exact check of the decomposition with the F6-free algebra at one point
    Ta, Tpd, Tpx, Tref, t0, tk = 30, 20, 7, 47, 5, 36
    fq = 2e6
    for act in ("upstream", "downstream"):
        xa = {"0": 0 if act == "upstream" else t0, "k": 0 if act == "upstream" else tk}
        tot = 0
        Lg = 3.0 * np.exp(-0.4j)  # same loop gain on both ports
        for j, tj, sgn in (("k", tk, 1), ("0", t0, -1)):
            S = e(fq, Ta + tj)
            P = e(fq, Tpd + tj + Tpx) - e(fq, Tref)
            Ap, As = e(fq, Tpx + xa[j]), e(fq, tj - xa[j])
            K = Lg / Ap
            tot += sgn * (S - As * K * P / (1 + K * Ap))
        ideal = H_general(fq, Ta, Tpd, Tpx, Tref, t0, tk, xa["0"], xa["k"])
        leakx = 0
        for j, tj, sgn in (("k", tk, 1), ("0", t0, -1)):
            P = e(fq, Tpd + tj + Tpx) - e(fq, Tref)
            leakx += sgn * e(fq, tj - xa[j]) / e(fq, Tpx + xa[j]) * P / (1 + Lg)
        assert abs(tot - (ideal + leakx)) < 1e-12

    print("D. linewidth requirement, white FM (S_φ = h0/f², Δν = π h0): σ² = 2π(τ_k − τ_0)Δν, in-band part ≈ 2(τ_k − τ_0)B")
    for L31, U_s in (("63 m (100 MHz)", 310e-9), ("630 m (10 MHz)", 3.1e-6)):
        dnu = 0.01**2 / (2 * np.pi * U_s)
        print(f"  {L31}: σ_Θ ≤ 10 mrad needs Δν ≤ {dnu:.1f} Hz; fraction below B = 10 kHz: {2 * U_s * 1e4:.3f}")
    # numeric check of ∫ 4 sin²(πfU) h0/f² df = 2π² U h0
    U_s = 310e-9
    ff = np.linspace(1e-3, 200 / U_s, 4_000_001)
    val = np.trapezoid(4 * np.sin(np.pi * ff * U_s) ** 2 / ff**2, ff)
    assert abs(val / (2 * np.pi**2 * U_s) - 1) < 0.01
    print("all checks passed")
