"""M4 review (WP-B8 demux): independent checks of ``demux_design`` / ``delay_line_k`` / ``demux_switch``.

  1  time-bin bookkeeping of a binary tree (per-depth shared drives) with N_b, n_guard: every slot is
     propagated stage by stage with explicit emission and arrival times; (a) the power map equals
     ``transfer_matrix``; (b) no leaked light ever lands in a target bin (crosstalk = pure loss for gated
     detection, T11 §3.3); (c) guard-slot photons land in non-target bins but are not in ``leak_in``.
  2  edge window between blocks = (1 + n_guard)·T_rep (D62-4) from the explicit slot timeline.
  3  sign of the rep-rate term of δt_k: explicit emission times vs ``delay_line_k``.
  4  shared stages: ports 0 and k pass the same physical stage-1 … switches; their common-phase bias
     drift is common to both, so Θ_k = θ_k − θ_0 cancels it. Two ``demux_switch`` instances draw
     independent walks (corr ≈ 0).

Run: uv run python reviews/scripts/m4_demux_checks.py
"""
import math

import numpy as np

import cranberry.blocks.dmx.testing as dt
from cranberry.blocks.dmx.demux import transfer_matrix
from cranberry.core.compile import run_reference


def tree_timebins(K, N_b, n_guard, t, eps, n_frames=3):
    """Propagate every slot of n_frames frames through a binary tree whose depth-s switches all follow
    bit s of the port addressed at that slot (guard slots: switch state of the next block = worst case).
    Returns {(port, arrival slot): power}, target bins, and per-slot origin flags."""
    S = int(math.log2(K))
    per = N_b + n_guard
    frame = K * per
    out = {}
    for j in range(n_frames * frame):
        blk, pos = divmod(j % frame, per)
        guard = pos >= N_b
        target = K - 1 - blk  # σ(k) = K − 1 − k
        if guard:
            target = K - 1 - ((blk + 1) % K)  # edge already moved to the next block
        bits = [(target >> (S - 1 - s)) & 1 for s in range(S)]
        # amplitude-free power bookkeeping: each stage passes the commanded way with t(1−ε), wrong way t·ε
        paths = {(): 1.0}
        for s in range(S):
            new = {}
            for p, w in paths.items():
                new[p + (bits[s],)] = new.get(p + (bits[s],), 0) + w * t * (1 - eps)
                new[p + (1 - bits[s],)] = new.get(p + (1 - bits[s],), 0) + w * t * eps
            paths = new
        for p, w in paths.items():
            port = int("".join(map(str, p)), 2)
            arr = j + port * per  # delay τ_k = k·(N_b + n_guard)·T_rep (L_0 = 0)
            key = (port, arr, "guard" if guard else "data")
            out[key] = out.get(key, 0.0) + w
    tgt = {}
    for j in range(n_frames * frame):
        blk, pos = divmod(j % frame, per)
        if pos < N_b:
            k = K - 1 - blk
            tgt.setdefault(k, set()).add(j + k * per)
    return out, tgt, frame, per


if __name__ == "__main__":
    K, t, eps = 8, 0.99, 0.01
    for N_b, n_guard in ((1, 0), (3, 1), (2, 2)):
        out, tgt, frame, per = tree_timebins(K, N_b, n_guard, t, eps)
        T = transfer_matrix(K, "tree", t, eps)
        # (a) power from one data slot targeted at port k into port k′ (middle frame)
        j0 = frame  # first slot of frame 1, targeted at port K − 1
        got = {kp: sum(w for (p, a, kind), w in out.items() if p == kp and a == j0 + kp * per and kind == "data") for kp in range(K)}
        np.testing.assert_allclose([got[kp] for kp in range(K)], T[:, K - 1], rtol=1e-12)
        # (b) target bins (middle frame) contain only their own data slot
        for k in range(K):
            for a in tgt[k]:
                if not frame <= a - k * per < 2 * frame:
                    continue
                contrib = [(kind, w) for (p, aa, kind), w in out.items() if p == k and aa == a]
                assert len(contrib) == 1 and contrib[0][0] == "data" and abs(contrib[0][1] - T[k, k]) < 1e-12
        # (c) photons (per unit sinh²r) per frame into non-target bins: data leak vs guard slots
        leak_data = sum(w for (p, a, kind), w in out.items() if kind == "data" and frame <= a - p * per < 2 * frame
                        and a not in tgt[p])
        leak_guard = sum(w for (p, a, kind), w in out.items() if kind == "guard" and frame <= a - p * per < 2 * frame)
        leak_in = (T.sum(axis=1) - np.diag(T)).sum() * N_b
        print(f"N_b = {N_b}, n_guard = {n_guard}: target bins clean; non-target photons/frame (×sinh²r): "
              f"data leak {leak_data:.4f} (Σ leak_in·N_b = {leak_in:.4f}), guard slots {leak_guard:.3f}")
        assert abs(leak_data - leak_in) < 1e-12
        if n_guard:
            assert leak_guard > 10 * leak_data  # guard pulses dominate the non-target-bin light

    print("2. gap between the last pulse of a block and the first of the next (in T_rep)")
    for N_b, n_guard in ((1, 0), (4, 0), (4, 2), (50, 3)):
        per = N_b + n_guard
        data = [j for j in range(3 * per) if j % per < N_b]
        gaps = sorted({b - a for a, b in zip(data[:-1], data[1:]) if a // per != b // per})
        print(f"  N_b = {N_b}, n_guard = {n_guard}: inter-block gap {gaps} T_rep (T11 §3.4 used N_b·T_rep = {N_b})")
        assert gaps == [1 + n_guard]

    print("3. rep-rate drift: arrival of port k minus port 0 (L_k fixed for T_rep)")
    T_rep, rel, k, Nb = 10e-9, 1e-6, 31, 1
    tau = {kk: kk * Nb * T_rep for kk in (0, k)}  # design delays
    T2 = T_rep * (1 + rel)
    emit = {kk: (31 - kk) * Nb * T2 for kk in (0, k)}  # slot σ(k) = 31 − k of a 32-port frame
    true = (emit[k] + tau[k]) - (emit[0] + tau[0])
    shift = {}
    for kk in (0, k):
        b, st, pr = dt.make(dt.delay(port=kk, rep_rel_drift=rel, toggles={"static_mismatch": False}), 1e6)
        run_reference(b, st, pr, 1)
        shift[kk] = float(st["timing.dl.dt"])
    code = shift[k] - shift[0]
    print(f"  true δt_31 − δt_0 = {true * 1e15:+.1f} fs (arrives early), delay_line_k: {code * 1e15:+.1f} fs")
    assert abs(true + 310e-15) < 1e-18 and abs(code - 310e-15) < 1e-18  # opposite signs
    b, st, pr = dt.make(dt.delay(port=k, in_dT=None, toggles={"static_mismatch": True, "rep_rate_drift": False},
                                 dtau_s=1e-12), 1e6)
    run_reference(b, st, pr, 1)
    assert float(st["timing.dl.dt"]) == 1e-12  # static term: + = longer delay = late

    print("4. two demux_switch instances (ports 0 and k) draw independent bias walks")
    from cranberry.blocks.err.testing import record
    sw = lambda name: dt.switch(name=name, n_stages=5, toggles={"drive_phase_noise": False, "bias_drift": True},  # noqa: E731
                                bias_rw_rad_per_rts=1.0, driver="eom", segment=f"dmx_{name}")
    rec = record([sw("p0"), sw("pk")], 1e4, 20000, ["act.p0.dpsi", "act.pk.dpsi"], external=())
    a, b2 = np.diff(rec["act.p0.dpsi"]), np.diff(rec["act.pk.dpsi"])
    c = float(np.corrcoef(a, b2)[0, 1])
    var_ind = np.var(np.diff(rec["act.pk.dpsi"] - rec["act.p0.dpsi"]))
    print(f"  corr of increments {c:+.3f}; Var of the Θ_k increment from the switches, independent model: {var_ind:.3e}")
    assert abs(c) < 0.05 and abs(var_ind / 1e-3 - 1) < 0.05  # 2·S·κ²σ²h = 10·1e-4
    S = 5
    for kk in (1, 2, 3, 8, 16, 31):
        m = S - int(math.log2(kk))  # highest differing bit (1 = MSB): ports share depths 1..m
        print(f"  port {kk:2d} vs 0 in a 32-port tree: shared stages {m}, physical Var = {2 * (S - m)}·κ²σ²h "
              f"(independent blocks: {2 * S}·κ²σ²h)")
    print("all checks passed")
