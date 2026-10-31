"""M6 review of the 32-port scale estimate (D74-6, future_support.scale_estimate).

The estimate apportions the 2-port loop time to a port by max(leaf share, state share) and extrapolates linearly.
Here a THIRD port is added for real (an overlay generated from F_common.yaml's port-1 text: segments DMX/D/F/G/PRS/
REF/LO, beams sqz/pr/ref/lo, lock, port entry, coefficients, the whole L-Liu lock stack; the fibre noise / thermal
blocks get the new segment; demux_switch / delay_line_k of the new port are left out — two light blocks), and the
measured increment of µs/sample from 2 to 3 ports is compared with the estimate's per-port cost.
Usage: uv run python reviews/scripts/m6_scale_three_ports.py
"""
import importlib.util
import re
import time
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
_sp = importlib.util.spec_from_file_location("future_support", ROOT / "tests/integration/future_support.py")
fs = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(fs)
from cranberry import config  # noqa: E402
from cranberry.runner.run import run  # noqa: E402

SRC = (ROOT / "presets/future/F_common.yaml").read_text().splitlines()
SEG = r"(DMX|D|F|G|PRS|REF|LO|sqz|pr|ref|lo|det|L|liu|adc|pid|sup|dac|dl|pzt|tec|stage|act|lock|spool)1(?![0-9])"


def two(line: str) -> str:
    return re.sub(SEG, r"\g<1>2", line).replace("{segment: EM1, direction: 1}, ", "")


def grab(pred, lo=0, hi=None):
    out, i, lines = [], lo, SRC[lo:hi]
    while i < len(lines):
        if pred(lines[i]):
            blk, ind0 = [lines[i]], len(lines[i]) - len(lines[i].lstrip())
            while i + 1 < len(lines) and lines[i + 1].startswith(" " * (ind0 + 2)):
                i += 1
                blk.append(lines[i])
            out += blk
        i += 1
    return out


segs = grab(lambda s: re.match(r"\s+- \{id: (DMX1|D1|F1|G1|PRS1|REF1|LO1),", s) is not None)
beams = grab(lambda s: re.match(r"\s+- \{id: (sqz1|pr1|ref1|lo1),", s) is not None)
coef = [s for s in SRC if re.match(r"\s+- \{id: (D1|F1|G1|PRS1|REF1|LO1)\.", s)]   # DMX2: no switch block (header)
i0 = next(i for i, s in enumerate(SRC) if "L-Liu, port k" in s)
i1 = next(i for i, s in enumerate(SRC) if "---- metrics" in s)
stack = SRC[i0 + 1:i1]
ind = lambda ls, n: [" " * n + s.lstrip() if s.lstrip().startswith("- ") else s for s in ls]  # noqa: E731
segs = [x.replace("drivers: {eom: uniform}", "drivers: {}") for x in segs]
y = ["offline:", "  params:", "    L_lo2_m: !expr \"L_lo1_m\"", "topology:", "  segments:", *map(two, segs), "  beams:", *map(two, beams),
     "  observers:", "    - {id: det2, kind: detector}",
     "  locks:", "    - {id: L2, combo: {pr2: 1.0, ref2: -1.0}, setpoint_rad: !pi_frac 0.5, observer: det2}",
     "  ports:",   # not a keyed list (no id/name): the overlay replaces it, so all three entries are written
     "    - {index: 0, squeezed: sqz0, lo: lo0, probe: pr0, gbs_observer: gbs, bhd_observer: gbs}",
     "    - {index: 1, squeezed: sqz1, lo: lo1, probe: pr1, gbs_observer: gbs, bhd_observer: gbs}",
     "    - {index: 2, squeezed: sqz2, lo: lo2, probe: pr2, gbs_observer: gbs, bhd_observer: gbs}",
     "  coef:", *map(two, coef),
     "pipeline:",
     "  - name: env.thermal.F1",
     "    params: {segments: [F1, F2], tau_th_s: {F1: !expr \"offline.params.tau_th_spool_s\", F2: !expr \"offline.params.tau_th_spool_s\"}}",
     "  - name: env.thermal.LO",
     "    params: {segments: [LO0, LO1, LO2], tau_th_s: {LO0: !expr \"offline.params.tau_th_spool_s\", LO1: !expr \"offline.params.tau_th_spool_s\", LO2: !expr \"offline.params.tau_th_spool_s\"}}",
     "  - name: env.mech.lo",
     "    params: {segments: [LO0, LO1, LO2]}",
     "  - name: env.mech.dist",
     "    params: {segments: [PRsrc, PRS0, PRS1, PRS2, REF, REF0, REF1, REF2]}",
     "  - name: env.mech.fs_lock",
     "    params: {segments: [D0, D1, D2]}",
     "  - name: env.mech.fs_sq",
     "    params: {segments: [PUMP, UP, G0, G1, G2]}",
     "  - name: env.mech.ports",
     "    params: {segments: [F0, F1, F2], spatial_coherence: {F0: 0.5, F1: 0.5, F2: 0.5}}",
     *map(two, stack)]
ov = Path(__file__).with_name("_m6_port2_overlay.yaml")
ov.write_text("\n".join(y) + "\n")


def us_per_sample(tree_or_path, overlays=()):
    walls = []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mk = lambda d: config.load({"extends": str(ROOT / fs.nominal("F_liu20")), "run": {"duration_s": d}}, overlays=list(overlays))  # noqa: E731
        run(mk(0.006))
        reps = {0.05: [], 0.3: []}
        for _ in range(3):           # interleaved repeats, min of each (shared machine: timing noise)
            for d in reps:
                t0 = time.perf_counter()
                res = run(mk(d))
                reps[d].append(time.perf_counter() - t0)
        walls = [min(reps[0.05]), min(reps[0.3])]
        print(f"  walls {({d: [round(x, 2) for x in v] for d, v in reps.items()})}")
    n = [int(d * res.spec.run.f_sim_hz) for d in (0.05, 0.3)]
    return (walls[1] - walls[0]) / (n[1] - n[0]) * 1e6, res


try:
    u2, r2 = us_per_sample(None)
    u3, r3 = us_per_sample(None, overlays=(str(ov),))
finally:
    ov.unlink(missing_ok=True)
assert set(r3.metrics.metric_ii) == {1, 2}, r3.metrics.metric_ii.keys()
est = fs.scale_estimate("F_liu20")
inc = u3 - u2
print(f"2 ports: {u2:.3f} µs/sample, 3 ports: {u3:.3f} µs/sample -> measured per-port increment {inc:.3f} µs/sample")
print(f"estimate: {est['us_per_sample']:.3f} µs/sample (2 ports), per port {est['us_per_sample_per_port']:.3f}, "
      f"32 ports {est['target_us_per_sample']:.2f} µs/sample")
print(f"32 ports from the measured increment: {u2 + 30 * inc:.2f} µs/sample "
      f"({(u2 + 30 * inc) * 1e-6 * 1e6:.1f} s wall per simulated s at 1 MHz)")
print(f"port-2 metric (ii) std {r3.metrics.metric_ii[2].stats.std:.4f} rad (sanity: port 2 duplicates port k)")
