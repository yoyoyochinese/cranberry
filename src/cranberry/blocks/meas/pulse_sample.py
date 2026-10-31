"""``pulse_sample``: per-pulse BHD sample variance (G16-f; T13 §1.8, T11 §3.7, D45). Owner: WP-B5.

Measurement stage. A frame has K slots (pulse bins) with normalized quadrature variances R_k
(shot-noise units, measured at the BHD's LO phase); the train repeats the frame. Slot charges
q_j ~ N(0, R_j) mix through the detector coefficients c_m (``pulsed_bhd_kernel``): Y_j = Σ_m c_m q_{j−m}.
Electronic noise of variance 1/C is added per estimate, in units of one pulse's own shot noise
(c_0 = 1; C = C_pulse = η_q N_LO/σ_el²; for ``integrating`` per frame-normalized estimate).

``mode``:

* ``pulse_resolved`` — LO in every slot, one estimate per target slot:
  R_meas = (Σ_m c_m² R_{j−m} + 1/C) / (Σ_m c_m² + 1/C).
* ``gated`` — LO only in the target slots (no homodyne charge elsewhere):
  same with R_{j−m} → R_{j−m}·[j − m is a target slot] in both numerator and reference.
* ``integrating`` — LO over the whole train, the detector integrates a frame (charge conserved, so
  c_m drops out): R_meas = (mean_k R_k + 1/C) / (1 + 1/C). T11: K = 32, one slot at −10 dB, the rest
  vacuum → −0.124 dB (D45 dilution).

Normalization to shot noise uses the same estimator on a vacuum train (as the experiment does).
:meth:`PulseSample.evaluate` returns the analytic expectation and, with toggle ``sampling_noise``,
a Monte-Carlo estimate from ``n_frames`` frames (sample variances; adjacent-estimate correlation).

Toggles: ``crosstalk`` (off ⇒ c = [1]), ``electronic_noise`` (off ⇒ 1/C = 0), ``sampling_noise``.
Params (all required, D34): mode, k_slots, target_slots, n_frames, toggles.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Any

import numba
import numpy as np

from cranberry.core.block import BaseBlock, Block, require
from cranberry.core.registry import register_block
from cranberry.core.units import lin_to_db

__all__ = ["PulseSample", "make_pulse_sample", "PULSE_MODES", "expected_r_meas"]

PULSE_MODES: tuple[str, ...] = ("pulse_resolved", "gated", "integrating")


def _lo_mask(mode: str, k: int, targets: Sequence[int]) -> np.ndarray:
    m = np.ones(k) if mode != "gated" else np.zeros(k)
    m[list(targets)] = 1.0
    return m


def expected_r_meas(r_slots: Any, mode: str, targets: Sequence[int], c: Any, inv_c: float) -> np.ndarray:
    """Analytic R_meas per target slot (module docstring); ``inv_c`` = 1/C (0 for no electronic noise)."""
    r = np.asarray(r_slots, dtype=float)
    k = r.size
    if mode == "integrating":
        return np.full(len(targets), (r.mean() + inv_c) / (1.0 + inv_c))
    lo = _lo_mask(mode, k, targets)
    c2 = np.asarray(c, dtype=float) ** 2
    out = []
    for j in targets:
        idx = (j - np.arange(c2.size)) % k
        num = np.dot(c2, r[idx] * lo[idx]) + inv_c
        out.append(num / (np.dot(c2, lo[idx]) + inv_c))
    return np.array(out)


@numba.njit
def _noop(S, Z, I, P):  # noqa: N803
    pass


class PulseSample(BaseBlock):
    """``pulse_sample`` (module docstring)."""

    KIND = "pulse_sample"
    toggle_names = ("crosstalk", "electronic_noise", "sampling_noise")

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or self.KIND
        require(params, "mode", "k_slots", "target_slots", "n_frames", "toggles", block=name)
        super().__init__(self.KIND, name, params)
        if params["mode"] not in PULSE_MODES:
            raise ValueError(f"block {name!r}: mode must be one of {PULSE_MODES}")
        self.k = int(params["k_slots"])
        self.targets = [int(j) for j in params["target_slots"]]
        if not self.targets or any(not 0 <= j < self.k for j in self.targets):
            raise ValueError(f"block {name!r}: target_slots must be non-empty and in [0, k_slots)")

    def bind(self, layout) -> None:
        self.kernel = _noop

    def __call__(self, state, params):
        return state

    def _estimates(self, q: np.ndarray, c: np.ndarray, el: np.ndarray, lo: np.ndarray) -> np.ndarray:
        """Estimates per (frame, target) from slot charges q (n_frames·K, train order)."""
        mode = self.params["mode"]
        if mode == "integrating":
            return (q.reshape(-1, self.k).sum(axis=1) / math.sqrt(self.k))[:, None] + el[:, :1]
        qq = q * np.tile(lo, q.size // self.k)
        y = np.convolve(qq, c)[: q.size]  # causal: Y_j = Σ_m c_m q_{j−m}
        return y.reshape(-1, self.k)[:, self.targets] + el

    def evaluate(self, r_slots: Any, *, c: Any, clearance_lin: float, rng: np.random.Generator | None = None) -> dict:
        """R_meas per target slot for slot variances ``r_slots`` (K values, shot-noise units).

        ``c``, ``clearance_lin``: from ``pulsed_bhd_kernel`` (``c``, ``clearance_pulse``). Returns
        ``expected`` / ``expected_db`` and, with ``sampling_noise``, ``sample`` / ``sample_db``,
        ``rho_adj_sample`` (correlation of successive estimates in the vacuum reference run) and
        ``n_samples``.
        """
        r = np.asarray(r_slots, dtype=float)
        if r.shape != (self.k,):
            raise ValueError(f"r_slots must have shape ({self.k},)")
        t = self._toggles
        cc = np.asarray(c, dtype=float) if t["crosstalk"] else np.ones(1)
        inv_c = (1.0 / clearance_lin if math.isfinite(clearance_lin) else 0.0) if t["electronic_noise"] else 0.0
        exp = expected_r_meas(r, self.params["mode"], self.targets, cc, inv_c)
        out: dict[str, Any] = {"mode": self.params["mode"], "targets": list(self.targets), "expected": exp,
                               "expected_db": lin_to_db(exp)}
        if not t["sampling_noise"]:
            return out
        if rng is None:
            raise ValueError(f"block {self.name!r}: sampling_noise needs an rng (core.rng.stream)")
        nf = int(self.params["n_frames"])
        lo = _lo_mask(self.params["mode"], self.k, self.targets)

        def run(rv: np.ndarray) -> np.ndarray:  # σ_el² = 1/C in units of one pulse's shot noise (c_0 = 1)
            q = rng.standard_normal(nf * self.k) * np.sqrt(np.tile(rv, nf))
            el = rng.standard_normal((nf, len(self.targets))) * math.sqrt(inv_c)
            return self._estimates(q, cc, el, lo)

        x, x0 = run(r), run(np.ones(self.k))
        sample = x.var(axis=0) / x0.var(axis=0)
        flat = x0.ravel()
        out.update({"sample": sample, "sample_db": lin_to_db(sample), "n_samples": nf,
                    "rho_adj_sample": float(np.corrcoef(flat[:-1], flat[1:])[0, 1])})
        return out


@register_block(PulseSample.KIND)
def make_pulse_sample(params: Mapping[str, Any], *, name: str | None = None) -> Block:
    """Factory of ``pulse_sample`` (module docstring)."""
    return PulseSample(params, name)
