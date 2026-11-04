"""``loop_identification``: open-loop gain L̂(f), S, T and margins from the recorded signals of one locked run
with an ``excitation`` (Phase 6 P7, D85-7; notes/lit-kashiwazaki2023.md §5.3 N7-2). Measurement-stage hook
``measure.loop_identification`` (runner, after the metrics); the time-loop kernel is a no-op.

Signals (slot names, params): ``x`` = ``ctrl.<excitation>.x`` (the instrument), ``reference`` and ``response``.
The runner adds them to ``measure.record`` itself. All ratios are formed against x (instrumental variable), so
loop noise uncorrelated with x does not bias them; it lowers the coherence.

Relations (``relation``):

* ``minus_l`` — injection at a summing point inside the loop (P7 (B), the Moku:Delta form): ``reference`` =
  the signal after the sum (what the controller sees, e + x), ``response`` = the signal before it (e).
  With e = −L(e + x): L̂ = −resp/ref, S = ref/x = 1/(1 + L), T = −resp/x = L/(1 + L). Every stage driven
  from that error is inside L (a stacked loop gives the full stack).
* ``plant`` — injection at the input of an actuator whose controller is replaced by the analyser while
  another stage holds the lock (P7 literal, (A)): ``reference`` null (x is the reference), ``response`` the
  error. P̂ = resp/x is the plant seen from the injection point (it carries 1/(1 + L_hold) of the holding
  loop, ≈ 1 above that loop's band); L̂ = −P̂·C(f) with ``controller`` = {kp, ki_per_s, kd_s,
  discretisation?, h_s?} (the ``loop_stability_check`` controller form) or ``model`` (C₁ of the model block).
  S = 1/(1 + L̂), T = L̂/(1 + L̂) are then computed, not measured.

Estimators (by the excitation's waveform):

* stepped sine: per step, single-bin lock-in over the samples of [t_k + settle_s, t_k + dwell_s) (an integer
  number of periods: the excitation snaps f_k; exact when the window is a multiple of the record step), H₁
  against x. Coherence from ``n_blocks`` contiguous sub-blocks of the window:
  γ² = |Σ R_m X_m*|² / (Σ|R_m|² Σ|X_m|²).
* multisine / chirp / prbs: Welch cross-spectra with a rectangular window of one excitation period, no overlap,
  from the first period boundary after t_start + settle_s (leakage-free for the periodic part); multisine at its
  tones, chirp / prbs at the bins in ``band_hz`` where S_xx ≥ 1e−6 max S_xx. H = S_{x,y}/S_{x,x}, γ²_{x,y} =
  |S_{x,y}|²/(S_xx S_yy).
* white: the same with a Hann window of ``nperseg_s`` and 50 % overlap, bins in ``band_hz``.

For ``minus_l`` the reported coherence is min(γ²_{x,ref}, γ²_{x,resp}); for ``plant`` γ²_{x,resp}.

Margins (:func:`loop_margins`, on the points with coherence ≥ ``coherence_min``): phase unwrapped along f and
put on the branch with the lowest-frequency phase in (−360°, 0°] (or nearest the model there, when a model is
given); f_u = first downward crossing of |L̂| = 1 (log–log interpolation), phase margin 180° + ∠L̂(f_u); f_180 =
first crossing of −180°, gain margin −20 log₁₀|L̂(f_180)|. NaN when the band holds no crossing.

Model (``model``: a ``loop_stability_check`` block name, or null): L_model = ``open_loop`` of its params on
the same frequencies (the D82 ``dsc1`` discretisation when the block has it); ``model_stages`` ``all`` (the full
stack, for ``minus_l``) or ``first`` (stages dropped: the fast path alone, for ``plant``). Reported as
|L̂|/|L_model| − 1 and ∠(L̂/L_model) [deg], and the model's own margins on a dense grid over the measured band.

Params: ``excitation``, ``x``, ``response``, ``relation``, ``model_stages``, ``coherence_min``; ``reference``
(minus_l), ``controller`` (plant), ``n_blocks`` (stepped sine), ``band_hz`` (chirp, prbs, white), ``nperseg_s``
(white). ``model`` and, for ``plant``, ``reference`` may be absent or null (YAML overlays delete a key set to
null, so absence is the overlay spelling of "none").
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from typing import Any

import numba
import numpy as np
import scipy.signal

from cranberry.blocks.chk.loop_stability import open_loop
from cranberry.blocks.ctrl.excitation import excitation_schedule
from cranberry.core.block import BaseBlock, Block, require
from cranberry.core.registry import register_block

__all__ = ["LoopIdentification", "make_loop_identification", "identify", "loop_margins", "RELATIONS"]

KIND = "loop_identification"
RELATIONS = ("minus_l", "plant")
MODEL_STAGES = ("all", "first")


def _demod(y: np.ndarray, t: np.ndarray, f: float) -> complex:
    return complex(np.sum(y * np.exp(-2j * np.pi * f * t)))


def _gamma2(a: np.ndarray, b: np.ndarray) -> float:
    """|Σ a b*|² / (Σ|a|² Σ|b|²) over sub-blocks (a: output, b: instrument)."""
    den = float(np.sum(np.abs(a) ** 2) * np.sum(np.abs(b) ** 2))
    return float(abs(np.sum(a * np.conj(b))) ** 2 / den) if den > 0 else math.nan


def _lockin(t, n_rec, sig, sc, f_sim, n_blocks):
    n0 = int(math.ceil(sc["t_start_s"] * f_sim - 1e-6))
    dwell_n = int(round(sc["dwell_s"] * f_sim))
    settle_n = int(round(sc["settle_s"] * f_sim))
    freqs = np.asarray(sc["freqs_hz"], float)
    h = {k: np.full(freqs.size, np.nan + 0j) for k in sig if k != "x"}
    coh = {k: np.full(freqs.size, np.nan) for k in sig if k != "x"}
    for k, f in enumerate(freqs):
        lo, hi = n0 + k * dwell_n + settle_n, n0 + (k + 1) * dwell_n
        sel = np.nonzero((n_rec >= lo) & (n_rec < hi))[0]
        if sel.size < 2:
            raise ValueError(f"loop_identification: step {k} ({f:g} Hz) has no recorded samples "
                             f"(record shorter than the excitation schedule?)")
        tt = t[sel]
        X = _demod(sig["x"][sel], tt, f)
        # sub-blocks on whole periods of f (the window holds m of them): leakage-free coherence
        m = max(1, int(round(f * (hi - lo) / f_sim)))
        nb = max(1, min(int(n_blocks), m))
        cuts = np.searchsorted(n_rec[sel], lo + np.round(np.round(np.arange(1, nb) * m / nb) / f * f_sim))
        parts = np.split(sel, cuts)
        Xm = np.array([_demod(sig["x"][p], t[p], f) for p in parts])
        for key in h:
            h[key][k] = _demod(sig[key][sel], tt, f) / X
            coh[key][k] = _gamma2(np.array([_demod(sig[key][p], t[p], f) for p in parts]), Xm)
    return freqs, h, coh


def _spectral(t, n_rec, sig, sc, f_sim, band_hz, nperseg_s):
    f_rec = 1.0 / float(t[1] - t[0])
    n0 = int(math.ceil(sc["t_start_s"] * f_sim - 1e-6))
    settle_n = int(round(sc["settle_s"] * f_sim))
    w = sc["waveform"]
    if w == "white":
        start = n0 + settle_n
        nper, window, nover = int(round(float(nperseg_s) * f_rec)), "hann", None
    else:
        per_n = int(round(sc["period_s"] * f_sim))
        start = n0 + int(math.ceil(settle_n / per_n)) * per_n
        nper_f = sc["period_s"] * f_rec
        nper = int(round(nper_f))
        if abs(nper_f - nper) > 1e-6 * nper_f:
            raise ValueError(f"loop_identification: excitation period × record rate = {nper_f:.9g} must be an integer")
        window, nover = "boxcar", 0
    sel = n_rec >= start
    n_use = int(np.sum(sel))
    if w != "white":
        n_use -= n_use % nper
    if nper < 2 or n_use < nper:
        raise ValueError("loop_identification: record too short for one analysis segment after t_start + settle_s")
    idx = np.nonzero(sel)[0][:n_use]
    kw = dict(fs=f_rec, window=window, nperseg=nper, noverlap=nover, detrend=False, scaling="spectrum")
    xx = sig["x"][idx]
    f, sxx = scipy.signal.welch(xx, **kw)
    h, coh = {}, {}
    for key in (k for k in sig if k != "x"):
        y = sig[key][idx]
        _, sxy = scipy.signal.csd(xx, y, **kw)
        _, syy = scipy.signal.welch(y, **kw)
        h[key] = sxy / sxx
        coh[key] = np.abs(sxy) ** 2 / (sxx * syy)
    if w == "multisine":
        pick = np.round(np.asarray(sc["freqs_hz"]) * sc["period_s"]).astype(int)
    else:
        lo, hi = (float(v) for v in band_hz)
        pick = np.nonzero((f >= lo) & (f <= hi) & (f > 0) & (sxx >= 1e-6 * np.max(sxx)))[0]
    return f[pick], {k: v[pick] for k, v in h.items()}, {k: v[pick] for k, v in coh.items()}


def _branch(ph_deg: np.ndarray, ref0: float | None) -> np.ndarray:
    ph = np.degrees(np.unwrap(np.radians(ph_deg)))
    if ph.size == 0:
        return ph
    if ref0 is None:
        return ph - 360.0 * math.ceil(ph[0] / 360.0)  # ph[0] in (−360°, 0°]
    return ph - 360.0 * round((ph[0] - ref0) / 360.0)


def loop_margins(f_hz: np.ndarray, l_f: np.ndarray, phase_ref0_deg: float | None = None) -> dict[str, float]:
    """f_u, phase margin, f_180, gain margin of L on the (ascending) grid ``f_hz`` (module docstring)."""
    f = np.asarray(f_hz, float)
    l_f = np.asarray(l_f, complex)
    out = {"f_u_hz": math.nan, "phase_margin_deg": math.nan, "f_180_hz": math.nan, "gain_margin_db": math.nan}
    if f.size < 2:
        return out
    mag = np.abs(l_f)
    ph = _branch(np.degrees(np.angle(l_f)), phase_ref0_deg)
    lf, lm = np.log(f), np.log(mag)
    i = np.nonzero((mag[:-1] >= 1.0) & (mag[1:] < 1.0))[0]
    if i.size:
        i = int(i[0])
        a = -lm[i] / (lm[i + 1] - lm[i])
        out["f_u_hz"] = float(math.exp(lf[i] + a * (lf[i + 1] - lf[i])))
        out["phase_margin_deg"] = float(180.0 + ph[i] + a * (ph[i + 1] - ph[i]))
    j = np.nonzero((ph[:-1] + 180.0) * (ph[1:] + 180.0) <= 0.0)[0]
    j = j[ph[j] != ph[j + 1]]
    if j.size:
        j = int(j[0])
        a = (-180.0 - ph[j]) / (ph[j + 1] - ph[j])
        out["f_180_hz"] = float(math.exp(lf[j] + a * (lf[j + 1] - lf[j])))
        out["gain_margin_db"] = float(-20.0 / math.log(10.0) * (lm[j] + a * (lm[j + 1] - lm[j])))
    return out


def _controller(spec: Any, model_params: Mapping[str, Any] | None, f: np.ndarray) -> np.ndarray:
    if spec == "model":
        if model_params is None:
            raise ValueError("loop_identification: controller 'model' needs a model block")
        spec = model_params
    q = {"k0": 1.0, "delay_s": 0.0, "poles_hz": [], "modes": [], "kp": spec["kp"], "ki_per_s": spec["ki_per_s"],
         "kd_s": spec["kd_s"]}
    if spec.get("discretisation", "continuous") == "dsc1":
        q.update(discretisation="dsc1", h_s=spec["h_s"])
    return open_loop(f, q)


def identify(series: Mapping[str, np.ndarray], t_s: np.ndarray, schedule: Mapping[str, Any], f_sim_hz: float, *,
             relation: str, controller: Any = None, model_params: Mapping[str, Any] | None = None,
             model_stages: str = "all", coherence_min: float = 0.0, n_blocks: int = 4,
             band_hz: Any = None, nperseg_s: float | None = None) -> dict[str, Any]:
    """L̂, S, T, coherence, margins and the model comparison (module docstring) from recorded series
    ``{"x", "response"[, "reference"]}`` on the common time axis ``t_s`` (samples at n/f_sim)."""
    if relation not in RELATIONS:
        raise ValueError(f"relation must be one of {RELATIONS}, got {relation!r}")
    if model_stages not in MODEL_STAGES:
        raise ValueError(f"model_stages must be one of {MODEL_STAGES}, got {model_stages!r}")
    t = np.asarray(t_s, float)
    n_rec = np.round(t * f_sim_hz).astype(np.int64)
    sig = {"x": np.asarray(series["x"], float), "response": np.asarray(series["response"], float)}
    if relation == "minus_l":
        sig["reference"] = np.asarray(series["reference"], float)
    if schedule["waveform"] == "stepped_sine":
        f, h, coh = _lockin(t, n_rec, sig, schedule, f_sim_hz, n_blocks)
    else:
        if schedule["waveform"] != "multisine" and band_hz is None:
            raise ValueError(f"loop_identification: waveform {schedule['waveform']!r} needs band_hz")
        if schedule["waveform"] == "white" and nperseg_s is None:
            raise ValueError("loop_identification: waveform 'white' needs nperseg_s")
        f, h, coh = _spectral(t, n_rec, sig, schedule, f_sim_hz, band_hz, nperseg_s)
    out: dict[str, Any] = {"relation": relation, "waveform": schedule["waveform"], "f_hz": f}
    if relation == "minus_l":
        s, tt = h["reference"], -h["response"]
        lh = tt / s
        c = np.minimum(coh["reference"], coh["response"])
        out.update(s_hat=s, t_hat=tt, coherence_ref=coh["reference"], coherence_resp=coh["response"])
    else:
        p = h["response"]
        lh = -p * _controller(controller, model_params, f)
        s, tt = 1.0 / (1.0 + lh), lh / (1.0 + lh)
        c = coh["response"]
        out.update(p_hat=p, p_mag=np.abs(p), p_phase_deg=np.degrees(np.angle(p)), s_hat=s, t_hat=tt)
    out.update(l_hat=lh, l_mag=np.abs(lh), coherence=c, s_mag=np.abs(s), t_mag=np.abs(tt))
    ok = np.isfinite(c) & (c >= float(coherence_min))
    lm = None
    ref0 = None
    if model_params is not None:
        mp = dict(model_params) if model_stages == "all" else {**model_params, "stages": []}
        lm = open_loop(f, mp)
        if ok.any():  # the model's own phase branch at the lowest used point
            ref0 = float(_branch(np.degrees(np.angle(lm[ok])), None)[0])
    ph = np.full(f.size, np.nan)
    if ok.any():
        ph[ok] = _branch(np.degrees(np.angle(lh[ok])), ref0)
    out["l_phase_deg"] = ph
    out["n_used"] = int(ok.sum())
    out.update(loop_margins(f[ok], lh[ok], ref0))
    if lm is not None:
        out.update(l_model=lm, l_model_mag=np.abs(lm), l_model_phase_deg=_branch(np.degrees(np.angle(lm)), None),
                   mag_err=np.abs(lh) / np.abs(lm) - 1.0, phase_err_deg=np.degrees(np.angle(lh / lm)))
        if f.size:
            fd = np.geomspace(float(f.min()), float(f.max()), 4000)
            mm = loop_margins(fd, open_loop(fd, mp))
            out.update({f"model_{k}": v for k, v in mm.items()})
    return out


@numba.njit
def _noop(S, Z, I, P):  # noqa: N803
    pass


class LoopIdentification(BaseBlock):
    """The ``loop_identification`` block (module docstring)."""

    def __init__(self, params: Mapping[str, Any], name: str | None):
        name = name or KIND
        require(params, "excitation", "x", "response", "relation", "model_stages", "coherence_min", block=name)
        super().__init__(KIND, name, params)
        if params["relation"] not in RELATIONS:
            raise ValueError(f"block {name!r}: relation must be one of {RELATIONS}, got {params['relation']!r}")
        if params["model_stages"] not in MODEL_STAGES:
            raise ValueError(f"block {name!r}: model_stages must be one of {MODEL_STAGES}")
        if params["relation"] == "minus_l" and params.get("reference") is None:
            raise ValueError(f"block {name!r}: relation minus_l needs a reference slot")
        if params["relation"] == "plant":
            require(params, "controller", block=name)

    def bind(self, layout) -> None:
        self.kernel = _noop

    def __call__(self, state, params):
        return state

    @property
    def signal_slots(self) -> list[str]:
        p = self.params
        return [s for s in (p["x"], p.get("reference"), p["response"]) if s is not None]

    def evaluate(self, series: Mapping[str, np.ndarray], t_s: np.ndarray, *, excitation: Any, model: Any,
                 f_sim_hz: float) -> dict[str, Any]:
        """Run :func:`identify` on the recorded ``series`` (``excitation``: the excitation block; ``model``:
        the ``loop_stability_check`` block or None)."""
        p = self.params
        sc = excitation_schedule(excitation.params, excitation.name)
        w = sc["waveform"]
        if w == "stepped_sine":
            require(p, "n_blocks", block=self.name)
        if w in ("chirp", "prbs", "white"):
            require(p, "band_hz", block=self.name)
        if w == "white":
            require(p, "nperseg_s", block=self.name)
        sig = {"x": series[p["x"]], "response": series[p["response"]]}
        if p.get("reference") is not None:
            sig["reference"] = series[p["reference"]]
        return identify(sig, t_s, sc, f_sim_hz, relation=p["relation"], controller=p.get("controller"),
                        model_params=None if model is None else model.params, model_stages=p["model_stages"],
                        coherence_min=float(p["coherence_min"]), n_blocks=int(p.get("n_blocks", 4)),
                        band_hz=p.get("band_hz"), nperseg_s=p.get("nperseg_s"))


@register_block(KIND)
def make_loop_identification(params: Mapping[str, Any], *, name: str | None = None) -> Block:
    """Factory of the ``loop_identification`` block (module docstring)."""
    return LoopIdentification(params, name)
