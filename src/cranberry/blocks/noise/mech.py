"""Segment phase-noise sources: ``segment_mech_noise``, ``segment_thermodynamic`` (opt), ``gawbs_phase_jitter``
(opt) (T08 §2.5–§2.8, T07). Owner: WP-B2.

``segment_mech_noise`` (T08 eq. 2.7, 2.12): acoustic / vibration strain on ``drv.<seg>.strain`` [1].
The phase PSD of segment i (reference axis, coefficient k_ε·L_i) is

    S_ψ,i(f) = (L_i/L_ref)^γ · [ h_ref f⁻² / ((1 + (f_lo/f)²)(1 + (f/f_c)²)) + Σ_k P_k(f) ]

(f⁻² region with the f_c roll-off; flattened below f_lo where the thermal model takes over; P_k =
band-pass peaks {f_hz, q, area_rad2}). It is written as strain: S_ε,i = S_ψ,i/(k_ε L_i)²
(``k_eps_rad_per_m``; the topology coefficient should be k_ε,a·L). γ ∈ {1, 2} (length scaling) is a
parameter, independent of the spatial placement.

* Inter-segment correlation (toggle ``correlation``): every coloured section is split into a part
  common to all segments of the block, shaped by √ρ(f), and an independent part shaped by
  √(1 − ρ(f)), with ρ(f) = ρ_0/(1 + (f/f_ρ)²) (rational form of T08's ρ_0/√(1 + (f/f_ρ)²) so it is a
  first-order section; ``f_rho_hz: .inf`` gives a constant ρ_0). Cross-PSD S_ij = ρ(f)√(S_i S_j),
  coherence ρ², difference of equal segments 2S(1 − ρ) (T08 eq. 2.13). Off: independent segments.
* Spatial structure (conv §9, T08 §2.7): a ``per_sub`` placement gets a part common to the N_sub
  subsegments (fraction κ = ``spatial_coherence[seg]`` of the segment-phase PSD; it carries the
  inter-segment correlation) and independent subsegment parts (PSD (1 − κ)·N·S each, so their average
  contributes (1 − κ)·S). κ = 0: distributed, 1: lumped. Scalar placements are lumped.
* Toggles {mech, peaks, correlation}; mech and peaks draw from separate RNG sub-streams.

``segment_thermodynamic`` (T08 §2.6, Dong 2016 Table I): Wanser thermal conduction (eq. 2.8) +
Duan thermomechanical (eq. 2.9) phase noise, ∝ L, independent per subsegment. Time domain:
conduction = OU at its DC level with corner k_min² D_th/(2π); thermomechanical = OU sum ≈ A/f from
``tm_f_min_hz`` to Nyquist. Written in driver units through ``k_rad_per_m_unit`` (phase per metre per
unit of ``driver``). With ``accumulate: true`` it adds to a driver slot another block writes (e.g. the
mechanical strain); otherwise it registers the slot. Toggles {conduction, thermomechanical}.
:func:`s_tc` / :func:`s_tm` are the offline PSDs (T08 V9).

``gawbs_phase_jitter`` (T07; GAWBS lines are at 10s–100s of MHz, above any loop band): nothing in the
time loop; declares ``variance_outside = sigma2_rad2`` (σ_G², from ``gawbs_table`` of B10) as a
``phase`` band on ``targets``. Toggle {gawbs}.

Ledger (G05): each segment declares its strain (or driver) PSD above the block Nyquist frequency.
"""

from __future__ import annotations

import math
from collections import namedtuple
from collections.abc import Mapping
from typing import Any

import numba
import numpy as np

from cranberry.blocks.noise._lgs import EMPTY_BANK, Bank, bank_step, bank_step_nb, cascade, flicker, lowpass, resonance, two_pole
from cranberry.blocks.noise.laser import resonance_psd
from cranberry.core.block import BaseBlock, Block, require
from cranberry.core.rates import block_dt_s
from cranberry.core.registry import register_block
from cranberry.offline.band_ledger import BandDecl
from cranberry.path.segments import select_ids
from cranberry.path.topology import as_topology

__all__ = ["SegmentMechNoise", "SegmentThermodynamic", "GawbsPhaseJitter", "make_segment_mech_noise",
           "make_segment_thermodynamic", "make_gawbs_phase_jitter", "s_tc", "s_tm", "two_pole_psd"]

_KB = 1.380649e-23


def two_pole_psd(f, h: float, f_lo_hz: float, f_c_hz: float) -> np.ndarray:
    f = np.asarray(f, dtype=float)
    return h / f_lo_hz**2 / ((1.0 + (f / f_lo_hz) ** 2) * (1.0 + (f / f_c_hz) ** 2))


def s_tc(f, L_m, *, T_k, wavelength_m, dn_dT_per_k, n, alpha_l_per_k, kappa_th_w_per_m_k, d_th_m2_per_s,
         k_max_per_m, k_min_per_m, **_):
    """Wanser thermal-conduction phase PSD [rad²/Hz] (T08 eq. 2.8)."""
    w = 2.0 * np.pi * np.asarray(f, dtype=float) / d_th_m2_per_s
    pre = 2.0 * np.pi * _KB * T_k**2 * L_m / (kappa_th_w_per_m_k * wavelength_m**2) * (dn_dT_per_k + n * alpha_l_per_k) ** 2
    return pre * np.log((k_max_per_m**4 + w * w) / (k_min_per_m**4 + w * w))


def s_tm(f, L_m, *, T_k, wavelength_m, n, e0_pa, area_m2, phi0, **_):
    """Duan thermomechanical phase PSD [rad²/Hz] (T08 eq. 2.9)."""
    return (2.0 * np.pi * n / wavelength_m) ** 2 * 2.0 * _KB * T_k * L_m * phi0 / (3.0 * np.pi * e0_pa * area_m2 * np.asarray(f, dtype=float))


BC = namedtuple("BC", "zero nzero ZO b1 t1 b2 t2")


def _build(jit: bool):
    bs = bank_step_nb if jit else bank_step

    def step(S, I, c):  # noqa: N803
        if c.zero:
            for j in range(c.nzero):
                S[c.ZO[j]] = 0.0
        if c.b1.on:
            bs(S, I, c.b1, c.t1)
        if c.b2.on:
            bs(S, I, c.b2, c.t2)

    return numba.njit(step) if jit else step


_STEP_PY, _STEP_NB = _build(False), _build(True)


class _BankBlock(BaseBlock):
    """Outputs = per-segment driver slots fed by up to two banks (``BANKS`` = toggle names)."""

    KIND = ""
    BANKS: tuple[str, str] = ("", "")

    def __init__(self, params, name, driver, toggles, accumulate=False):
        name = name or self.KIND
        require(params, "topology", "segments", "toggles", block=name)
        topo = as_topology(params["topology"])
        self.segs = []
        for s in select_ids(params["segments"], topo.segments, "segments", name):
            seg = topo.segments[s]
            if driver not in seg.drivers:
                raise ValueError(f"block {name!r}: segment {s!r} declares no {driver} driver")
            per_sub = seg.drivers[driver][0] == "per_sub"
            self.segs.append((s, seg.length_m, seg.n_sub if per_sub else 0))
        self.driver, self.accumulate = driver, accumulate
        self.slots = [f"drv.{s}.{driver}" for s, _, _ in self.segs]
        super().__init__(self.KIND, name, params, reads=tuple(self.slots) if accumulate else (),
                         writes=tuple(self.slots), toggle_names=toggles)
        self.slot_units = dict.fromkeys(self.slots, self.UNIT)
        self.banks: dict[str, Bank] = {}
        self._h: float | None = None

    def _fill(self, h: float, f_n: float) -> None:  # pragma: no cover - overridden
        raise NotImplementedError

    def init(self, layout, rng):
        super().init(layout, rng)
        self._h = h = block_dt_s(layout.f_sim_hz, self.rate_div)
        if not self.accumulate:
            for slot, (_, _, n) in zip(self.slots, self.segs):
                layout.add(slot, (n,) if n else (), np.float64, self.UNIT)
        self._fill(h, 0.5 / h)
        extra: list[str] = []
        for b in self.banks.values():
            extra += b.init(layout, rng, self.rate_div, h)
        self.writes = (*self.writes, *extra)

    def bind(self, layout):
        zo = [layout.slot(s).offset + j for s, (_, _, n) in zip(self.slots, self.segs) for j in range(max(n, 1))]
        bk = [(self.banks[k].bind(layout), self.banks[k].scratch_off(layout)) if k in self.banks else (EMPTY_BANK, -1)
              for k in self.BANKS]
        self._c = BC(not self.accumulate, len(zo), np.array(zo, dtype=np.int64), bk[0][0], bk[0][1], bk[1][0], bk[1][1])
        c, nb = self._c, _STEP_NB

        @numba.njit
        def kernel(S, Z, I, P):  # noqa: N803
            nb(S, I, c)

        self.kernel = kernel

    def prepare(self, state, params, n0, n_samples):
        for b in self.banks.values():
            b.refill(state, n0, n_samples)

    def __call__(self, state, params):
        _STEP_PY(state.S, state.I, self._c)
        return state

    def _outputs(self, i: int, gain: float, kappa: float) -> list[tuple[str, int, float]]:
        """Entries of segment i's slot fed by its (spatially) common part."""
        _, _, n = self.segs[i]
        return [(self.slots[i], j, gain * math.sqrt(kappa)) for j in range(max(n, 1))]

    def _bank(self, label: str) -> Bank:
        return self.banks.setdefault(label, Bank(self.name, label))

    def _psd(self, i: int):  # pragma: no cover - overridden
        raise NotImplementedError

    def bands(self) -> list[BandDecl]:
        if self._h is None:
            raise RuntimeError(f"block {self.name!r}: bands() needs init (block rate unknown)")
        if not any(self.toggles.values()):
            psds = [None] * len(self.segs)
        else:
            psds = [self._psd(i) for i in range(len(self.segs))]
        return [BandDecl(self.name, self.driver, (s,), 0.0, 0.5 / self._h, self.UNIT, psd_outside=p,
                         variance_outside=None if p is not None else 0.0)
                for (s, _, _), p in zip(self.segs, psds)]


class SegmentMechNoise(_BankBlock):
    KIND, UNIT, BANKS = "segment_mech_noise", "1", ("mech", "peaks")

    def __init__(self, params, name):
        name = name or self.KIND
        require(params, "L_ref_m", "gamma", "k_eps_rad_per_m", "toggles", block=name)
        tg = params["toggles"]
        if tg.get("mech"):
            require(params, "h_ref_rad2_hz", "f_c_hz", "f_lo_hz", block=name)
        if tg.get("peaks"):
            require(params, "peaks", block=name)
        if tg.get("correlation"):
            require(params, "rho_0", "f_rho_hz", block=name)
        super().__init__(params, name, "strain", ("mech", "peaks", "correlation"))
        for s, _, n in self.segs:
            if n:
                require(params, "spatial_coherence", block=name)
                if not 0.0 <= float(params["spatial_coherence"][s]) <= 1.0:
                    raise ValueError(f"block {name!r}: spatial_coherence[{s!r}] must be in [0, 1]")

    def _amp(self, i: int) -> float:
        p = self.params
        _, L, _ = self.segs[i]
        return math.sqrt((L / float(p["L_ref_m"])) ** float(p["gamma"])) / (float(p["k_eps_rad_per_m"]) * L)

    def _kappa(self, i: int) -> float:
        s, _, n = self.segs[i]
        return float(self.params["spatial_coherence"][s]) if n else 1.0

    def _units(self) -> dict[str, list]:
        p, tg = self.params, self.toggles
        units: dict[str, list] = {}
        if tg["mech"]:
            units["mech"] = [two_pole(float(p["h_ref_rad2_hz"]), float(p["f_lo_hz"]), float(p["f_c_hz"]))]
        if tg["peaks"] and p["peaks"]:
            units["peaks"] = [resonance(float(k["f_hz"]), float(k["q"]), float(k["area_rad2"])) for k in p["peaks"]]
        return units

    def _fill(self, h, f_n):
        p, tg = self.params, self.toggles
        rho0 = float(p["rho_0"]) if tg["correlation"] else 0.0
        f_rho = float(p["f_rho_hz"]) if tg["correlation"] else math.inf
        idx = range(len(self.segs))
        for label, secs in self._units().items():
            b = self._bank(label)
            for sec in secs:
                if tg["correlation"]:
                    b.add(cascade(sec, f_rho, math.sqrt(rho0), 0.0),
                          [o for i in idx for o in self._outputs(i, self._amp(i), self._kappa(i))])
                    ind = cascade(sec, f_rho, math.sqrt(1.0 - rho0), 1.0)
                else:
                    ind = sec
                for i in idx:
                    b.add(ind, self._outputs(i, self._amp(i), self._kappa(i)))
                    _, _, n = self.segs[i]
                    k = self._kappa(i)
                    if n and k < 1.0:
                        g = self._amp(i) * math.sqrt((1.0 - k) * n)
                        for m in range(n):
                            b.add(sec, [(self.slots[i], m, g)])

    def _psd(self, i):
        p, tg, a2 = self.params, self.toggles, self._amp(i) ** 2

        def psd(f):
            s = np.zeros_like(np.asarray(f, dtype=float))
            if tg["mech"]:
                s = s + two_pole_psd(f, float(p["h_ref_rad2_hz"]), float(p["f_lo_hz"]), float(p["f_c_hz"]))
            if tg["peaks"]:
                for k in p["peaks"]:
                    s = s + resonance_psd(f, float(k["f_hz"]), float(k["q"]), float(k["area_rad2"]))
            return a2 * s

        return psd


_MATERIAL = ("T_k", "wavelength_m", "dn_dT_per_k", "n", "alpha_l_per_k", "kappa_th_w_per_m_k", "d_th_m2_per_s",
             "k_max_per_m", "k_min_per_m", "e0_pa", "area_m2", "phi0")


class SegmentThermodynamic(_BankBlock):
    KIND, BANKS = "segment_thermodynamic", ("conduction", "thermomechanical")

    def __init__(self, params, name):
        name = name or self.KIND
        require(params, "driver", "k_rad_per_m_unit", "accumulate", "material", "tm_f_min_hz", "tm_per_decade",
                block=name)
        require(params["material"], *_MATERIAL, block=name)
        if params["driver"] not in ("dT", "strain"):
            raise ValueError(f"block {name!r}: driver must be dT or strain")
        self.UNIT = "K" if params["driver"] == "dT" else "1"
        super().__init__(params, name, params["driver"], ("conduction", "thermomechanical"), bool(params["accumulate"]))

    def _gain2(self, i: int) -> float:
        """(1/(k·L))²: phase PSD → driver PSD of the segment average."""
        _, L, _ = self.segs[i]
        return 1.0 / (float(self.params["k_rad_per_m_unit"]) * L) ** 2

    def _fill(self, h, f_n):
        m, tg = dict(self.params["material"]), self.toggles
        for i, (_, L, n) in enumerate(self.segs):
            g = math.sqrt(self._gain2(i) * max(n, 1))  # independent per subsegment (∝ L)
            if tg["conduction"]:
                f_k = m["k_min_per_m"] ** 2 * m["d_th_m2_per_s"] / (2.0 * math.pi)
                # ponytail: conduction PSD above its corner approximated by f⁻² (true: logarithmic); negligible noise
                sec = lowpass(f_k, float(s_tc(1e-9, L, **m)))
                for j in range(max(n, 1)):
                    self._bank("conduction").add(sec, [(self.slots[i], j, g)])
            if tg["thermomechanical"]:
                a_tm = float(s_tm(1.0, L, **m))
                for sec in flicker(a_tm, float(self.params["tm_f_min_hz"]), f_n, int(self.params["tm_per_decade"])):
                    for j in range(max(n, 1)):
                        self._bank("thermomechanical").add(sec, [(self.slots[i], j, g)])

    def _psd(self, i):
        m, tg, L, g2 = dict(self.params["material"]), self.toggles, self.segs[i][1], self._gain2(i)
        return lambda f: g2 * ((s_tc(f, L, **m) if tg["conduction"] else 0.0) + (s_tm(f, L, **m) if tg["thermomechanical"] else 0.0))


class GawbsPhaseJitter(BaseBlock):
    """``gawbs_phase_jitter`` (module docstring): ledger only."""

    toggle_names = ("gawbs",)

    def __init__(self, params, name):
        name = name or "gawbs_phase_jitter"
        require(params, "targets", "sigma2_rad2", "toggles", block=name)
        super().__init__("gawbs_phase_jitter", name, params)
        self._h: float | None = None

    def init(self, layout, rng):
        super().init(layout, rng)
        self._h = block_dt_s(layout.f_sim_hz, self.rate_div)

    def bind(self, layout):
        @numba.njit
        def kernel(S, Z, I, P):  # noqa: N803
            pass

        self.kernel = kernel

    def __call__(self, state, params):
        return state

    def bands(self) -> list[BandDecl]:
        if self._h is None:
            raise RuntimeError(f"block {self.name!r}: bands() needs init (block rate unknown)")
        v = float(self.params["sigma2_rad2"]) if self.toggles["gawbs"] else 0.0
        return [BandDecl(self.name, "phase", tuple(self.params["targets"]), 0.0, 0.5 / self._h, "rad",
                         variance_outside=v)]


def _register(kind, cls):
    @register_block(kind)
    def factory(params: Mapping[str, Any], *, name: str | None = None) -> Block:
        return cls(params, name)

    factory.__name__ = factory.__qualname__ = f"make_{kind}"
    return factory


make_segment_mech_noise = _register("segment_mech_noise", SegmentMechNoise)
make_segment_thermodynamic = _register("segment_thermodynamic", SegmentThermodynamic)
make_gawbs_phase_jitter = _register("gawbs_phase_jitter", GawbsPhaseJitter)
