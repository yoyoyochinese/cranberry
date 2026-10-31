"""Stateful coloured-noise generators shared by the NOISE blocks (G11). Owner: WP-B2.

A coloured process is a sum of independent **linear Gaussian sections**: continuous-time, strictly
proper SISO systems x' = A x + B ξ, y = C x, driven by white noise ξ of one-sided PSD 1, so the
output has the one-sided PSD |H(i2πf)|² (conv §7) with H(s) = C (s − A)⁻¹ B. Building blocks:

* :func:`lowpass`    H = √S0 / (1 + s/ω_c)                   S = S0 / (1 + (f/f_c)²)  (OU process)
* :func:`two_pole`   H = (√h/f_lo) / ((1 + s/ω_lo)(1 + s/ω_c)) S ≈ h/f² for f_lo ≪ f ≪ f_c (T08 eq. 2.7)
* :func:`resonance`  band pass, peak at f0, quality q, area ∫S df = area (technical / mechanical peaks)
* :func:`flicker`    OU sum ≈ h/f on [f_min, f_max]: corners f_k = f_min·r^k, r = 10^(1/m), each with
  S0_k = h·(2 ln r/π)/f_k (Σ_k → h/f in the continuum limit; ripple set by m OUs per decade, G11)
* :func:`cascade`    series connection with a first-order biproper filter (correlation shaping ρ(f))

Discretisation is exact (Van Loan): x_{n+1} = Φ x_n + L w_n, Φ = e^{Ah}, L Lᵀ = ∫_0^h e^{As}B q Bᵀ e^{Aᵀs} ds
(q = 1/2: two-sided intensity of a one-sided-PSD-1 white noise), w_n ~ N(0, 1). The samples are
samples of the continuous process (no discretisation error at any step; only aliasing of the part
above Nyquist — a section whose corner lies above the block Nyquist frequency folds its power into the
loop band, which the ledger then also declares above f_N; run such blocks at f_N ≫ corner). The initial state is drawn from the stationary covariance (continuous Lyapunov
equation), so statistics hold from n = 0.

A :class:`Bank` groups sections that share one ``DrawBuffer`` (one RNG sub-stream, so toggling one
bank off leaves the others' realisations unchanged) and adds each section's output, with gains,
to scalar outputs (slot entries). :func:`bank_step` is written once and run both as pure Python
(reference) and numba-compiled (kernel), so the two agree to rounding (VT-44).
"""

from __future__ import annotations

import math
from collections import namedtuple
from dataclasses import dataclass

import numba
import numpy as np
import scipy.linalg

from cranberry.core.rng import DrawBuffer, substream

__all__ = ["Section", "lowpass", "two_pole", "resonance", "flicker", "cascade", "section_psd", "Bank",
           "BankC", "bank_step", "bank_step_nb", "EMPTY_BANK"]

_Q = 0.5


@dataclass(frozen=True)
class Section:
    A: np.ndarray  # (d, d)
    B: np.ndarray  # (d,)
    C: np.ndarray  # (d,)

    @property
    def dim(self) -> int:
        return len(self.B)


def _w(f_hz: float) -> float:
    return 2.0 * math.pi * f_hz


def lowpass(f_c_hz: float, s0: float) -> Section:
    w = _w(f_c_hz)
    return Section(np.array([[-w]]), np.array([w]), np.array([math.sqrt(s0)]))


def two_pole(h: float, f_lo_hz: float, f_c_hz: float) -> Section:
    a, b = _w(f_lo_hz), _w(f_c_hz)
    return Section(np.array([[-a, 0.0], [b, -b]]), np.array([a, 0.0]), np.array([0.0, math.sqrt(h) / f_lo_hz]))


def resonance(f0_hz: float, q: float, area: float) -> Section:
    """H = K (ω0/q) s / (s² + (ω0/q) s + ω0²), |H(ω0)| = K, ∫S df = K² (f0/q) π/2 = area."""
    w0 = _w(f0_hz)
    k = math.sqrt(2.0 * area * q / (math.pi * f0_hz))
    # balanced states (y, y'/ω0) keep A well scaled for the Lyapunov / expm solves
    return Section(np.array([[0.0, w0], [-w0, -w0 / q]]), np.array([0.0, w0]), np.array([0.0, k / q]))


def flicker(h: float, f_min_hz: float, f_max_hz: float, per_decade: int) -> list[Section]:
    # ponytail: plain OU sum; within ~1 decade of f_min / f_max the PSD falls below h/f (−6 % one decade
    # in, −50 % at the edge); ≥ 2 per decade gives < 1 % ripple inside. Add end-corner weights if needed.
    r = 10.0 ** (1.0 / per_decade)
    n = int(math.floor(math.log(f_max_hz / f_min_hz) / math.log(r) + 1e-9)) + 1
    c = 2.0 * math.log(r) / math.pi
    return [lowpass(f_min_hz * r**k, h * c / (f_min_hz * r**k)) for k in range(n)]


def cascade(sec: Section, f_hz: float, num0: float, num1: float) -> Section:
    """Series ``sec`` → H2(s) = (num0 + num1·s/ω)/(1 + s/ω) (``f_hz = inf``: static gain num0)."""
    if math.isinf(f_hz):
        return Section(sec.A, sec.B, num0 * sec.C)
    w, d = _w(f_hz), sec.dim
    A = np.zeros((d + 1, d + 1))
    A[:d, :d] = sec.A
    A[d, :d] = w * sec.C
    A[d, d] = -w
    B = np.append(sec.B, 0.0)
    # H2 = num1 + (num0 − num1)/(1 + s/ω): state z' = −ω z + ω y_in, out = num1·y_in + (num0 − num1)·z
    C = np.append(num1 * sec.C, num0 - num1)
    return Section(A, B, C)


def section_psd(sec: Section, f_hz) -> np.ndarray:
    """One-sided PSD |C (i2πf − A)⁻¹ B|² of a section (analytic, for checks and ledger rows)."""
    f = np.atleast_1d(np.asarray(f_hz, dtype=float))
    eye = np.eye(sec.dim)
    h = np.array([sec.C @ np.linalg.solve(1j * _w(x) * eye - sec.A, sec.B) for x in f])
    return np.abs(h) ** 2


def _psd_sqrt(m: np.ndarray) -> np.ndarray:
    w, v = np.linalg.eigh(0.5 * (m + m.T))
    return v * np.sqrt(np.clip(w, 0.0, None))


def discretize(sec: Section, h_s: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(Φ, L, P∞^{1/2}) of the exact discretisation (module docstring)."""
    d = sec.dim
    G = _Q * np.outer(sec.B, sec.B)
    M = np.zeros((2 * d, 2 * d))
    M[:d, :d], M[:d, d:], M[d:, d:] = -sec.A, G, sec.A.T
    E = scipy.linalg.expm(M * h_s)
    phi = E[d:, d:].T
    Q = phi @ E[:d, d:]
    P = scipy.linalg.solve_continuous_lyapunov(sec.A, -G)
    return phi, _psd_sqrt(Q), _psd_sqrt(P)


#: compiled-side description of one bank (all ints / flat arrays; see :func:`bank_step`)
BankC = namedtuple("BankC", "on nsec x buf idx D SO PO E0 E1 OUT CO PHI LQ CF")
EMPTY_BANK = BankC(False, 0, -1, -1, -1, *(np.zeros(1, np.int64) for _ in range(7)), *(np.zeros(1) for _ in range(3)))


class Bank:
    """Sections sharing one DrawBuffer; outputs are (slot, flat index, gain) triples."""

    def __init__(self, owner: str, label: str):
        self.owner, self.label = owner, label
        self.items: list[tuple[Section, list[tuple[str, int, float]]]] = []

    def add(self, sec: Section, outputs: list[tuple[str, int, float]]) -> None:
        self.items.append((sec, outputs))

    @property
    def dims(self) -> int:
        return sum(s.dim for s, _ in self.items)

    def init(self, layout, rng, rate_div: int, h_s: float) -> tuple[str, ...]:
        """Discretise, register the state slot (stationary initial draw) and the DrawBuffer."""
        self.state_slot = f"rng.{self.owner}.{self.label}.x"
        self.disc = [discretize(s, h_s) for s, _ in self.items]
        g = substream(rng, f"{self.label}.x0")
        x0 = np.concatenate([p @ g.standard_normal(len(p)) for _, _, p in self.disc])
        layout.add(self.state_slot, (len(x0),), np.float64, "1", init=x0)
        scratch = f"rng.{self.owner}.{self.label}.tmp"
        layout.add(scratch, (2 * max(s.dim for s, _ in self.items),), np.float64, "1")
        self.scratch = scratch
        self.draws = DrawBuffer(self.owner, self.dims, "standard_normal", label=self.label)
        self.draws.init(layout, substream(rng, self.label), rate_div)
        return (self.state_slot, scratch, *self.draws.slot_names)

    def bind(self, layout) -> BankC:
        self.draws.bind(layout)
        D, SO, PO, E0, E1, OUT, CO = ([] for _ in range(7))
        PHI, LQ, CF = [], [], []
        so = 0
        for (sec, outs), (phi, lq, _) in zip(self.items, self.disc):
            d = sec.dim
            D.append(d)
            SO.append(so)
            PO.append(len(PHI))
            PHI += phi.reshape(-1).tolist()
            LQ += lq.reshape(-1).tolist()
            E0.append(len(OUT))
            for slot, j, gain in outs:
                OUT.append(layout.slot(slot).offset + j)
                CO.append(len(CF))
                CF += (gain * sec.C).tolist()
            E1.append(len(OUT))
            so += d
        i64 = lambda v: np.array(v, dtype=np.int64)  # noqa: E731
        f64 = lambda v: np.array(v, dtype=np.float64)  # noqa: E731
        return BankC(True, len(D), layout.slot(self.state_slot).offset, self.draws.buf_off, self.draws.idx_off,
                     i64(D), i64(SO), i64(PO), i64(E0), i64(E1), i64(OUT), i64(CO), f64(PHI), f64(LQ), f64(CF),
                     ) if D else EMPTY_BANK

    def scratch_off(self, layout) -> int:
        return layout.slot(self.scratch).offset

    def refill(self, state, n0: int, n_samples: int) -> None:
        self.draws.refill(state, n0, n_samples)


def bank_step(S, I, b, tmp):  # noqa: N803
    """Add y_n = C x_n of every section to its outputs, then x ← Φ x + L w (``tmp``: scratch offset)."""
    for k in range(b.nsec):
        d = b.D[k]
        xo = b.x + b.SO[k]
        po = b.PO[k]
        for e in range(b.E0[k], b.E1[k]):
            y = 0.0
            for i in range(d):
                y += b.CF[b.CO[e] + i] * S[xo + i]
            S[b.OUT[e]] += y
        for i in range(d):
            S[tmp + d + i] = S[b.buf + I[b.idx]]
            I[b.idx] += 1
        for i in range(d):
            v = 0.0
            for j in range(d):
                v += b.PHI[po + i * d + j] * S[xo + j] + b.LQ[po + i * d + j] * S[tmp + d + j]
            S[tmp + i] = v
        for i in range(d):
            S[xo + i] = S[tmp + i]


bank_step_nb = numba.njit(bank_step)
