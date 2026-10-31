"""Units, dB conversions and the parameter-name unit-suffix convention (G16-i, conv §1, §2, §7).

Rules (G16-i):

* Internal values are SI (conv §1); quantum quadratures are dimensionless with ħ = 1 and vacuum
  variance V_vac = 1/2 (conv §2).
* Parameter names carry a unit suffix (``f_sim_hz``, ``tau_s``, ``L_m``, ``il_db``, ...). The
  recognised suffixes are :data:`UNIT_SUFFIXES`; :func:`unit_of` parses them (longest match wins).
* Values in dB are distinguished by name (``*_db``, ``*_dbc_per_hz``, ``*_db_per_hz``) and converted
  to linear when a config is loaded (WP-F2). Code below the config layer works in linear units.
* State slots carry a free-form unit tag (``SlotLayout.add(..., unit=...)``); the canonical tag
  vocabulary is :data:`UNIT_TAGS`. The wiring check (WP-F1) compares tags.

dB conventions:

* power ratio:      x_dB = 10·log10(x)            (:func:`db_to_lin`, :func:`lin_to_db`)
* amplitude ratio:  x_dB = 20·log10(x)            (:func:`db_to_amp`, :func:`amp_to_db`)
* insertion loss:   IL_dB ≥ 0 for a passive element, η = 10^(−IL_dB/10)   (:func:`il_db_to_eta`)
* noise level:      S_dB = 10·log10(V/V_SN), V_SN = 1/2 (conv §2)          (:func:`noise_db`)
* phase noise:      ℒ(f) [dBc/Hz] with ℒ = S_φ/2 (one-sided S_φ, conv §7)  (:func:`dbc_per_hz_to_s_phi`)
* RIN:              S_RIN [dB/Hz] → 10^(x/10) [1/Hz], one-sided             (:func:`db_per_hz_to_lin`)
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "V_VAC",
    "V_SN",
    "UNIT_SUFFIXES",
    "UNIT_TAGS",
    "unit_of",
    "is_db_name",
    "db_to_lin",
    "lin_to_db",
    "db_to_amp",
    "amp_to_db",
    "il_db_to_eta",
    "eta_to_il_db",
    "noise_db",
    "noise_db_to_variance",
    "dbc_per_hz_to_s_phi",
    "s_phi_to_dbc_per_hz",
    "db_per_hz_to_lin",
    "lin_to_db_per_hz",
]

#: Vacuum quadrature variance (conv §2, ħ = 1, â = (q̂ + i p̂)/√2).
V_VAC: float = 0.5
#: Shot-noise reference variance V_SN = V_vac (conv §2).
V_SN: float = V_VAC

#: Parameter-name suffix → unit (G16-i). Longest suffix wins in :func:`unit_of`.
#: Extending this table is an additive, non-breaking change (propose it in the WP report).
UNIT_SUFFIXES: dict[str, str] = {
    "_hz": "Hz",
    "_s": "s",
    "_m": "m",
    "_w": "W",
    "_v": "V",
    "_a": "A",
    "_k": "K",
    "_rad": "rad",
    "_db": "dB",
    "_per_s": "1/s",
    "_rad_per_s": "rad/s",
    "_rad_per_v": "rad/V",
    "_v_per_a": "V/A",
    "_a_per_w": "A/W",
    "_m_per_v": "m/V",
    "_rad_per_k": "rad/K",
    "_rad_per_m_k": "rad/(m·K)",
    "_k_per_s": "K/s",
    "_s2_per_m": "s²/m",
    "_s3_per_m": "s³/m",
    "_per_w_m": "1/(W·m)",
    "_per_hz": "1/Hz",
    "_rad2_per_hz": "rad²/Hz",
    "_hz2_per_hz": "Hz²/Hz",
    "_db_per_hz": "dB/Hz",
    "_dbc_per_hz": "dBc/Hz",
    "_ohm": "Ω",
}

#: Canonical unit tags for State slots ("1" = dimensionless, e.g. quadrature variance, η, counts).
UNIT_TAGS: frozenset[str] = frozenset(
    {"1", "s", "Hz", "rad", "rad/s", "V", "A", "W", "K", "m", "rad²", "V²", "1/Hz", "rad²/Hz"}
)

_DB_SUFFIXES = ("_db", "_db_per_hz", "_dbc_per_hz")


def unit_of(name: str) -> str | None:
    """Unit implied by a parameter name's suffix (G16-i), or None if it has no recognised suffix.

    >>> unit_of("f_sim_hz"), unit_of("gain_v_per_a"), unit_of("dL")
    ('Hz', 'V/A', None)
    """
    lname = name.lower()
    best: str | None = None
    for suffix in UNIT_SUFFIXES:
        if lname.endswith(suffix) and len(lname) > len(suffix):
            if best is None or len(suffix) > len(best):
                best = suffix
    return None if best is None else UNIT_SUFFIXES[best]


def is_db_name(name: str) -> bool:
    """True if the parameter name marks a logarithmic (dB) value that must be linearised on load."""
    lname = name.lower()
    return any(lname.endswith(s) and len(lname) > len(s) for s in _DB_SUFFIXES)


def db_to_lin(x_db):
    """Power ratio from dB: 10^(x/10)."""
    return np.power(10.0, np.asarray(x_db, dtype=float) / 10.0)


def lin_to_db(x):
    """dB from a power ratio: 10·log10(x). Raises for negative input; 0 maps to −inf."""
    x = np.asarray(x, dtype=float)
    if np.any(x < 0):
        raise ValueError("lin_to_db: power ratio must be >= 0")
    with np.errstate(divide="ignore"):
        return 10.0 * np.log10(x)


def db_to_amp(x_db):
    """Amplitude ratio from dB: 10^(x/20)."""
    return np.power(10.0, np.asarray(x_db, dtype=float) / 20.0)


def amp_to_db(x):
    """dB from an amplitude ratio: 20·log10(|x|)."""
    with np.errstate(divide="ignore"):
        return 20.0 * np.log10(np.abs(np.asarray(x, dtype=float)))


def il_db_to_eta(il_db):
    """Power transmission η = 10^(−IL/10) from insertion loss IL [dB] (IL ≥ 0 for passive elements)."""
    return np.power(10.0, -np.asarray(il_db, dtype=float) / 10.0)


def eta_to_il_db(eta):
    """Insertion loss IL = −10·log10(η) [dB] from power transmission η ∈ (0, 1]."""
    eta = np.asarray(eta, dtype=float)
    if np.any(eta <= 0):
        raise ValueError("eta_to_il_db: transmission must be > 0")
    return -10.0 * np.log10(eta)


def noise_db(variance, v_sn: float = V_SN):
    """Noise level S_dB = 10·log10(V/V_SN) relative to shot noise (conv §2, V_SN = 1/2)."""
    return lin_to_db(np.asarray(variance, dtype=float) / v_sn)


def noise_db_to_variance(s_db, v_sn: float = V_SN):
    """Quadrature variance V = V_SN·10^(S_dB/10) (inverse of :func:`noise_db`)."""
    return v_sn * db_to_lin(s_db)


def dbc_per_hz_to_s_phi(l_dbc_per_hz):
    """One-sided phase-noise PSD S_φ [rad²/Hz] from single-sideband ℒ(f) [dBc/Hz]: S_φ = 2·10^(ℒ/10)."""
    return 2.0 * db_to_lin(l_dbc_per_hz)


def s_phi_to_dbc_per_hz(s_phi):
    """ℒ(f) [dBc/Hz] = 10·log10(S_φ/2) from one-sided S_φ [rad²/Hz] (conv §7)."""
    return lin_to_db(np.asarray(s_phi, dtype=float) / 2.0)


def db_per_hz_to_lin(x_db_per_hz):
    """Linear one-sided PSD [1/Hz] from a dB/Hz value (e.g. RIN −140 dB/Hz → 1e−14 /Hz)."""
    return db_to_lin(x_db_per_hz)


def lin_to_db_per_hz(x):
    """dB/Hz from a linear one-sided PSD [1/Hz]."""
    return lin_to_db(x)
