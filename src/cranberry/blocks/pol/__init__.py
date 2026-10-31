"""POL: polarization elements (PLAN §4.3 PATH, polarization part; G04, T09). Owner: WP-B7.

Block kinds (factories ``make_<kind>``, registered with ``@register_block``):

* ``segment_jones``: U for ``J.<seg>`` (PM: birefringence + end crosstalk; SMF: static U0, D_pol
  diffusion or OU spectrum, temperature). M1 + M4.
* ``polarizer_projection``: ER projection, power drift δp, output state; IL through ``losses()``. M1.
* ``per_fringe``: Jones → (x, Δχ) slots for ``dc_fringe`` (Ref D adapter). M1.
* ``pdl_element`` (toggle ``pol_pdl``): PDL Jones form in a J slot; averaged IL through ``losses()``. M4.
* ``faraday_nonreciprocal`` (optional, toggle ``pol_faraday``): non-reciprocal rotation in a
  ``J_nr_slot`` (not transposed for d = −1). M4.
* ``pmd_mode_overlap`` (offline, toggle ``pol_pmd``): SMF PMD mode-overlap penalty for pulses;
  ``pulse_projection`` for any U(Ω). M4.

The PER accounting blocks ``per_budget`` and ``coil_birefringence`` live in ``blocks.loss``
(PLAN §4.3 LOSS).

Contract: core.block (factory interface, reads/writes, losses(), bands(), toggles, D34 no defaults).
"""

from cranberry.blocks.pol import (  # noqa: F401  (registers kinds)
    faraday_nonreciprocal,
    pdl_element,
    per_fringe,
    pmd_mode_overlap,
    polarizer_projection,
    segment_jones,
)
from cranberry.blocks.pol.per_fringe import mu_to_x_dchi

__all__ = ["mu_to_x_dchi"]
