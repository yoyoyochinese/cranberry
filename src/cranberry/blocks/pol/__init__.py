"""POL: polarization elements (PLAN §4.3 PATH, polarization part; G04, T09). Owner: WP-B7.

Block kinds (factories ``make_<kind>``, registered with ``@register_block``):

* ``segment_jones``: U for ``J.<seg>`` (PM: birefringence + end crosstalk; SMF: D_pol diffusion +
  temperature). M1.
* ``polarizer_projection``: ER projection, power drift δp, output state; IL through ``losses()``. M1.
* ``per_fringe``: Jones → (x, Δχ) slots for ``dc_fringe`` (Ref D adapter). M1.
* ``pdl_element`` (toggle ``pol_pdl``): M4, not implemented yet.
* ``faraday_nonreciprocal`` (optional, toggle ``pol_faraday``): M4, not implemented yet.

The PER accounting blocks ``per_budget`` and ``coil_birefringence`` live in ``blocks.loss``
(PLAN §4.3 LOSS).

Contract: core.block (factory interface, reads/writes, losses(), bands(), toggles, D34 no defaults).
"""

from cranberry.blocks.pol import per_fringe, polarizer_projection, segment_jones  # noqa: F401  (registers kinds)
from cranberry.blocks.pol.per_fringe import mu_to_x_dchi

__all__ = ["mu_to_x_dchi"]
