"""FSR: OPO free-spectral-range comb filters (PLAN §4.3 FSR).

Owner: WP-B11.

Block kinds (factories ``make_<kind>``, registered with ``@register_block``):

* ``fsr_comb_filter``

Contract: core.block (factory interface, reads/writes, losses(), bands(), toggles, D34 no defaults).
"""
