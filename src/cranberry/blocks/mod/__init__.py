"""MOD: modulators (PLAN §4.3 MOD). Squeezed-path modulators register their insertion loss (G06).

Owner: WP-B4 (``probe_pulser``: WP-B8).

Block kinds (factories ``make_<kind>``, registered with ``@register_block``):

* ``eom_phase``, ``dither_modulator``, ``aom_shift`` (``modulators.py``)
* ``probe_pulser`` (WP-B8, Liu-type probe; not implemented here)

Contract: core.block (factory interface, reads/writes, losses(), bands(), toggles, D34 no defaults).
"""
