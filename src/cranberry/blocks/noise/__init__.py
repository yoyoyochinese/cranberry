"""NOISE: environmental and electronic noise sources (PLAN §4.3 NOISE; stateful filters, G11).

Owner: WP-B2.

Block kinds (factories ``make_<kind>``, registered with ``@register_block``):

* ``ambient_temperature``
* ``segment_thermal``
* ``segment_mech_noise``
* ``segment_thermodynamic`` (optional)
* ``gawbs_phase_jitter`` (optional)
* ``pd_offset_drift``

Contract: core.block (factory interface, reads/writes, losses(), bands(), toggles, D34 no defaults).
"""
