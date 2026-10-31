"""DET: detectors (PLAN §4.3 DET).

Owner: WP-B5.

Block kinds (factories ``make_<kind>``, registered with ``@register_block``):

* ``detector_efficiency``
* ``bhd_frontend``
* ``pulsed_bhd_kernel``

Contract: core.block (factory interface, reads/writes, losses(), bands(), toggles, D34 no defaults).
"""
