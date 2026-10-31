"""ERR: error-signal generators #1-#6 (PLAN §4.3 ERR; baseband API at loop rate, G16-b).

Owner: WP-B3.

Block kinds (factories ``make_<kind>``, registered with ``@register_block``):

* ``parametric_beat_detector``
* ``bhd_probe_beat``
* ``quantum_noise_power``
* ``dc_fringe``
* ``liu_probe_interferometer``
* ``demodulator``

Contract: core.block (factory interface, reads/writes, losses(), bands(), toggles, D34 no defaults).
"""
