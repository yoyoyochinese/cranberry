"""MEAS: measurement and metrics (PLAN §4.3 MEAS). Metric containers: meas.metrics; phase → squeezing map: meas.phase_to_quadrature.

Owner: WP-F4 (metrics, phase_to_quadrature frame), WP-B5 (ESA internals).

Block kinds (factories ``make_<kind>``, registered with ``@register_block``):

* ``phase_to_quadrature``
* ``hf_phase_noise_integral``
* ``vpm_table``
* ``esa_zero_span``
* ``esa_swept_spectrum``
* ``pulse_sample``
* ``metric_i``
* ``metric_ii``
* ``required_theta_contour``

Contract: core.block (factory interface, reads/writes, losses(), bands(), toggles, D34 no defaults).
"""
