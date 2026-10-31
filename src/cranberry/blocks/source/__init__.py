"""SRC: sources and pump (PLAN §4.3 SRC).

Owner: WP-B1 (single-mode, waveguide OPA), WP-B6 (CW OPO), WP-B11 (pulsed), WP-B2 (laser).

Block kinds (factories ``make_<kind>``, registered with ``@register_block``):

* ``laser_source``
* ``single_mode_reference``
* ``opa_waveguide``
* ``opa_supermodes``
* ``opo_cw_spectrum``
* ``pulsed_opo_source``
* ``opo_detuning_noise``
* ``opo_length_lock``
* ``opo_angle_map``
* ``pump_power_stabilization``
* ``crystal_thermal``

Contract: core.block (factory interface, reads/writes, losses(), bands(), toggles, D34 no defaults).
"""
