"""ACT/TEC: actuators (PLAN §4.3 ACT; thermal_actuator writes temperature, G08).

Owner: WP-B4.

Block kinds (factories ``make_<kind>``, registered with ``@register_block``):

* ``pzt_fiber_stretcher``
* ``hv_amp_914``
* ``pzt_mirror``
* ``thermal_actuator``

Contract: core.block (factory interface, reads/writes, losses(), bands(), toggles, D34 no defaults).
"""
