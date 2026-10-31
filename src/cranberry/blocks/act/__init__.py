"""ACT/TEC: actuators (PLAN §4.3 ACT; thermal_actuator writes temperature, G08).

Owner: WP-B4.

Block kinds (factories ``make_<kind>``, registered with ``@register_block``):

* ``pzt_fiber_stretcher``, ``pzt_mirror``, ``hv_amp_914`` (``pzt.py``)
* ``thermal_actuator`` (``thermal.py``)

Actuators write ``act.<name>.*`` and the segment driver ``drv.<segment>.<driver>`` (F6 names; see
``blocks/ctrl/_kit.py``). Contract: core.block (factory interface, reads/writes, losses(), bands(),
toggles, D34 no defaults).
"""
