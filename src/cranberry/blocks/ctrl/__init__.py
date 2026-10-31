"""CTRL: controllers and digital chain (PLAN §4.3 CTRL; composition is core.pipeline).

Owner: WP-B4.

Block kinds (factories ``make_<kind>``, registered with ``@register_block``):

* ``pid``
* ``adc``
* ``dac``
* ``iir``
* ``hpf``
* ``lpf``
* ``adder``
* ``ramp``
* ``lock_supervisor``

Contract: core.block (factory interface, reads/writes, losses(), bands(), toggles, D34 no defaults).
"""
