"""CTRL: controllers and digital chain (PLAN §4.3 CTRL; composition is core.pipeline).

Owner: WP-B4.

Block kinds (factories ``make_<kind>``, registered with ``@register_block``):

* ``pid`` (``pid.py``): canonical gains, discretisation ``dsc1`` / ``tustin`` / ``forward_euler``,
  anti-windup ``none`` / ``clamp`` / ``back_calculation``
* ``adc``, ``dac`` (``digital.py``): quantisation, range, averaging, noise, output divider
* ``iir``, ``lpf``, ``hpf``, ``band_split`` (``filters.py``): biquad cascades
* ``adder``, ``ramp``, ``sine_injection`` (``signals.py``)
* ``lock_supervisor`` (``supervisor.py``): lock-loss detection, fringe slips, re-acquisition

Shared plumbing (one step function for the reference and the kernel, driver slots, losses):
``_kit.py``. Contract: core.block (factory interface, reads/writes, losses(), bands(), toggles,
D34 no defaults).
"""
