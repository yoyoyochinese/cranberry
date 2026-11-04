"""CTRL: controllers and digital chain (PLAN §4.3 CTRL; composition is core.pipeline).

Owner: WP-B4.

Block kinds (factories ``make_<kind>``, registered with ``@register_block``):

* ``pid`` (``pid.py``): canonical gains, discretisation ``dsc1`` / ``tustin`` / ``forward_euler``,
  anti-windup ``none`` / ``clamp`` / ``back_calculation``
* ``adc``, ``dac`` (``digital.py``): quantisation, range, averaging, noise, output divider
* ``iir``, ``lpf``, ``hpf``, ``band_split`` (``filters.py``): biquad cascades
* ``adder``, ``ramp``, ``sine_injection`` (``signals.py``)
* ``excitation`` (``excitation.py``): stepped sine / multisine / chirp / PRBS / white, in-lock identification (D85-7)
* ``lock_supervisor`` (``supervisor.py``): lock-loss detection, fringe slips, re-acquisition
* ``gated_hold`` (``gated_hold.py``, WP6-G): error passed only inside a gate window; zero / hold / noise_only outside,
  ``settle_s`` blanking after each gate edge

Shared plumbing (one step function for the reference and the kernel, driver slots, losses):
``_kit.py``. Contract: core.block (factory interface, reads/writes, losses(), bands(), toggles,
D34 no defaults).
"""
