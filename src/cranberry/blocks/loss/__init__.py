"""LOSS: loss and mode-matching elements feeding the LossRegistry (PLAN §4.3 LOSS, G06).

Owner: WP-F3 (registry), WP-B5 (detection), WP-B7 (PER, coil), WP-B1 (tap).

Block kinds (factories ``make_<kind>``, registered with ``@register_block``):

* ``loss_budget``
* ``pure_loss``
* ``mode_overlap``
* ``per_budget``
* ``coil_birefringence``
* ``tap``
* ``etalon_noise`` (optional)

Contract: core.block (factory interface, reads/writes, losses(), bands(), toggles, D34 no defaults).
"""
