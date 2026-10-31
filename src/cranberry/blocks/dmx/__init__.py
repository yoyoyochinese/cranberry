"""DMX: time demultiplexing to ports (PLAN §4.3 DMX).

Owner: WP-B8.

Block kinds (factories ``make_<kind>``, registered with ``@register_block``):

* ``demux_design``
* ``demux_switch``
* ``delay_line_k``

Contract: core.block (factory interface, reads/writes, losses(), bands(), toggles, D34 no defaults).
"""
