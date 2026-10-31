"""POL: polarization elements (PLAN §4.3 PATH, polarization part; G04).

Owner: WP-B7.

Block kinds (factories ``make_<kind>``, registered with ``@register_block``):

* ``segment_jones``
* ``polarizer_projection``
* ``pdl_element``
* ``faraday_nonreciprocal`` (optional)

Contract: core.block (factory interface, reads/writes, losses(), bands(), toggles, D34 no defaults).
"""
