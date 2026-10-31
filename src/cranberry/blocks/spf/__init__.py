"""SPF: spectral filters and dispersion compensation (PLAN §4.3 SPF; LO filter cavity placed here, see note).

Owner: WP-B11 (filter_cavity: WP-B2).

Block kinds (factories ``make_<kind>``, registered with ``@register_block``):

* ``spectral_filter``
* ``dispersion_comp``
* ``filter_cavity``

Contract: core.block (factory interface, reads/writes, losses(), bands(), toggles, D34 no defaults).
"""
