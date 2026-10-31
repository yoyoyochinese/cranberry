"""FIB: fiber offline computations, GNLSE (term toggles), quantum layer, scattering tables, mode overlap (PLAN §4.3 FIB).

Owner: WP-B9 (classical), WP-B10 (quantum, scattering), WP-B11 (mode overlap).

Block kinds (factories ``make_<kind>``, registered with ``@register_block``):

* ``fiber_gnlse``
* ``fiber_linear_table``
* ``lo_nonlinear_table``
* ``xpm_counterprop_coeff``
* ``spm_rin_to_phase``
* ``fiber_quantum_linearized``
* ``fiber_quantum_tw`` (optional)
* ``raman_noise_table``
* ``gawbs_table``
* ``brillouin_table``
* ``stray_light_table``
* ``path_mode_transfer``
* ``mode_overlap_table``
* ``timing_state``
* ``mode_overlap_lookup``

Contract: core.block (factory interface, reads/writes, losses(), bands(), toggles, D34 no defaults).
"""
