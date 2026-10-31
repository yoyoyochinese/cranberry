"""NOISE: environmental and electronic noise sources, and the laser (PLAN §4.3 NOISE / SRC; stateful
filters, G11). Owner: WP-B2.

Block kinds (factories ``make_<kind>``, registered with ``@register_block``):

* ``laser_source``           (``laser``: white / flicker FM, peaks, drift, dither, SBS, RIN; ledger)
* ``segment_thermal``        (``thermal``: zone ambient temperature = the former ``ambient_temperature``,
  first-order lag per segment)
* ``crystal_thermal``        (``thermal``: crystal ΔT and Δk for the OPA mismatch rotation)
* ``pd_offset_drift``        (``thermal``: PD / TIA DC offset drift)
* ``segment_mech_noise``     (``mech``: h/f², peaks, inter-segment ρ(f), spatial structure)
* ``segment_thermodynamic``  (``mech``, optional: Wanser + thermomechanical)
* ``gawbs_phase_jitter``     (``mech``, optional: ledger only)

Coloured noise: ``_lgs`` (exactly discretised linear Gaussian sections, one RNG sub-stream per bank).
Contract: core.block (factory interface, reads/writes, losses(), bands(), toggles, D34 no defaults).
"""
