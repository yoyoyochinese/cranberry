"""DET: detectors (PLAN §4.3 DET; T13).

Owner: WP-B5.

Block kinds (factories ``make_<kind>``, registered with ``@register_block``):

* ``detector_efficiency`` (det.bhd, offline): QE, circuit noise (clearance), imbalance leak → LossRegistry
* ``bhd_frontend`` (det.bhd, measurement): RF-output PSDs, DC / MONITOR voltages, clipping, coupler sign
* ``pulsed_bhd_kernel`` (det.pulsed, offline): inter-pulse crosstalk c_m, per-pulse clearance

Device presets: ``devices/pdb450c.yaml`` (Ref D, gain 10³), ``devices/bhd_generic_highqe.yaml`` (F, Open).
The Ref D DC fringe itself is ``err.dc_fringe`` (G16-d); its PDB450C parameters use the same key names.

Contract: core.block (factory interface, reads/writes, losses(), bands(), toggles, D34 no defaults).
"""
