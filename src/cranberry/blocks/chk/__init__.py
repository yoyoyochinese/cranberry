"""CHK: consistency warnings that never stop a run (PLAN §4.3 CHK, G16-j; core.checks.ChkLog).

Owner: WP-B10 (sbs), WP-B11 (rep rate), WP-B8 (demux), WP-B9 (grid edge), Phase 5 M1D-02 (loop stability; stacked loops M2C-02).

Block kinds (factories ``make_<kind>``, registered with ``@register_block``):

* ``sbs_check``
* ``rep_rate_tail_check``
* ``demux_feasibility``
* ``grid_edge_monitor``
* ``loop_stability_check``
* ``gated_loop_stability_check`` (WP6-G): Floquet stability of a periodically gated loop, residual in the measurement
  window, Dick (frequency-conversion) term, backscatter round-trip limit

Contract: core.block (factory interface, reads/writes, losses(), bands(), toggles, D34 no defaults).
"""
