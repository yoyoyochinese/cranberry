"""CHK: consistency warnings that never stop a run (PLAN §4.3 CHK, G16-j; core.checks.ChkLog).

Owner: WP-B10 (sbs), WP-B11 (rep rate), WP-B8 (demux), WP-B9 (grid edge), Phase 5 M1D-02 (loop stability).

Block kinds (factories ``make_<kind>``, registered with ``@register_block``):

* ``sbs_check``
* ``rep_rate_tail_check``
* ``demux_feasibility``
* ``grid_edge_monitor``
* ``loop_stability_check``

Contract: core.block (factory interface, reads/writes, losses(), bands(), toggles, D34 no defaults).
"""
