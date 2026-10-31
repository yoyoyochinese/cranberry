"""DMX: time demultiplexing to ports (PLAN §4.3 DMX).

Owner: WP-B8.

Block kinds (factories ``make_<kind>``, registered with ``@register_block``):

* ``demux_design`` (offline table, loss registration, CHK ``demux_feasibility`` = ``check()``)
* ``demux_switch`` (time loop: switch drive phase noise, bias drift)
* ``delay_line_k`` (offline length + time-loop timing error δt_k)

all in ``demux.py``; Liu-type lock helpers (actuator position → loop delay, probe through the demux)
in ``liu.py``. Per-port stray-light gating: ``blocks.mod.probe_pulser``.

Contract: core.block (factory interface, reads/writes, losses(), bands(), toggles, D34 no defaults).
"""
