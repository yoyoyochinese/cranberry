"""Block families (PLAN §4.1, §4.3), one subpackage per family. Each family module registers its
factories with ``@register_block`` at import; ``core.registry.ensure_loaded`` imports them all.

Families: source, noise, pol, mod, err, ctrl, act, det, dmx, loss, fsr, spf, los, fib, meas, chk.
Framework blocks live elsewhere: ``delay_buffer`` (core.delay), ``segment_phase`` (path.segments),
``beam_assembler`` (path.beams).
"""
