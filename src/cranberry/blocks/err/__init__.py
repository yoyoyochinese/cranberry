"""ERR: error-signal generators #1-#6, demodulator, RF reference (PLAN §4.3 ERR; baseband API at loop rate, G16-b).

Owner: WP-B3. Notes: T14 (``notes/error-signal-generation.md``), T15 (DC fringe superset).

Block kinds (factories ``make_<kind>``, registered with ``@register_block``):

=============================  ==========================================================  ===============
kind                           error signal                                                module
=============================  ==========================================================  ===============
``parametric_beat_detector``   #1a/#1b Ref A lock 1 (2 f_sh), tap variant, #2 Ref B lock 1   ``coherent``
``bhd_probe_beat``             #1c Ref A lock 2 (f_sh), #3 Ref B lock 2 (dither), DC        ``coherent``
``liu_probe_interferometer``   #6 Liu-type counter-propagating probe vs reference           ``coherent``
``quantum_noise_power``        #4 Ref C quantum noise lock + radiometer noise               ``qnl``
``dc_fringe``                  #5 DC fringe (T15 eq. (3) superset, Ref D)                    ``dc_fringe``
``demodulator``                I/Q demodulation + LPF + offset                              ``demod``
``rf_reference``               RF phase noise of AOM drives / demod references              ``demod``
=============================  ==========================================================  ===============

Signal representation (baseband vs ``rf_rate_mode``), sensor noise and dither harmonics: ``_common``.
Contract: core.block (factory interface, reads/writes, losses(), bands(), toggles, D34 no defaults).
"""

from cranberry.blocks.err import coherent, dc_fringe, demod, qnl  # noqa: F401  (registers the kinds)
from cranberry.blocks.err._common import amplitude_gains_from_table, amplitude_gains_opa, amplitude_gains_opo
from cranberry.blocks.err.dc_fringe import fringe_lock_point, gamma_coh_from_ledger
from cranberry.blocks.err.qnl import qnl_band_integrals, radiometer_mc

__all__ = [
    "amplitude_gains_opa",
    "amplitude_gains_opo",
    "amplitude_gains_from_table",
    "fringe_lock_point",
    "gamma_coh_from_ledger",
    "qnl_band_integrals",
    "radiometer_mc",
]
