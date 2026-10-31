"""Per-path phase model (protected list; PLAN §2.3, conv §9): Topology, segments/subsegments/drivers,
``segment_phase``, ``beam_assembler`` (emission time, counter-propagation), locks, coefficients.
Owner: WP-F6.

Typical pipeline order within one sample (wiring by slot name)::

    [noise / TEC / actuator blocks → drv.<seg>.<driver>]  [laser_source → laser.<s>.phi]
    → segment_phase (psi.<seg>.<beam>, tau_g.<seg>.<beam>)
    → beam_assembler (beam.<b>.phi, .power, .jones, .t_emit)
    → lock_error (lock.<id>.err) / error-signal blocks → controllers → delay_buffer → actuators

Importing the package registers the block kinds ``segment_phase``, ``beam_assembler`` and
``lock_error`` (``core.registry.ensure_loaded`` imports ``cranberry.path.segments``, which runs this
``__init__`` first).
"""

from cranberry.path import beams, locks, segments  # noqa: F401  (register block kinds)
from cranberry.path.topology import Topology, as_topology  # noqa: F401
