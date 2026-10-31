"""cranberry: modular simulation of phase-stabilized squeezed light (OPA/OPO → fiber → BHD / GBS).

Two stages (prompt §5.1, PLAN §2): an offline stage (``cranberry.offline``, GNLSE / quantum layer /
tables) and a time-loop stage (``cranberry.core`` pipeline of blocks over the per-path phase model
``cranberry.path``), driven by the runner (``cranberry.runner``) from a ``RunSpec``
(``cranberry.config``). Conventions: notes/conventions.md. Decisions: DECISIONS.md.

Importing the package is cheap; subpackages are imported on demand.
"""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("cranberry")
except PackageNotFoundError:  # running from a source tree without installation
    __version__ = "0+unknown"

__all__ = ["__version__"]
