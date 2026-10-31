"""Command-line entry: ``python -m cranberry <command> <config.yaml> [options]`` (PLAN §4.1).

Commands: ``run``, ``sweep``, ``ablate``, ``decompose``, ``convergence``. Arguments are validated here;
the work is delegated to ``cranberry.config.load`` (WP-F2) and the runner functions (WP-F1 / WP-F4).
Until those land, the delegated call raises ``NotImplementedError`` naming the responsible WP.
"""

from __future__ import annotations

import argparse
from collections.abc import Sequence
from pathlib import Path

__all__ = ["main", "build_parser"]

_YAML_SUFFIXES = (".yaml", ".yml")


def _yaml_file(value: str) -> Path:
    p = Path(value)
    if p.suffix.lower() not in _YAML_SUFFIXES:
        raise argparse.ArgumentTypeError(f"{value!r}: expected a .yaml/.yml file")
    if not p.is_file():
        raise argparse.ArgumentTypeError(f"{value!r}: file not found")
    return p


def _seed(value: str) -> int:
    try:
        s = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"{value!r}: seed must be an integer") from None
    if s < 0:
        raise argparse.ArgumentTypeError("seed must be >= 0")
    return s


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="cranberry", description="Squeezed-light phase-stabilization simulator")
    sub = parser.add_subparsers(dest="command", required=True)
    helps = {
        "run": "run one configuration",
        "sweep": "run the sweep axes (!sweep/!sweep_log/!choice) of a configuration",
        "ablate": "switch toggles off one at a time",
        "decompose": "four-component decomposition of metric (ii) (G13)",
        "convergence": "numerical convergence checks (PLAN §6.4)",
    }
    for name, text in helps.items():
        p = sub.add_parser(name, help=text, description=text)
        p.add_argument("config", type=_yaml_file, help="preset / config YAML file")
        p.add_argument("--overlay", "-o", type=_yaml_file, action="append", default=[], metavar="YAML",
                       help="overlay YAML applied in order (G12); repeatable")
        p.add_argument("--seed", type=_seed, default=None, help="override run.seed (root of all RNG streams)")
        if name == "ablate":
            p.add_argument("--toggle", action="append", default=[], metavar="BLOCK.TOGGLE",
                           help="toggle to switch off (default: all toggles); repeatable")
        if name == "convergence":
            p.add_argument("--what", action="append", default=[], choices=["f_sim", "chunk", "seed", "n_sub"],
                           help="check to run (default: all); repeatable")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    from cranberry import config

    spec = config.load(args.config, overlays=args.overlay)  # WP-F2
    if args.seed is not None:
        from dataclasses import replace

        spec = replace(spec, run=replace(spec.run, seed=args.seed))

    if args.command == "run":
        from cranberry.runner.run import run

        run(spec)  # WP-F1
    elif args.command == "sweep":
        from cranberry.runner.sweep import sweep

        sweep(spec)
    elif args.command == "ablate":
        from cranberry.runner.ablate import ablate

        ablate(spec, args.toggle or sorted(spec.toggles))
    elif args.command == "decompose":
        from cranberry.runner.decompose import decompose_metric_ii

        decompose_metric_ii(spec)  # WP-F4
    elif args.command == "convergence":
        from cranberry.runner.convergence import CONVERGENCE_CHECKS, convergence

        convergence(spec, args.what or CONVERGENCE_CHECKS)  # WP-F4
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
