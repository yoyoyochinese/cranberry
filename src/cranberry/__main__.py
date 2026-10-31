"""Command-line entry: ``python -m cranberry <command> <config.yaml> [options]`` (PLAN §4.1).

Commands: ``run``, ``sweep``, ``ablate``, ``decompose``, ``convergence``. Arguments are validated here;
the config is loaded, overlaid and validated by ``cranberry.config`` (WP-F2; an invalid config prints
the reason and exits with status 2), then the work is delegated to the runner functions (WP-F1 /
WP-F4). Until those land, the delegated call raises ``NotImplementedError`` naming the responsible WP.

``run`` / ``sweep`` write run records (``runner.record``) only when ``--out DIR`` is given (e.g.
``--out runs/``): one directory per run, with the input as written (``config.input.yaml``), the
resolved config, seeds, git, env, results and ``invocation.json`` (argv, ``--seed`` override, sweep
point). There is no silent default directory (D34).
"""

from __future__ import annotations

import argparse
import sys
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
        if name in ("run", "sweep"):
            p.add_argument("--out", type=Path, default=None, metavar="DIR",
                           help="write a run record per run under DIR (e.g. runs/); without it nothing is written")
        if name == "ablate":
            p.add_argument("--toggle", action="append", default=[], metavar="BLOCK.TOGGLE",
                           help="toggle to switch off (default: all toggles); repeatable")
        if name == "convergence":
            p.add_argument("--what", action="append", default=[], choices=["f_sim", "chunk", "seed", "n_sub"],
                           help="check to run (default: all); repeatable")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    args = build_parser().parse_args(argv)

    from dataclasses import replace

    from cranberry import config

    try:
        spec = config.load(args.config, overlays=args.overlay)  # WP-F2
        if args.seed is not None:
            spec = replace(spec, run=replace(spec.run, seed=args.seed))
        if args.command == "sweep":
            from cranberry.runner.sweep import expand

            points = expand(spec)
            for p in points:
                config.validate(_as_run_validates(p.spec))
        else:
            config.validate(_as_run_validates(spec))
    except config.ConfigError as e:
        print(f"cranberry: invalid configuration {args.config}: {e}", file=sys.stderr)
        return 2

    record: dict = {}
    if getattr(args, "out", None) is not None:  # D34: records only where asked (no silent default dir)
        record = {
            "record_root": args.out,
            "config_input": config.read_input(args.config, args.overlay),
            "invocation": {"argv": argv, "overrides": {} if args.seed is None else {"run.seed": args.seed}},
        }
    if args.command == "run":
        from cranberry.runner.run import run

        run(spec, **record)  # WP-F1
    elif args.command == "sweep":
        from cranberry.runner.sweep import sweep

        sweep(spec, **record)  # each point via runner.run (WP-F1)
    elif args.command == "ablate":
        from cranberry.runner.ablate import ablate, list_toggles

        ablate(spec, args.toggle or list_toggles(spec))
    elif args.command == "decompose":
        from cranberry.runner.decompose import decompose_metric_ii

        decompose_metric_ii(spec)  # WP-F4
    elif args.command == "convergence":
        from cranberry.runner.convergence import CONVERGENCE_CHECKS, convergence

        convergence(spec, args.what or CONVERGENCE_CHECKS)  # WP-F4
    return 0


def _as_run_validates(spec):
    """The spec as ``runner.run`` validates it: ``RunSpec.topology`` injected into the topology-taking
    leaves and toggle overrides of ``offline.blocks`` set aside (WP-I1 fix: before it, every preset
    relying on topology injection, e.g. refA, was rejected by the CLI although ``run`` accepts it)."""
    from dataclasses import replace

    from cranberry.runner.run import inject_topology, offline_block_specs

    off = {b.name for b in offline_block_specs(spec)}
    toggles = {k: v for k, v in spec.toggles.items() if k.rpartition(".")[0] not in off}
    return replace(spec, pipeline=inject_topology(spec.pipeline, spec.topology), toggles=toggles)


if __name__ == "__main__":
    raise SystemExit(main())
