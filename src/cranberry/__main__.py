"""Command-line entry: ``python -m cranberry <command> <config.yaml> [options]`` (PLAN §4.1).

Commands: ``run``, ``sweep``, ``ablate``, ``decompose``, ``convergence``, ``diagram``. Arguments are validated here;
the config is loaded, overlaid and validated by ``cranberry.config`` (WP-F2; an invalid config prints
the reason and exits with status 2), then the work is delegated to the runner functions (WP-F1 /
WP-F4). Until those land, the delegated call raises ``NotImplementedError`` naming the responsible WP.

``run`` / ``sweep`` write run records (``runner.record``) only when ``--out DIR`` is given (e.g.
``--out runs/``): one directory per run, with the input as written (``config.input.yaml``), the
resolved config, seeds, git, env, results and ``invocation.json`` (argv, ``--seed`` override, sweep
point). There is no silent default directory (D34).

``ablate`` / ``decompose`` / ``convergence`` print their result as a table on stdout; ``--save FILE``
also writes it as JSON (``ablate``: ``ablation_impact`` plus each run's validity; ``decompose``: the
components per port; ``convergence``: the report's checks).

``diagram <preset or path> [--overlay REF ...] --out DIR`` writes the signal-flow and optical-layout
diagrams (``viz.diagram``, Graphviz ``dot``). The preset is a name (the unique ``presets/**/<name>.yaml``)
or a path; overlays are paths or ``<axis>/<value>`` references (``demux/none``), as in ``extends:`` /
``overlays:``.
"""

from __future__ import annotations

import argparse
import json
import math
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
        if name == "ablate":
            p.add_argument("--quantity", "-q", action="append", default=[], metavar="KEY",
                           help="RunResult.summary() key to tabulate (default: every key); repeatable")
        if name == "convergence":
            p.add_argument("--what", action="append", default=[], choices=["f_sim", "chunk", "seed", "n_sub"],
                           help="check to run (default: all); repeatable")
        if name in ("ablate", "decompose", "convergence"):
            p.add_argument("--save", type=Path, default=None, metavar="FILE",
                           help="also write the result as JSON to FILE (the table is always printed)")
    p = sub.add_parser("diagram", help="signal-flow and optical-layout diagrams (Graphviz)",
                       description="signal-flow and optical-layout diagrams of one configuration (viz.diagram)")
    p.add_argument("config", help="preset name (presets/**/<name>.yaml) or YAML path")
    p.add_argument("--overlay", "-o", action="append", default=[], metavar="REF",
                   help="overlay path or <axis>/<value> reference (e.g. demux/none); repeatable")
    p.add_argument("--out", type=Path, required=True, metavar="DIR", help="output directory")
    p.add_argument("--name", default=None, help="file stem (default: the preset file name, plus the overlays)")
    p.add_argument("--expand", action="append", default=[], choices=["noise", "bookkeeping", "all"],
                   help="draw a collapsed cluster block by block; repeatable")
    p.add_argument("--format", "-f", action="append", default=[], choices=["svg", "png", "pdf", "dot"],
                   help="output format (default: svg and png); repeatable")
    return parser


def _json_default(x):
    import dataclasses

    import numpy as np

    if dataclasses.is_dataclass(x):
        return dataclasses.asdict(x)
    if isinstance(x, (np.generic, np.ndarray)):
        return x.tolist()
    if isinstance(x, Path):
        return str(x)
    raise TypeError(f"not JSON serialisable: {type(x).__name__}")


def _save(path: Path | None, obj) -> None:
    if path is not None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(obj, indent=1, default=_json_default, allow_nan=True) + "\n")
        print(f"saved {path}")


def _report_ablate(results, quantities, save: Path | None) -> None:
    from cranberry.runner.ablate import ablation_impact

    impact = ablation_impact(results, quantities or None)
    first = quantities[0] if quantities else None

    def key(t):
        row = impact[t].get(first) if first else None
        return -abs(row["delta"]) if row and math.isfinite(row["delta"]) else 0.0

    print(f"{'toggle off':<40} {'quantity':<36} {'baseline':>12} {'ablated':>12} {'delta':>12}")
    for t in sorted(impact, key=key):
        for q, row in impact[t].items():
            print(f"{t:<40} {q:<36} {row['baseline']:>12.5g} {row['ablated']:>12.5g} {row['delta']:>+12.4g}")
    invalid = [k for k, r in results.items() if not r.valid]
    if invalid:
        print(f"invalid runs (saturation / lock loss / unstable loop): {invalid}")
    _save(save, {"impact": impact, "valid": {k: r.valid for k, r in results.items()}})


def _report_decompose(comps, save: Path | None) -> None:
    if not comps:
        print("no metric (ii) port in this configuration")
    for k, c in comps.items():
        print(f"port {k}: total {math.sqrt(max(c.total_var_rad2, 0.0)) * 1e3:.4g} mrad rms, "
              f"Σ/total {c.linearity_ratio:.4f}, {c.status}" + (f" ({c.reason})" if c.reason else ""))
        for name, v in c.variances_rad2.items():
            frac = v / c.total_var_rad2 if c.total_var_rad2 > 0 else float("nan")
            print(f"  {name:<26} {math.sqrt(max(v, 0.0)) * 1e3:>10.4g} mrad  {frac:>8.1%}")
    _save(save, {str(k): c for k, c in comps.items()})


def _report_convergence(report, save: Path | None) -> None:
    for name, c in report.checks.items():
        print(f"{name:<8} passed={c.get('passed')}  ({c.get('criterion')})" + (f"  error: {c['error']}" if "error" in c else ""))
        vals = c.get("values") or {}
        if name == "seed":
            for q, v in vals.items():
                if isinstance(v, dict) and "mean" in v:
                    print(f"  {q:<40} mean {v['mean']:.5g}  std {v['std']:.3g}")
        elif name in ("f_sim", "n_sub"):
            rows = vals.items() if name == "f_sim" else [("x2", vals)]
            for tag, cmp in rows:
                for q, r in cmp.items():
                    print(f"  {tag:<4} {q:<36} ref {r['reference']:.5g}  change {r['change']:+.3g}  tol {r['tolerance']:.3g}")
        elif name == "chunk":
            print(f"  {vals}")
    print(f"all passed: {report.passed}")
    _save(save, report.checks)


def _diagram(args) -> int:
    from cranberry import config
    from cranberry.config.overlay import find_base
    from cranberry.viz.diagram import write_diagrams

    try:
        path = Path(args.config) if Path(args.config).is_file() else find_base(args.config)
        spec = config.load(path, overlays=args.overlay)
    except config.ConfigError as e:
        print(f"cranberry: invalid configuration {args.config}: {e}", file=sys.stderr)
        return 2
    name = args.name or "+".join([path.stem, *(Path(o).with_suffix("").as_posix().replace("/", "-")
                                                                for o in args.overlay)])
    for p in write_diagrams(spec, args.out, name, expand=args.expand, formats=args.format or ("svg", "png")):
        print(p)
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    args = build_parser().parse_args(argv)
    if args.command == "diagram":
        return _diagram(args)

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

        _report_ablate(ablate(spec, args.toggle or list_toggles(spec)), args.quantity, args.save)
    elif args.command == "decompose":
        from cranberry.runner.decompose import decompose_metric_ii

        _report_decompose(decompose_metric_ii(spec), args.save)  # WP-F4
    elif args.command == "convergence":
        from cranberry.runner.convergence import CONVERGENCE_CHECKS, convergence

        _report_convergence(convergence(spec, args.what or CONVERGENCE_CHECKS), args.save)  # WP-F4
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
