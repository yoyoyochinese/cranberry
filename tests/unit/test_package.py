"""Package skeleton (G15): every subpackage imports, the layout matches PLAN §4.1, CLI delegates."""

import importlib
import pkgutil
from pathlib import Path

import pytest

import cranberry
from cranberry import __main__ as cli

ROOT = Path(__file__).resolve().parents[2]

SUBPACKAGES = ["core", "config", "runner", "offline", "path", "blocks", "hardware"]
FAMILIES = ["source", "noise", "pol", "mod", "err", "ctrl", "act", "det", "dmx", "loss", "fsr", "spf", "los",
            "fib", "meas", "chk"]


def test_every_module_imports():
    names = [m.name for m in pkgutil.walk_packages(cranberry.__path__, prefix="cranberry.")]
    for sub in SUBPACKAGES:
        assert f"cranberry.{sub}" in names
    for fam in FAMILIES:
        assert f"cranberry.blocks.{fam}" in names
    for name in names:
        importlib.import_module(name)


def test_repository_layout():
    for d in ["presets/reference", "presets/future", "presets/axes", "presets/bad", "presets/locks", "devices",
              "tests/unit", "tests/vt", "tests/integration"]:
        assert (ROOT / d).is_dir(), d
    assert not any(p.name == "refs" for p in (ROOT / "src").rglob("*")), "no directory named refs (G15)"
    assert "runs/" in (ROOT / ".gitignore").read_text().splitlines()


def test_cli_validates_arguments(tmp_path, capsys):
    with pytest.raises(SystemExit):
        cli.main(["run", str(tmp_path / "missing.yaml")])
    txt = tmp_path / "x.txt"
    txt.write_text("a: 1\n")
    with pytest.raises(SystemExit):
        cli.main(["run", str(txt)])
    good = tmp_path / "x.yaml"
    good.write_text("a: 1\n")
    with pytest.raises(SystemExit):
        cli.main(["run", str(good), "--seed", "-3"])
    with pytest.raises(SystemExit):
        cli.main([])


@pytest.mark.parametrize(
    ("cmd", "wp"), [("ablate", "WP-F"), ("decompose", "WP-F4"), ("convergence", "WP-F4")]
)
def test_cli_delegates_and_names_wp(tmp_path, cmd, wp):
    # WP-F2: the CLI now loads and validates the config; an incomplete one exits with status 2
    bad = tmp_path / "refD.yaml"
    bad.write_text("preset: refD\n")
    assert cli.main([cmd, str(bad)]) == 2
    good = tmp_path / "ok.yaml"
    good.write_text("preset: ok\nrun: {duration_s: 1.0, f_sim_hz: 1.0e3, chunk_samples: 8, seed: 0}\npipeline: []\n")
    # run/sweep are implemented (WP-F1/F2) and tested in tests/unit/config/test_record_cli.py
    with pytest.raises(NotImplementedError, match=wp):
        cli.main([cmd, str(good)])
