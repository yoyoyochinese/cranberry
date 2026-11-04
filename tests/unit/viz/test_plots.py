"""``viz.plots``: the Phase 6 report helpers draw what they are given (D85-13)."""

import numpy as np
import pytest

mpl = pytest.importorskip("matplotlib")
mpl.use("Agg")

from cranberry.runner.explore import sensitivity_s2, servo_bump  # noqa: E402
from cranberry.viz import plots  # noqa: E402


def test_psd_with_bands_and_sensitivity(tmp_path):
    f = np.logspace(0, 5, 200)
    ax = plots.plot_psd({"free": (f, 1 / f**2), "locked": (f, 1e-3 / f)}, bands=[(10, 1e3, "loop band")])
    assert len(ax.get_lines()) == 2 and len(ax.patches) == 1
    ax = plots.plot_sensitivity({"nominal": (f, np.ones_like(f))}, limit=2.0)
    assert any("limit 2" in t.get_text() for t in ax.get_legend().get_texts())
    ax.figure.savefig(tmp_path / "s.png")
    assert (tmp_path / "s.png").stat().st_size > 0


def test_component_bars_and_length_plot():
    var = {"liu": {"locked_path_residual": 4e-6, "common_points": 1e-6, "laser_decorrelation_oob": 9e-6,
                   "squeezed_only": 0.0},
           "p2": {"locked_path_residual": 1e-6, "common_points": 1e-6, "laser_decorrelation_oob": 9e-6,
                  "squeezed_only": 1e-6}}
    ax = plots.plot_metric_ii_components(var, spread_rad2={"liu": {"locked_path_residual": 1e-6}})
    heights = [p.get_height() for p in ax.patches]
    assert len(heights) == 8 and heights[0] == pytest.approx(2.0) and heights[2] == pytest.approx(3.0)  # mrad
    L = [1, 10, 100]
    ax = plots.plot_vs_length(L, {"a": np.array([[1, 2], [3, 4], [5, 6]]), "b": ([1, 2, 3], [0.1, 0.1, 0.1])},
                              ylabel="ii [mrad]")
    assert np.allclose(ax.get_lines()[0].get_ydata(), [1.5, 3.5, 5.5])
    with pytest.raises(ValueError):
        plots.plot_vs_length(L, {"bad": np.zeros((2, 2))}, ylabel="x")


def test_sensitivity_estimator_recovers_a_known_filter():
    rng = np.random.default_rng(0)
    x = rng.standard_normal(1 << 15)
    y = 0.5 * x  # |S|² = 0.25 at every frequency
    f, s2 = sensitivity_s2(y, x, 1e4, nperseg=256)
    assert np.allclose(s2, 0.25)
    peak, fp = servo_bump(f, s2 * (1 + (np.abs(f - 2e3) < 30)), 100, 4e3)
    assert peak == pytest.approx(0.5) and abs(fp - 2e3) < 40
    with pytest.raises(ValueError):
        servo_bump(f, s2, 1e9)
