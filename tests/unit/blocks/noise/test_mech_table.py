"""``segment_mech_noise`` with ``mech_shape: table`` (D85-1): the OU-sum synthesis follows the tabulated PSD
within the documented tolerance, the ledger rows are the table, the time loop reproduces the table, and the
RNG sub-streams stay deterministic and independent (G11)."""

import numpy as np
import pytest
import scipy.signal

import cranberry.blocks.noise.testing as nt
from cranberry.blocks.err.testing import record, rig
from cranberry.blocks.noise._lgs import loglog_interp, section_psd, table_fit

pytestmark = pytest.mark.numba
TOPO = nt.seg_topo({"a": (5.0, 1, "uniform"), "b": (9.0, 4, "per_sub")})
# presets/axes/lab/acoustic.yaml shape (Ref D free-running phase PSD, COMPARISON.md §3 A, + its extrapolation)
REFD = [[10.0, 6.32e-4], [100.0, 2.0e-5], [300.0, 3.6e-6], [1000.0, 6.7e-7], [3000.0, 9.0e-8], [3.0e4, 9.0e-10]]
L_REF = 13.2


def table_mech(segments=("a",), on=("mech",), **over):
    p = {"mech_shape": "table", "psd_table": REFD, "table_f_min_hz": 10.0, "table_f_max_hz": 3.0e4,
         "table_per_decade": 4, "L_ref_m": L_REF}
    p.update(over)
    return nt.mech(TOPO, list(segments), on=on, **p)


def test_loglog_interp():
    s = loglog_interp(REFD)
    np.testing.assert_allclose(s([100.0, 300.0, 3.0e4]), [2e-5, 3.6e-6, 9e-10], rtol=1e-12)
    np.testing.assert_allclose(s(np.sqrt(100.0 * 300.0)), np.sqrt(2e-5 * 3.6e-6), rtol=1e-12)  # log-log midpoint
    np.testing.assert_allclose(s(3.0e5), 9e-12, rtol=1e-9)  # last segment's f⁻² extended
    with pytest.raises(ValueError):
        loglog_interp([[100.0, 1e-5], [100.0, 2e-6]])


@pytest.mark.parametrize("per_decade, tol", [(3, 0.08), (4, 0.06), (6, 0.06)])
def test_synthesis_matches_table_in_band(per_decade, tol):
    """Documented tolerance (mech.py): ≤ 10 % on [3 f_min, f_max/1.5]; this table: ≤ 6–8 %."""
    s = loglog_interp(REFD)
    secs = table_fit(s, 10.0, 3.0e4, per_decade)
    f = np.geomspace(30.0, 2.0e4, 400)
    rel = sum(section_psd(x, f) for x in secs) / s(f) - 1.0
    assert np.abs(rel).max() < tol
    # above f_max the OU roll-off is f⁻², the table's last segment: the two agree (ledger continuity)
    fh = np.geomspace(5.0e4, 5.0e5, 50)
    np.testing.assert_allclose(sum(section_psd(x, fh) for x in secs) / s(fh), 1.0, atol=0.06)


def test_ledger_rows_are_the_table():
    root, _, _ = rig([table_mech(segments=["a", "b"], spatial_coherence={"b": 0.3})], 1e5, external={})
    decls = {d.targets: d for d in root.bands()}
    f = np.geomspace(5.0e4, 1e6, 7)
    for seg, L in (("a", 5.0), ("b", 9.0)):
        d = decls[(seg,)]
        assert d.quantity == "strain" and d.f_lo_hz == 0.0 and d.f_hi_hz == pytest.approx(5.0e4)
        np.testing.assert_allclose(d.psd_outside(f) * (nt.K_EPS * L) ** 2, (L / L_REF) * loglog_interp(REFD)(f),
                                   rtol=1e-12)


def test_time_loop_reproduces_table():
    fs, n = 1.0e5, 1 << 19
    x = record([table_mech()], fs, n, ["drv.a.strain"], seed=11, external={})["drv.a.strain"]
    phase = x * nt.K_EPS * 5.0 / np.sqrt(5.0 / L_REF)  # back to the table's reference length (γ = 1)
    f, p = scipy.signal.welch(phase, fs, nperseg=1 << 14)
    for f0 in (100.0, 300.0, 1000.0, 3000.0, 1.0e4):
        band = (f > f0 / 1.15) & (f < f0 * 1.15)
        assert np.mean(p[band] / loglog_interp(REFD)(f[band])) == pytest.approx(1.0, abs=0.15), f0


def test_gamma_zero_is_length_independent():
    root, _, _ = rig([table_mech(segments=["a", "b"], gamma=0.0, spatial_coherence={"b": 1.0})], 1e5, external={})
    f = np.array([6e4])
    d = {x.targets: x.psd_outside(f)[0] for x in root.bands()}
    assert d[("a",)] * (nt.K_EPS * 5.0) ** 2 == pytest.approx(d[("b",)] * (nt.K_EPS * 9.0) ** 2, rel=1e-12)


def test_rng_substreams_deterministic():
    n, s = 4000, ["drv.a.strain"]
    a = record([table_mech()], 1e5, n, s, seed=5, external={})["drv.a.strain"]
    b = record([table_mech()], 1e5, n, s, seed=5, external={})["drv.a.strain"]
    c = record([table_mech()], 1e5, n, s, seed=6, external={})["drv.a.strain"]
    np.testing.assert_array_equal(a, b)
    assert np.abs(a - c).max() > 0.5 * np.std(a)
    # the peaks bank is a separate sub-stream: switching it on adds its own part, the table part is unchanged
    both = record([table_mech(on=("mech", "peaks"))], 1e5, n, s, seed=5, external={})["drv.a.strain"]
    pk = record([table_mech(on=("peaks",))], 1e5, n, s, seed=5, external={})["drv.a.strain"]
    np.testing.assert_allclose(both, a + pk, rtol=1e-12, atol=1e-25)


def test_two_pole_explicit_equals_absent():
    n, s = 3000, ["drv.a.strain"]
    absent = record([nt.mech(TOPO, ["a"])], 1e5, n, s, seed=2, external={})["drv.a.strain"]
    explicit = record([nt.mech(TOPO, ["a"], mech_shape="two_pole")], 1e5, n, s, seed=2, external={})["drv.a.strain"]
    np.testing.assert_array_equal(absent, explicit)
    with pytest.raises(ValueError):
        rig([nt.mech(TOPO, ["a"], mech_shape="spline")], 1e5, external={})
