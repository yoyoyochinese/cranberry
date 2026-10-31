"""Test-only blocks for the config tests (registered under ``cfgtest_*`` kinds) and a reference loop.

``cfgtest_noise``: out = sigma · N(0, 1) from the block's own RNG stream (toggle ``enabled``).
``cfgtest_shaper``: out = curve(in), ``curve`` a callable parameter (number, !expr, !fn, Python fn).
"""

from types import SimpleNamespace

import numpy as np
import pytest

from cranberry.config.resolve import as_callable
from cranberry.config.spec import apply_toggles
from cranberry.core import registry
from cranberry.core.block import BaseBlock, require
from cranberry.core.layout import SlotLayout
from cranberry.core.pipeline import from_spec, iter_leaves
from cranberry.core.rng import stream
from cranberry.core.state import Params, State


class Noise(BaseBlock):
    toggle_names = ("enabled",)

    def __init__(self, params, name):
        require(params, "out", block=name or "cfgtest_noise")
        super().__init__("cfgtest_noise", name, params, writes=(params["out"],), runtime_params=("sigma",))

    def init(self, layout, rng):
        super().init(layout, rng)
        layout.add(self.params["out"])

    def __call__(self, state, params):
        x = self.rng.standard_normal()
        state[self.params["out"]] = float(params[self.param_name("sigma")]) * x if self.toggles["enabled"] else 0.0
        return state


class Shaper(BaseBlock):
    def __init__(self, params, name):
        require(params, "in", "out", "curve", block=name or "cfgtest_shaper")
        super().__init__("cfgtest_shaper", name, params, reads=(params["in"],), writes=(params["out"],))
        self.curve = as_callable(params["curve"], ("x",))

    def init(self, layout, rng):
        super().init(layout, rng)
        layout.add(self.params["out"])

    def __call__(self, state, params):
        state[self.params["out"]] = self.curve(float(state[self.params["in"]]))
        return state


@registry.register_block("cfgtest_noise")
def make_cfgtest_noise(params, *, name=None):
    return Noise(params, name)


@registry.register_block("cfgtest_shaper")
def make_cfgtest_shaper(params, *, name=None):
    return Shaper(params, name)


@registry.register_fn("cfgtest_affine")
def cfgtest_affine(x, a=1.0, b=0.0):
    return a * x + b


def simulate_block(block, f_sim_hz, seed, n_samples):
    """Reference loop over a built block tree: {written slot: time series}."""
    lay = SlotLayout(f_sim_hz)
    for leaf in iter_leaves(block):
        leaf.init(lay, stream(seed, leaf.name))
    lay.freeze()
    block.bind(lay)
    st, p = State(lay), Params(lay)
    out = {s: np.empty(n_samples) for s in block.writes}
    for n in range(n_samples):
        st["n"] = n
        st["t"] = n * lay.dt_s
        block(st, p)
        for s in out:
            out[s][n] = float(st[s])
    return out


def simulate(spec):
    block = from_spec(apply_toggles(spec.pipeline, spec.toggles))
    return simulate_block(block, spec.run.f_sim_hz, spec.run.seed, spec.run.n_samples)


@pytest.fixture
def blocks():
    return SimpleNamespace(
        noise=make_cfgtest_noise, shaper=make_cfgtest_shaper, affine=cfgtest_affine,
        simulate=simulate, simulate_block=simulate_block,
    )
