"""``decompose_metric_ii(spec)``: four-component split of metric (ii) by reruns (G13). Owner: WP-F4.

Runs all-on plus one rerun per component noise-source set with the same RNG streams (G11), then checks
linearity: the component variances must sum to the all-on variance within 5 % (T04 V12). Saturation or
lock loss marks the decomposition invalid.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from cranberry.blocks.meas.metrics import MetricIIComponents
    from cranberry.config.spec import RunSpec

__all__ = ["decompose_metric_ii"]


def decompose_metric_ii(spec: RunSpec) -> dict[int, MetricIIComponents]:
    """{port k: components of Θ_k}; ports without metric (ii) are omitted. WP-F4."""
    raise NotImplementedError("runner.decompose_metric_ii is implemented by WP-F4")
