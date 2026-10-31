"""Shared pytest setup. Markers are registered in pyproject.toml ([tool.pytest.ini_options]).

Tests under tests/vt/ are auto-marked ``vt`` and tests under tests/integration/ ``integration``, so
``pytest -m "not vt and not integration"`` runs the fast unit tests only.
"""

from pathlib import Path

import pytest

_TESTS = Path(__file__).parent


def pytest_collection_modifyitems(config, items):
    for item in items:
        rel = Path(item.path).resolve().relative_to(_TESTS.resolve())
        top = rel.parts[0] if rel.parts else ""
        if top == "vt":
            item.add_marker(pytest.mark.vt)
        elif top == "integration":
            item.add_marker(pytest.mark.integration)
