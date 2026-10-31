"""Restricted ``!expr`` evaluator (G12). YAML is external input: the negative cases are security tests."""

import math

import numpy as np
import pytest

from cranberry.config.expr import ExprError, ExprFunction, compile_expr, free_names, safe_eval


@pytest.mark.parametrize(
    ("src", "env", "expected"),
    [
        ("1 + 2 * 3", {}, 7),
        ("2 ** 10", {}, 1024),
        ("-x + +x - (-1)", {"x": 4.0}, 1.0),
        ("7 / 2", {}, 3.5),
        ("0.00249 + 0.00222 * x", {"x": 2.0}, 0.00249 + 0.00222 * 2.0),
        ("arm1_915b.L_m - dL", {"arm1_915b.L_m": 8.4, "dL": 9.0}, 8.4 - 9.0),
        ("sqrt(x) + exp(0) + log(e) + log10(100)", {"x": 16.0}, 4.0 + 1.0 + 1.0 + 2.0),
        ("sin(pi / 2) + cos(0) + arctan2(1, 1) * 4 / pi", {}, 3.0),
        ("hypot(3, 4) + abs(-2) + min(1, 2) + max(1, 2)", {}, 10.0),
        ("floor(2.7) + ceil(2.1) + deg2rad(180) / pi", {}, 6.0),
        ("1 if x > 0 else -1", {"x": 0.5}, 1),
        ("1 if x > 0 else -1", {"x": -0.5}, -1),
        ("0 < x <= 1", {"x": 1.0}, True),
        ("0 < x <= 1", {"x": 2.0}, False),
        ("x == 'pm' and y != 3", {"x": "pm", "y": 4}, True),
        ("not (x or y)", {"x": False, "y": 0}, True),
        ("'slow' if a >= b else 'fast'", {"a": 1, "b": 2}, "fast"),
        ("pi", {"pi": 3.0}, 3.0),  # env shadows constants
    ],
)
def test_allowed_grammar(src, env, expected):
    got = safe_eval(src, env)
    assert got == pytest.approx(expected) if isinstance(expected, float) else got == expected
    assert type(got) in (int, float, bool, str)


def test_results_are_python_scalars_and_ints_stay_ints():
    assert type(safe_eval("2 * n", {"n": 4})) is int
    assert type(safe_eval("sqrt(4)", {})) is float


def test_free_names():
    assert free_names("a.b.c + sqrt(x) * pi - a.b.c") == ("a.b.c", "x", "pi")


@pytest.mark.parametrize(
    "src",
    [
        # imports / dunder escapes
        "__import__('os').system('true')",
        "().__class__.__bases__[0].__subclasses__()",
        "x.__class__",
        "_private",
        "a._b",
        # attribute access on non-names, arbitrary calls
        "(1).real",
        "'abc'.upper()",
        "x.y()",
        "eval('1')",
        "exec('1')",
        "open('/etc/passwd')",
        "getattr(x, 'y')",
        "globals()",
        "type(x)",
        "np.sqrt(4)",
        "math.sqrt(4)",
        "sqrt",  # unknown name when not called (no function objects as values)
        # calls with wrong form
        "sqrt(x=4)",
        "sqrt(*xs)",
        "sqrt(1, 2)",
        "min(1)",
        "sqrt(4)(1)",
        # containers, subscripts, comprehensions, lambdas
        "[1, 2]",
        "(1, 2)",
        "{1: 2}",
        "{1, 2}",
        "x[0]",
        "[i for i in range(3)]",
        "{i for i in x}",
        "(i for i in x)",
        "lambda: 1",
        # other statements / expressions
        "x := 1",
        "f'{x}'",
        "1 % 2",
        "7 // 2",
        "1 << 2",
        "1 & 2",
        "~1",
        "x is None",
        "1 in x",
        "None",
        "b'abc'",
        "1j",
        "...",
        "import os",
        "x = 1",
        "yield 1",
        "await x",
        "",
        "   ",
        "1 +",
    ],
)
def test_rejected_constructs(src):
    with pytest.raises(ExprError):
        safe_eval(src, {"x": 1.0, "xs": [1.0], "y": 2.0, "np": np, "math": math})


def test_no_python_attribute_lookup_on_env_values():
    class Evil:
        @property
        def boom(self):
            raise AssertionError("attribute was looked up")

    # "x.boom" is a dotted config reference, looked up as the key "x.boom", never getattr(x, "boom")
    with pytest.raises(ExprError, match="unknown name 'x.boom'"):
        safe_eval("x.boom", {"x": Evil()})


@pytest.mark.parametrize(
    "src",
    [
        "9 ** 9 ** 9",  # huge integer power
        "2 ** 5000",
        "(10 ** 1000) ** 5",  # integer longer than the bound
        "'a' * 1000000000",  # no string repetition
        "'a' + 'b'",
        "-'a'",
        "sqrt('a')",
        "1 / 0",
        "1.0 / 0",
        "log(0)",
        "sqrt(-1)",
        "exp(1000)",
        "10.0 ** 400",
        "(-8) ** (1 / 3)",  # complex result
        "x < 'a'",
        "undefined_name + 1",
    ],
)
def test_resource_and_evaluation_errors(src):
    with pytest.raises(ExprError):
        safe_eval(src, {"x": 1.0})


def test_size_limits():
    with pytest.raises(ExprError, match="longer than"):
        safe_eval("1+" * 1500 + "1", {})
    with pytest.raises(ExprError):
        safe_eval("-" * 100 + "1", {})  # nesting depth
    with pytest.raises(ExprError):
        safe_eval("+".join(["x"] * 400), {"x": 1})


def test_compile_expr_callable():
    f = compile_expr("a * x + b", ("x",), {"a": 2.0, "b": 1.0})
    assert isinstance(f, ExprFunction)
    assert f(3.0) == 7.0
    np.testing.assert_array_equal(f(np.array([0.0, 1.0])), [1.0, 3.0])  # works on grids (tables)
    assert f == compile_expr("a * x + b", ("x",), {"a": 2.0, "b": 1.0})
    assert f != compile_expr("a * x + b", ("x",), {"a": 3.0, "b": 1.0})
    with pytest.raises(TypeError):
        f(1.0, 2.0)
    with pytest.raises(ExprError, match="unknown name"):
        compile_expr("a * x + c", ("x",), {"a": 1.0})
    with pytest.raises(ExprError):
        compile_expr("x.__class__", ("x",))
    with pytest.raises(ExprError):
        compile_expr("x", ("__x",))
