from __future__ import annotations

import unittest

from pypospack.expressions import (
    ParameterExpressionError,
    evaluate_constraint_expression,
    evaluate_numeric_expression,
)


class ParameterExpressionTest(unittest.TestCase):
    def test_evaluates_closed_arithmetic_with_explicit_names(self) -> None:
        self.assertEqual(
            evaluate_numeric_expression(
                "-chrg_Mg + 2 * offset",
                {"chrg_Mg": 2.0, "offset": 0.25},
            ),
            -1.5,
        )
        self.assertEqual(
            evaluate_numeric_expression(
                "a + gamma",
                {"a": 1.0, "gamma": 2.0, "latticetype": "fcc"},
            ),
            3.0,
        )

    def test_evaluates_comparisons_and_boolean_composition(self) -> None:
        self.assertTrue(
            evaluate_constraint_expression(
                "left < value and value <= right",
                {"left": 0.0, "value": 1.0, "right": 1.0},
            )
        )

    def test_rejects_code_execution_constructs(self) -> None:
        for expression in (
            "__import__('os').system('true')",
            "value.__class__",
            "[value][0]",
            "lambda: value",
        ):
            with self.subTest(expression=expression):
                with self.assertRaisesRegex(
                    ParameterExpressionError,
                    "not allowed",
                ):
                    evaluate_numeric_expression(expression, {"value": 1.0})

    def test_rejects_unknown_names_and_non_finite_values(self) -> None:
        with self.assertRaisesRegex(
            ParameterExpressionError,
            "unknown parameter name",
        ):
            evaluate_numeric_expression("missing + 1", {})
        with self.assertRaisesRegex(ParameterExpressionError, "must be finite"):
            evaluate_numeric_expression("value", {"value": float("inf")})

    def test_rejects_unbounded_exponents(self) -> None:
        with self.assertRaisesRegex(
            ParameterExpressionError,
            "exponent exceeds",
        ):
            evaluate_numeric_expression("2 ** 65", {})


if __name__ == "__main__":
    unittest.main()
