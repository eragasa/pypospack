"""Bounded arithmetic for historical PyPosPack parameter expressions."""

from __future__ import annotations

import ast
import math
import operator
from collections.abc import Callable, Mapping
from numbers import Real

_MAX_EXPRESSION_CHARACTERS = 512
_MAX_AST_NODES = 64
_MAX_ABS_EXPONENT = 64.0

_BINARY_OPERATORS: dict[type[ast.operator], Callable[[float, float], float]] = {
    ast.Add: operator.add,
    ast.Sub: operator.sub,
    ast.Mult: operator.mul,
    ast.Div: operator.truediv,
    ast.Pow: operator.pow,
}
_UNARY_OPERATORS: dict[type[ast.unaryop], Callable[[float], float]] = {
    ast.UAdd: operator.pos,
    ast.USub: operator.neg,
}
_COMPARISON_OPERATORS: dict[type[ast.cmpop], Callable[[float, float], bool]] = {
    ast.Eq: operator.eq,
    ast.NotEq: operator.ne,
    ast.Lt: operator.lt,
    ast.LtE: operator.le,
    ast.Gt: operator.gt,
    ast.GtE: operator.ge,
}


class ParameterExpressionError(ValueError):
    """Raised when an expression violates the bounded arithmetic contract."""


def evaluate_numeric_expression(
    expression: str,
    parameters: Mapping[str, object],
) -> float:
    """Evaluate closed arithmetic over explicitly supplied finite parameters."""
    value = _evaluate(expression, parameters)
    if isinstance(value, bool):
        raise ParameterExpressionError("numeric expression returned a boolean")
    return _finite_number(value, label="expression result")


def evaluate_constraint_expression(
    expression: str,
    parameters: Mapping[str, object],
) -> bool:
    """Evaluate a closed comparison over explicitly supplied parameters."""
    value = _evaluate(expression, parameters)
    if not isinstance(value, bool):
        raise ParameterExpressionError("constraint expression must return bool")
    return value


def _evaluate(expression: str, parameters: Mapping[str, object]) -> float | bool:
    if not isinstance(expression, str) or not expression.strip():
        raise ParameterExpressionError("expression must be a non-empty string")
    if len(expression) > _MAX_EXPRESSION_CHARACTERS:
        raise ParameterExpressionError("expression exceeds its character limit")
    try:
        tree = ast.parse(expression, mode="eval")
    except SyntaxError as error:
        raise ParameterExpressionError("expression syntax is invalid") from error
    if sum(1 for _node in ast.walk(tree)) > _MAX_AST_NODES:
        raise ParameterExpressionError("expression exceeds its AST node limit")
    return _evaluate_node(tree.body, parameters)


def _evaluate_node(
    node: ast.expr,
    parameters: Mapping[str, object],
) -> float | bool:
    if isinstance(node, ast.Constant):
        return _finite_number(node.value, label="numeric literal")
    if isinstance(node, ast.Name):
        try:
            value = parameters[node.id]
        except KeyError as error:
            raise ParameterExpressionError(
                f"unknown parameter name: {node.id}"
            ) from error
        return _finite_number(value, label=f"parameter {node.id!r}")
    if isinstance(node, ast.BinOp):
        operation = _BINARY_OPERATORS.get(type(node.op))
        if operation is None:
            raise ParameterExpressionError("binary operator is not allowed")
        left = _finite_number(
            _evaluate_node(node.left, parameters),
            label="left operand",
        )
        right = _finite_number(
            _evaluate_node(node.right, parameters),
            label="right operand",
        )
        if isinstance(node.op, ast.Pow) and abs(right) > _MAX_ABS_EXPONENT:
            raise ParameterExpressionError("exponent exceeds its fixed limit")
        try:
            result = operation(left, right)
        except (ArithmeticError, OverflowError) as error:
            raise ParameterExpressionError("arithmetic operation failed") from error
        return _finite_number(result, label="arithmetic result")
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.Not):
        value = _evaluate_node(node.operand, parameters)
        if not isinstance(value, bool):
            raise ParameterExpressionError("not requires a comparison operand")
        return not value
    if isinstance(node, ast.UnaryOp):
        operation = _UNARY_OPERATORS.get(type(node.op))
        if operation is None:
            raise ParameterExpressionError("unary operator is not allowed")
        operand = _finite_number(
            _evaluate_node(node.operand, parameters),
            label="unary operand",
        )
        return _finite_number(operation(operand), label="unary result")
    if isinstance(node, ast.Compare):
        left = _finite_number(
            _evaluate_node(node.left, parameters),
            label="comparison operand",
        )
        for operator_node, comparator in zip(
            node.ops,
            node.comparators,
            strict=True,
        ):
            operation = _COMPARISON_OPERATORS.get(type(operator_node))
            if operation is None:
                raise ParameterExpressionError("comparison operator is not allowed")
            right = _finite_number(
                _evaluate_node(comparator, parameters),
                label="comparison operand",
            )
            if not operation(left, right):
                return False
            left = right
        return True
    if isinstance(node, ast.BoolOp):
        values = [_evaluate_node(value, parameters) for value in node.values]
        if not all(isinstance(value, bool) for value in values):
            raise ParameterExpressionError(
                "boolean operators require comparison operands"
            )
        if isinstance(node.op, ast.And):
            return all(values)
        if isinstance(node.op, ast.Or):
            return any(values)
        raise ParameterExpressionError("boolean operator is not allowed")
    raise ParameterExpressionError(
        f"expression construct is not allowed: {type(node).__name__}"
    )


def _finite_number(value: object, *, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ParameterExpressionError(f"{label} must be a real number")
    normalized = float(value)
    if not math.isfinite(normalized):
        raise ParameterExpressionError(f"{label} must be finite")
    return normalized
