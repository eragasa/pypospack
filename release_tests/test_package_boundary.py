from __future__ import annotations

import ast
import re
import tomllib
import unittest
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
PACKAGE_ROOT = REPOSITORY_ROOT / "pypospack"
_SHA256_LINE = re.compile(r"[0-9a-f]{64}  ([^\r\n]+)")


class PackageBoundaryTest(unittest.TestCase):
    def test_all_packaged_python_modules_parse(self) -> None:
        failures = []
        for path in sorted(PACKAGE_ROOT.rglob("*.py")):
            try:
                ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
            except SyntaxError as error:
                failures.append(
                    f"{path.relative_to(PACKAGE_ROOT)}:{error.lineno}: {error.msg}"
                )
        self.assertEqual(failures, [])

    def test_package_tree_contains_python_source_only(self) -> None:
        unexpected = sorted(
            str(path.relative_to(PACKAGE_ROOT))
            for path in PACKAGE_ROOT.rglob("*")
            if path.is_file()
            and path.suffix != ".py"
            and "__pycache__" not in path.parts
        )
        self.assertEqual(unexpected, [])

    def test_build_metadata_defines_explicit_boundaries(self) -> None:
        configuration = tomllib.loads(
            (REPOSITORY_ROOT / "pyproject.toml").read_text(encoding="utf-8")
        )
        project = configuration["project"]
        self.assertEqual(project["version"], "0.1.0")
        self.assertEqual(project["requires-python"], ">=3.11")
        self.assertEqual(project["license"], "MIT")
        self.assertEqual(
            project["dependencies"],
            ["numpy>=1.24", "PyYAML>=6", "scipy>=1.10"],
        )
        package_find = configuration["tool"]["setuptools"]["packages"]["find"]
        self.assertEqual(package_find["include"], ["pypospack*"])
        self.assertNotIn("tests*", package_find["include"])

    def test_parameter_modules_do_not_call_builtin_eval(self) -> None:
        violations = []
        for path in sorted(PACKAGE_ROOT.rglob("*.py")):
            module = ast.parse(path.read_text(encoding="utf-8"))
            for node in ast.walk(module):
                if (
                    isinstance(node, ast.Call)
                    and isinstance(node.func, ast.Name)
                    and node.func.id == "eval"
                ):
                    violations.append(f"{path.relative_to(PACKAGE_ROOT)}:{node.lineno}")
        self.assertEqual(violations, [])

    def test_yaml_loaders_are_safe(self) -> None:
        filesystem = ast.parse(
            (PACKAGE_ROOT / "io/filesystem.py").read_text(encoding="utf-8")
        )
        loader_bases = {
            ast.unparse(base)
            for node in filesystem.body
            if isinstance(node, ast.ClassDef) and node.name == "OrderedDictYAMLLoader"
            for base in node.bases
        }
        self.assertEqual(loader_bases, {"yaml.SafeLoader"})

        violations = []
        for path in sorted(PACKAGE_ROOT.rglob("*.py")):
            module = ast.parse(path.read_text(encoding="utf-8"))
            for node in ast.walk(module):
                if not (
                    isinstance(node, ast.Call)
                    and isinstance(node.func, ast.Attribute)
                    and isinstance(node.func.value, ast.Name)
                    and node.func.value.id == "yaml"
                    and node.func.attr == "load"
                ):
                    continue
                has_ordered_safe_loader = (
                    len(node.args) >= 2
                    and isinstance(node.args[1], ast.Name)
                    and node.args[1].id == "OrderedDictYAMLLoader"
                )
                if not has_ordered_safe_loader:
                    violations.append(f"{path.relative_to(PACKAGE_ROOT)}:{node.lineno}")
        self.assertEqual(violations, [])

    def test_offline_checksum_manifest_is_complete_and_paths_are_absent(
        self,
    ) -> None:
        lines = (
            (REPOSITORY_ROOT / "OFFLINE_ARTIFACT_SHA256SUMS")
            .read_text(encoding="utf-8")
            .splitlines()
        )
        paths = []
        for line in lines:
            match = _SHA256_LINE.fullmatch(line)
            self.assertIsNotNone(match)
            assert match is not None
            paths.append(match.group(1))
        self.assertEqual(len(paths), 5_396)
        self.assertEqual(paths, sorted(paths))
        self.assertEqual(len(paths), len(set(paths)))
        self.assertFalse(any((REPOSITORY_ROOT / path).exists() for path in paths))


if __name__ == "__main__":
    unittest.main()
