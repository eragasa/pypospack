from __future__ import annotations

import sys
import tempfile
import unittest
import zipfile
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]

sys.path.insert(0, str(REPOSITORY_ROOT / "tools"))
from verify_wheel import WheelBoundaryError, verify_wheel  # noqa: E402


class WheelBoundaryTest(unittest.TestCase):
    def test_accepts_recursive_package_and_distribution_metadata(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            wheel = Path(directory) / "pypospack-0.1.0-py3-none-any.whl"
            self._write_wheel(wheel)

            names = verify_wheel(REPOSITORY_ROOT, wheel)

            self.assertIn("pypospack/__init__.py", names)
            self.assertIn("pypospack/pyposmat/data/datafile.py", names)
            self.assertFalse(any(name.startswith("examples/") for name in names))

    def test_rejects_file_outside_package_boundary(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            wheel = Path(directory) / "pypospack-0.1.0-py3-none-any.whl"
            self._write_wheel(wheel, extra_name="examples/output.dat")

            with self.assertRaisesRegex(
                WheelBoundaryError,
                "outside its boundary",
            ):
                verify_wheel(REPOSITORY_ROOT, wheel)

    def test_rejects_non_python_package_data(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            wheel = Path(directory) / "pypospack-0.1.0-py3-none-any.whl"
            self._write_wheel(wheel, extra_name="pypospack/output.dat")

            with self.assertRaisesRegex(
                WheelBoundaryError,
                "non-source package data",
            ):
                verify_wheel(REPOSITORY_ROOT, wheel)

    @staticmethod
    def _write_wheel(wheel: Path, *, extra_name: str | None = None) -> None:
        with zipfile.ZipFile(wheel, mode="w") as archive:
            for source in sorted((REPOSITORY_ROOT / "pypospack").rglob("*.py")):
                archive.writestr(
                    source.relative_to(REPOSITORY_ROOT).as_posix(),
                    source.read_bytes(),
                )
            archive.writestr("pypospack-0.1.0.dist-info/METADATA", "metadata")
            if extra_name is not None:
                archive.writestr(extra_name, "unexpected")


if __name__ == "__main__":
    unittest.main()
