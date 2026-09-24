from __future__ import annotations

import argparse
import tarfile
from pathlib import Path, PurePosixPath

_FORBIDDEN_ROOTS = frozenset(
    {
        ".offline",
        "data",
        "dev",
        "examples",
        "src",
        "tests",
        "tests_old",
        "visualization_apps",
    }
)
_REQUIRED_PATHS = frozenset(
    {
        "LICENSE",
        "MANIFEST.in",
        "OFFLINE_ARTIFACTS.md",
        "OFFLINE_ARTIFACT_SHA256SUMS",
        "README.rst",
        "docs/release-readiness.md",
        "pyproject.toml",
        "release_tests/test_expressions.py",
        "release_tests/test_package_boundary.py",
        "release_tests/test_wheel_boundary.py",
        "requirements.txt",
        "setup.py",
        "tools/verify_sdist.py",
        "tools/verify_wheel.py",
    }
)


class SourceDistributionBoundaryError(ValueError):
    """Raised when a source distribution violates its release boundary."""


def verify_sdist(archive_path: Path) -> tuple[str, ...]:
    """Verify source-distribution paths without extracting the archive."""
    source = archive_path.resolve(strict=True)
    with tarfile.open(source, mode="r:gz") as archive:
        members = archive.getmembers()
    names = [member.name for member in members]
    if len(names) != len(set(names)):
        raise SourceDistributionBoundaryError(
            "source distribution contains duplicate paths"
        )
    roots = {PurePosixPath(name).parts[0] for name in names if name}
    if len(roots) != 1:
        raise SourceDistributionBoundaryError(
            "source distribution must have one root directory"
        )
    root = next(iter(roots))
    relative_files: set[str] = set()
    for member in members:
        path = PurePosixPath(member.name)
        if (
            path.is_absolute()
            or "." in path.parts
            or ".." in path.parts
            or not path.parts
            or path.parts[0] != root
        ):
            raise SourceDistributionBoundaryError(
                f"source distribution contains an unsafe path: {member.name}"
            )
        if member.issym() or member.islnk():
            raise SourceDistributionBoundaryError(
                f"source distribution contains a link: {member.name}"
            )
        if len(path.parts) > 1 and path.parts[1] in _FORBIDDEN_ROOTS:
            raise SourceDistributionBoundaryError(
                f"source distribution contains excluded content: {member.name}"
            )
        if member.isfile():
            relative_files.add(PurePosixPath(*path.parts[1:]).as_posix())
    missing = sorted(_REQUIRED_PATHS - relative_files)
    if missing:
        raise SourceDistributionBoundaryError(
            f"source distribution omits required files: {missing}"
        )
    package_files = {path for path in relative_files if path.startswith("pypospack/")}
    if any(not path.endswith(".py") for path in package_files):
        raise SourceDistributionBoundaryError(
            "source distribution contains non-source package data"
        )
    if any(path.endswith((".pyc", ".pyo")) for path in relative_files):
        raise SourceDistributionBoundaryError(
            "source distribution contains interpreter output"
        )
    return tuple(sorted(relative_files))


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(
        description="Verify the PyPosPack source-distribution boundary."
    )
    result.add_argument("archive", type=Path)
    return result


def main(arguments: list[str] | None = None) -> int:
    args = parser().parse_args(arguments)
    for name in verify_sdist(args.archive):
        print(name)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
