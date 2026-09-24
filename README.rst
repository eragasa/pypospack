PyPosPack
=========

PyPosPack is a collection of research libraries for parameterizing parametric
interatomic potentials. Version 0.1.0 is an alpha source baseline: it defines a
bounded Python distribution and preservation boundary, but does not claim
scientific validation, calculator portability, workflow reproducibility, or
production readiness.

Installation
------------

The core distribution requires Python 3.11 or newer, NumPy, SciPy, and PyYAML::

    python -m pip install pypospack-0.1.0-py3-none-any.whl

Optional dependency groups are ``data``, ``structure``, ``plot``, and ``mpi``.
The ``all`` extra installs all of them. External calculators, scheduler clients,
and site-specific executables are intentionally not installed or discovered by
the package.

Release boundary
----------------

Only Python modules under ``pypospack`` are installed. Historical tests,
examples, development trees, potential and structure collections, generated
results, calculator output, private offline artifacts, and alternate source
trees are excluded from wheels and source distributions.

The release gate validates package syntax, safe YAML loading, bounded parameter
expression evaluation, archive boundaries, wheel installation, and selected
imports. It never launches LAMMPS, VASP, GULP, PHONTS, MPI jobs, schedulers, or
shell workflows.

Safety and support status
-------------------------

Calculator and workflow modules remain historical research interfaces. They may
write files or invoke external programs when explicitly called and are not part
of the validated v0.1.0 surface. Inputs from untrusted sources must not be
passed to those interfaces.

Historical unrestricted parameter-expression evaluation has been replaced by a
closed AST evaluator allowing only finite numeric constants, explicitly named
finite parameters, bounded arithmetic, comparisons, and boolean composition.
YAML readers use safe loaders, including the ordered mapping loader.

Generated and restricted artifacts removed from the maintained tip remain
recoverable from a private, checksum-bound offline archive. See
``OFFLINE_ARTIFACTS.md``. The private archive is not a package input or public
release asset.

Provenance and audit
--------------------

The release-readiness audit is documented in
``docs/release-readiness.md``. Existing Git history is retained and no
scientific acceptance claim is implied by a release tag.

Legacy environment variables
----------------------------

Historical workflows may refer to ``PYPOSPACK_ROOT_DIR`` and
``LAMMPS_SERIAL_BIN``. The package does not configure, validate, or consume
these values during installation or ordinary release tests.
