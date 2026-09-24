# Release-readiness audit

This audit begins from PyPosPack commit
`21cdecaf3b05c87acc532d992be2c04d85bfbc22` and Git tree
`a5a4cee972e487512275c34f308251e6cc38c2fa` on
`dev/release-readiness`. Inspection is static except for an isolated package
build attempt. No PyPosPack module, calculator, scheduler, retained script, or
configuration expression was executed.

## Definition of ready

A source baseline is release-ready when it has an explicit package boundary,
reproducible metadata and dependency declarations, buildable distributions, a
bounded non-calculator test suite, CI enforcement, documented limitations, and
no restricted or generated calculator artifacts in the maintained tip or
release packages. An annotated version tag remains the canonical release
identity.

## Evidence summary

- Tracked inventory: 11,516 files and 810,361,268 bytes.
- Largest roots: `tests/` 343,165,225 bytes, `dev/` 314,333,029 bytes,
  `examples/` 123,714,633 bytes, and `docs/` 12,725,310 bytes.
- Maintained package candidate: 149 files and 1,318,571 bytes under
  `pypospack/`.
- Python inventory: 1,275 files; 97 fail static parsing, including two beneath
  `pypospack/`.
- CI configuration: absent.
- `requirements.txt` and `setup.py`: both empty.
- Isolated wheel build: failed because the empty `setup.py` produced no
  distribution.
- External package imports include NumPy, SciPy, pandas, matplotlib,
  scikit-learn, PyYAML, ASE, mpi4py, and Bokeh.
- Effectful package paths include shell-enabled subprocess calls, calculator
  launch boundaries, unrestricted expression evaluation, and historical YAML
  loading APIs.

## Findings

### Package boundary is absent

Disposition: `NO_ACTION_REQUIRED`

`pyproject.toml` now defines version 0.1.0, Python 3.11 or newer, MIT licensing,
core and optional dependencies, and recursive selection limited to
`pypospack`. Wheel and source-distribution verifiers reject content outside the
maintained boundary. Development, examples, historical tests, local data,
offline artifacts, alternate source trees, and calculator output are absent
from both distributions.

### Maintained package contains incomplete syntax

Disposition: `NO_ACTION_REQUIRED`

The two incomplete package modules and their unusable dependent sampler modules
were removed from the maintained tree rather than assigned unverified
semantics. Historical Git retains them. The remaining 142 package modules parse
under Python 3.14 with `SyntaxWarning` promoted to an error. The other 95
historical failures remain outside release distributions.

### Generated and restricted artifacts are tracked

Disposition: `NO_ACTION_REQUIRED`

The reviewed declarative policy selected 5,396 files totaling 637,072,563
bytes. The selection includes VASP pseudopotentials and outputs,
LAMMPS outputs and restarts, PHONTS result datasets and archives, optimization
logs and result tables, scheduler output, generated documentation, editor
state, platform metadata, and interpreter caches.

The policy is tracked at `devtools/offline-artifacts.toml`. Selection is
repository-owned; generic staging, verification, archive, and recovery are
provided by `projectkoios-bootstrap` commit
`c0f406685a49838a690c1519a6174564bbaa057b`. The selected bytes were copied and
verified against the original Git commit, archived privately, recovered and
reverified in full, and removed from the maintained tip without rewriting
history. `OFFLINE_ARTIFACTS.md` records the bounded result and limitations.

### Runtime and effect boundaries are undeclared

Disposition: `NO_ACTION_REQUIRED`

Core numerical and YAML dependencies and optional data, structure, plotting,
and MPI groups are explicit. Unrestricted parameter-expression evaluation was
replaced with bounded AST arithmetic over explicitly supplied names; package
YAML readers now use safe loaders. Calculator, scheduler, shell, and workflow
interfaces are documented as explicit, unverified historical research
surfaces. Installation and release validation do not invoke them. Revisit these
interfaces before presenting calculator-backed behavior as supported.

### Portable test and CI contract is absent

Disposition: `NO_ACTION_REQUIRED`

`release_tests/` provides sanitized static and in-memory checks for package
syntax, metadata, offline evidence, safe loading, bounded expressions, and wheel
contents. CI on Python 3.11 and 3.14 builds both distributions, validates their
boundaries, installs the wheel with bounded optional dependencies, runs
`pip check`, and imports selected maintained modules outside the checkout. It
does not discover the historical test tree or execute calculator workflows.

### Potential and structure collections need a separate ownership decision

Disposition: `SAFE_TO_DEFER`

Potential definitions, EAM tables, and structure inputs remain outside the
initial generated-artifact selection because filename extensions alone do not
establish authorship, redistribution rights, or generated status. They can
remain in the Git source baseline while all are excluded from release packages.

Revisit this classification before distributing any collection as package
data, copying it to a new repository, or making a compatibility or scientific
claim about it.

## Recommended sequence

1. ~~Commit the declarative selection policy and audit evidence.~~
2. ~~Apply the generic bootstrap staging command from a clean checkout with
   explicit original removal.~~
3. ~~Verify selected bytes against the original Git commit and create a private
   Dropbox archive under `pypospack/v0.1/data/`.~~
4. ~~Perform and document a complete recovery round trip.~~
5. ~~Add explicit package and source-distribution boundaries.~~
6. ~~Repair or exclude the two invalid package modules.~~
7. ~~Declare core and optional dependencies from the maintained import
   surface.~~
8. ~~Add bounded tests and CI without calculator or scheduler execution.~~
9. ~~Build, inspect, install, and import the distributions in isolation.~~
10. Reassess findings after remote CI before merging and creating an annotated
    release tag.

## Review result

Review outcome: `NO_BLOCKING_FINDINGS`

Operator request: `NONE`

Local release validation found no current technical blocker. This outcome does
not imply scientific validation, calculator portability, publication, or
release authorization. Remote CI and exact release-object evidence remain the
next checks before tagging.
