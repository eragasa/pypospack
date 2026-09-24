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

Disposition: `MUST_FIX`

The repository does not define a distribution, package selection, version,
supported Python versions, or dependencies. The smallest correction is modern
metadata that explicitly selects only `pypospack` packages. Development,
examples, historical tests, generated documentation, local data, offline
artifacts, and calculator output must remain outside wheel and source-
distribution boundaries.

### Maintained package contains incomplete syntax

Disposition: `MUST_FIX`

Static parsing fails in:

- `pypospack/optimization/loss_function.py:39`; and
- `pypospack/pyposmat/engines/base_sampler.py:61`.

The package modules must be repaired or excluded explicitly. The other 95
failures are under development, examples, data, tests, historical tests, and
an alternate `src/version_1` tree. Those trees must not silently enter release
packages.

### Generated and restricted artifacts are tracked

Disposition: `FIXED`

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

Disposition: `MUST_FIX`

The package has numerical, plotting, analysis, MPI, structure-tooling, and YAML
dependencies but no metadata. It also contains calculator and shell process
boundaries and historical expression evaluation. Dependencies must be divided
into core and optional surfaces. Installation and ordinary CI must never run a
calculator, scheduler, shell command, retained expression, or workflow.

This source baseline may preserve effectful APIs only if they are clearly
classified as explicit, unverified research interfaces and are never exercised
by release validation. Revisit them before presenting any calculator-backed
feature as supported.

### Portable test and CI contract is absent

Disposition: `MUST_FIX`

The historical test and example trees mix useful source, machine-specific
paths, generated outputs, incomplete scripts, and calculator execution. A
release gate must instead begin with sanitized static and in-memory tests for
package syntax, metadata, selection policy, and distribution boundaries. CI
must import maintained modules only after installing declared dependencies and
must not discover the historical test tree automatically.

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
5. Add explicit package and source-distribution boundaries.
6. Repair or exclude the two invalid package modules.
7. Declare core and optional dependencies from the maintained import surface.
8. Add bounded tests and CI without calculator or scheduler execution.
9. Build, inspect, install, and import the distributions in isolation.
10. Reassess findings before merging and creating an annotated release tag.

## Review result

Review outcome: `CHANGES_REQUIRED`

Artifact externalization, packaging, syntax, dependency, test, CI, and support-
boundary work remains before a release is appropriate.
