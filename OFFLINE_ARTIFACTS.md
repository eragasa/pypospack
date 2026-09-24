# Offline artifact policy

Generated calculator results, restricted VASP inputs, scheduler output, caches,
and generated documentation are not maintained release-source inputs. The
source-specific selection policy is `devtools/offline-artifacts.toml`; generic
staging, verification, archive, and recovery mechanics are owned by
`projectkoios-bootstrap`.

From original commit `21cdecaf3b05c87acc532d992be2c04d85bfbc22`
(tree `a5a4cee972e487512275c34f308251e6cc38c2fa`), the reviewed policy selected
5,396 files totaling 637,072,563 bytes. They were copied into ignored local
staging, byte-verified against that Git object, and removed from the maintained
tip without rewriting history.

The private preservation bundle is stored locally under:

```text
~/Library/CloudStorage/Dropbox/pypospack/v0.1/data/
```

Its archive is:

```text
pypospack-offline-artifacts-git-21cdecaf3b05.tar
SHA-256 386e854f42616966ae31fe7d266077bed71d58a911c5395d62a4e625e60acfbf
```

A complete recovery into a new directory verified all 5,396 paths, sizes, and
checksums against both the canonical manifest and the original Git commit. The
bundle contains its manifest, source identity, checksum evidence, and recovery
instructions. Dropbox remote synchronization and an independent second backup
must be confirmed separately.

`OFFLINE_ARTIFACT_SHA256SUMS` records conventional per-file checksums without
including private bytes. Historical Git objects and the private archive remain
recovery evidence only. Package code and tests must not discover, restore,
import, or execute them. In particular, no redistribution authority is claimed
for removed pseudopotentials.
