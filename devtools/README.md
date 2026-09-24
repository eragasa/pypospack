# Release-maintenance tools

`offline-artifacts.toml` is the reviewed PyPosPack selection policy for generated
and restricted files that must leave the maintained release tip. It contains no
copying, removal, archive, or recovery implementation.

The generic implementation is pinned to `projectkoios-bootstrap` commit
`c0f406685a49838a690c1519a6174564bbaa057b` (tree
`3be0101c8a582e9cdcaa943197176770d52c7f78`). From a clean PyPosPack checkout,
apply it explicitly:

```bash
PYTHONPATH=/explicit/projectkoios-bootstrap/python python3.14 -m \
  projectkoios.bootstrap.harness.offline_artifacts stage \
  --repository-root . \
  --policy devtools/offline-artifacts.toml \
  --destination .offline/release-readiness \
  --checksum-manifest OFFLINE_ARTIFACT_SHA256SUMS \
  --remove-originals
```

The destination is ignored and must not exist. The command refuses a dirty
repository, copies each selected tracked regular file into staging, computes the
canonical manifest and checksums, publishes the complete staging tree, and only
then removes originals. Source-specific selection remains reviewable here while
the safety-sensitive mechanics stay in their owning repository.

Do not execute any staged Python, calculator, scheduler, shell, or expression
content. The staged bytes are private preservation evidence and never package
inputs.
