# Main workflow update

The main ASK and ASK-search code now includes the current breakpoint evidence,
circle selection, single-cell Q1 assignment, and ALL breakpoint-matrix workflow.

- `ask` uses the combined evidence workflow by default. Replace `--method both`
  with the default invocation; `--no-sa` disables SA candidates while retaining
  clip-site and PE/DS evidence. Recompute breakpoint pairs when changing this option.
- Use `-sc --fragments PATH` for single-cell assignment, optionally with
  `--sc-barcodes PATH`. Fragment generation is also available via
  `ask-make-fragments` or `--make-fragments`; see each command's `--help` for its
  external tool requirements and reuse rules.
- The final single-cell breakpoint support/normal matrices use raw ALL pairs
  passing the Count threshold (default 5). ALL/ecDNA annotation tables identify
  the column sets; binary cell/ecDNA assignment remains a separate Q1 result.
- The main workflow skips the unused legacy normal matrix. The standalone
  assignment and matrix APIs retain their default behavior; opting out with
  `compute_normal=False` returns `None` for the uncomputed legacy normal matrix.
- Step4 retains its four-item pickle structure and now stores circle scores in
  versioned metadata. Older caches without scores require normalization data to
  recompute nonempty circle scores. Failure to recover or write scores is explicit.
- The multi-Group breakpoint matcher now compares the two endpoint sequences.

Validation used synthetic indexed BAM fixtures and targeted cache/assignment
regressions. Production sample reruns and migration of historical cache versions
outside the tested paths have not been performed.
