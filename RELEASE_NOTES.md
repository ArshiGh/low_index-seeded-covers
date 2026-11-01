# Release 1.3+seeded.covers.2

This is a downstream seeded-cover extension of `low_index` 1.3.

It is not an upstream `low_index` 1.3.2 release.

## Added in the downstream seeded-cover series

- `permutation_reps_from_seed`
- prescribed partial coset-table seeds through `initial_edges`
- optional seeded canonical minimality control
- seed-preserving symmetry pruning
- generic exact full-cycle constraints
- single- and multi-threaded support for the seeded API

## Performance work included

- relator endpoint worklists
- rollback-maintained relator wake indexing
- root-fixed elimination of already-determined full-cycle constraints

## Corrected in 1.3+seeded.covers.2

The seeded-root initializer now repeats broad short-relator propagation until
a complete pass adds no new edge before the rollback wake index is built.
This establishes a synchronized initial wake index for arbitrary accepted
seeded roots, rather than relying on a single initialization pass.

## Compatibility

The ordinary `permutation_reps` API remains available.  Its output
compatibility is regression-tested against low_index 1.3.

## Validation

The package contains 34 regression tests covering seeded enumeration, seed
validation, symmetry pruning, exact full-cycle constraints, invalid inputs,
ordinary-API compatibility, and threaded equivalence.

The root-synchronization correction was additionally validated on MSI against
the frozen degree-24 constrained outputs for O333_2, O333_3, and O333_4.  The
patched build returned 46, 46, and 51 representations respectively, with exact
byte-for-byte equality to the frozen reference files.  The validation Slurm
job completed successfully with exit code 0:0.
