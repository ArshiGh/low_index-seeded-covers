# Release 1.3+seeded.covers.1

This is a downstream seeded-cover extension of `low_index` 1.3.

It is not an upstream `low_index` 1.3.1 release.

## Added

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

## Compatibility

The ordinary `permutation_reps` API remains available and follows the
low_index 1.3 ordinary search path.

## Validation

The release candidate contains 34 regression tests, including tests for seeded
enumeration, symmetry pruning, exact full-cycle constraints, invalid inputs,
and threaded equivalence.
