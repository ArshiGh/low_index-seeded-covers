# low_index 1.3 seeded-cover extension

This source tree is a downstream extension of released `low_index` 1.3 for
seeded and constrained cover enumeration.

The release version is:

`1.3+seeded.covers.1`

It is not intended to represent an upstream `low_index` 1.3.1 release.

## Public API

The ordinary `permutation_reps` API is retained.

The extension adds:

`permutation_reps_from_seed(...)`

with the additional controls:

- `initial_edges`
- `minimal_pruning`
- `seed_symmetry_pruning`
- `full_cycle_constraints`

### Initial edges

`initial_edges` is a list of triples

`[generator, from_vertex, to_vertex]`.

Generators are positive and numbered from 1. Vertices are numbered from 0.

The seed must describe a valid connected partial covering containing vertex 0.
Conflicting edges, invalid generator numbers, invalid vertex indices, and
incompatible short relators are rejected before recursive enumeration.

### Minimal pruning

Seeded searches default to:

`minimal_pruning=False`

because a user-specified seed may use a structural labeling that is not
compatible with the canonical labeling assumed by the ordinary conjugacy-class
minimality test.

### Seed-symmetry pruning

`seed_symmetry_pruning=True` quotients a connected exact-degree seeded search
by automorphisms preserving the prescribed seed.

It is independent of ordinary canonical minimality pruning and defaults to
`False`.

### Full-cycle constraints

`full_cycle_constraints` is a list of pairs

`(word, order)`.

For a complete representation, each constraint requires the permutation
induced by `word` to consist entirely of cycles of exactly `order`.

During recursive enumeration:

- an open path is still potentially valid;
- a cycle closing before the required order is rejected;
- if all `order` repetitions are known and do not return to the starting
  sheet, the branch is rejected.

Signed generator words are supported.

## Performance implementation

The final seeded-cover implementation includes the following accepted
performance improvements.

### Relator endpoint worklist

Relator propagation revisits affected lift states through an endpoint-oriented
worklist rather than repeatedly scanning all states.

### Rollback-maintained relator wake index

A reverse wake index is maintained incrementally during depth-first search and
rolled back with the search state. This avoids rebuilding the reverse index for
every recursive call.

### Root-fixed full-cycle elimination

For an exact-degree seed root, full-cycle constraints that are already fully
determined and valid at the root are removed from descendant checks.

This preserves the constraint semantics while avoiding repeated evaluation of
conditions that cannot change below that root.

## Compatibility

Calling the ordinary `permutation_reps` entry point continues to use the
ordinary low_index search path.

An empty seeded search with canonical minimality enabled is regression-tested
against the ordinary API.

Both single-threaded and multi-threaded seeded searches are regression-tested.

## Validation

The release contains 34 regression tests covering:

- ordinary low-index enumeration;
- seeded completion;
- seed preservation and validation;
- empty-seed compatibility;
- single- versus multi-threaded enumeration;
- seed-symmetry pruning;
- full-cycle constraints of several orders;
- signed full-cycle words;
- open-path behavior;
- premature-cycle rejection;
- invalid full-cycle constraints;
- public API exposure.

The constrained-cover workflow was additionally validated against frozen exact
output references during development before this release candidate was
prepared.
