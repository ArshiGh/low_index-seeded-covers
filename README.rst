Low Index Subgroups
===================

The low_index project provides a Python module which implements a variant
of Charles Sims' Low Index Subgroups algorithm for enumerating all of
the conjugacy classes of subgroups of a finitely presented group with
finite index less than a given bound.

The upstream ``low_index`` package is available on PyPI.  This repository is
a downstream extension of low_index 1.3 for seeded and constrained cover
enumeration, released as ``1.3+seeded.covers.2``.

For this release, install a wheel from the corresponding GitHub release or
install directly from a source checkout:

.. code-block:: bash

    python -m pip install .

The ordinary ``permutation_reps`` API remains available.  Its output
compatibility is regression-tested against low_index 1.3.

Here is a sample computation:

.. code-block:: python

    >>> from low_index import *
    >>>  # Conjugacy classes of subgroups of F_3 with index at most 4: 
    >>> reps = permutation_reps(rank = 3, short_relators = [], long_relators = [], max_degree = 4)
    >>> len(reps)
    653
    >>> from snappy import *
    >>> G = Manifold('K11n34').fundamental_group(); G
    Generators:
       a,b,c
    Relators:
       aaBcbbcAc
       aacAbCBBaCAAbbcBc
    >>> # Degree at most 7 covers of the exterior of the Conway knot:
    >>> reps = permutation_reps(G.num_generators(), G.relators()[:1], G.relators()[1:], 7)
    >>> len(reps)
    52
    >>> reps[25]
    [[1, 0, 3, 2, 5, 6, 4], [1, 4, 0, 6, 2, 5, 3], [3, 0, 2, 6, 4, 1, 5]]


Seeded cover enumeration extension
==================================

This release extends low_index 1.3 with
``permutation_reps_from_seed`` for searches in which part of the desired
covering action is already prescribed.  The ordinary ``permutation_reps``
entry point is retained.

A seed is supplied through ``initial_edges`` as triples

``[generator, from_vertex, to_vertex]``.

Generators in ``initial_edges`` are positive and numbered from 1.  Vertices
are numbered from 0, matching the returned permutation representations.  The
seed must describe a connected partial covering containing vertex 0.

Seeded searches often use a labeling chosen to express the desired geometric
or combinatorial structure rather than the canonical labeling used by the
ordinary Sims search.  Therefore ``minimal_pruning`` defaults to ``False`` for
``permutation_reps_from_seed``.

Seed symmetry pruning
---------------------

Setting

``seed_symmetry_pruning=True``

allows a connected exact-degree seeded search to quotient by automorphisms
that preserve the prescribed seed.  This is independent of the ordinary
canonical minimality test.  The option defaults to ``False``.

Exact full-cycle constraints
----------------------------

``full_cycle_constraints`` accepts pairs

``(word, order)``

where ``word`` is a signed generator word and ``order`` is a positive integer.

For a complete permutation representation, a constraint requires the
permutation induced by ``word`` to consist entirely of cycles of exactly the
specified length.  During a partial search, an open path is still allowed
because later edges may complete it.  A path that closes before the required
order is rejected immediately.

For example::

    from low_index import permutation_reps_from_seed

    # seed_edges describes a prescribed connected degree-24 seed.
    reps = permutation_reps_from_seed(
        rank=4,
        short_relators=relators,
        long_relators=[],
        max_degree=24,
        initial_edges=seed_edges,
        num_threads=8,
        minimal_pruning=False,
        seed_symmetry_pruning=True,
        full_cycle_constraints=[
            ([1], 3),
            ([-2, 1], 2),
        ],
    )

Performance implementation
--------------------------

The seeded search reuses the low_index Sims enumeration engine and includes
the performance work developed for constrained cover enumeration: endpoint
worklists for relator propagation, a rollback-maintained relator wake index,
and elimination of full-cycle constraints that are already completely
determined by an exact-degree root seed.

Testing
-------

The release contains 34 regression tests covering the original API, seeded
completion, seed validation, single- and multi-threaded agreement,
seed-symmetry pruning, and exact full-cycle constraints.

After installation, run::

    python -m low_index.test

Additional implementation details are recorded in
``SEEDED_MODIFICATIONS.md`` and release-specific notes in
``RELEASE_NOTES.md``.


Credits
=======

Primarily developed by `Marc Culler <https://marc-culler.info>`_,
`Nathan Dunfield <http://dunfield.info>`_, and `Matthias Goerner
<http://www.unhyperbolic.org/>`_

License
=======

Copyright 2022 by Marc Culler, Nathan Dunfield, Matthias Goerner
and others.

This code is released under the `GNU General Public License, version 2
<http://www.gnu.org/licenses/gpl-2.0.txt>`_ or (at your option) any
later version as published by the Free Software Foundation.

