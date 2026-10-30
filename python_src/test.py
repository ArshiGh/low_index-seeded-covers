import unittest

from collections import Counter

from low_index import *

class TestSimsNode(unittest.TestCase):
    def test_doc(self):
        self.assertIn('A non-abstract SimsNode',
                      SimsNode.__doc__)
        self.assertIn('Create SimsNode for a covering graph ',
                      SimsNode.__init__.__doc__)
        self.assertIn('Copy a SimsNode',
                      SimsNode.__init__.__doc__)

    def test_add_edge(self):
        node = SimsNode(8, 4)
        self.assertEqual(node.rank, 8)

        node.add_edge(2, 1, 1)
        node.add_edge(3, 1, 2)
        node.add_edge(4, 2, 3)
        self.assertTrue(node.verified_add_edge(2, 1, 2))
        self.assertTrue(node.verified_add_edge(2, 2, 3))
        self.assertFalse(node.verified_add_edge(2, 3, 1))

        self.assertFalse(node.is_complete())
        self.assertEqual(node.degree, 3)
        self.assertIn('1--( 3)->2', str(node))

        with self.assertRaises(ValueError):
            node.permutation_rep()

    def test_permutation_rep(self):
        node = SimsNode(2, 2)
        node.add_edge(1, 1, 1)
        node.add_edge(2, 1, 2)
        node.add_edge(2, 2, 1)
        node.add_edge(1, 2, 2)
        self.assertTrue(node.is_complete())
        self.assertEqual(node.degree, 2)
        self.assertEqual(node.permutation_rep(), [[0, 1], [1, 0]])

class TestSimsTree(unittest.TestCase):
    def test_basic(self):
        t = SimsTree(2, 2, [], [])

        perm_reps = [ node.permutation_rep() for node in t.list() ]

        self.assertEqual(
            perm_reps,
            [[[0], [0]],
             [[0, 1], [1, 0]],
             [[1, 0], [0, 1]],
             [[1, 0], [1, 0]]])

    def test_figure_eight(self):
        #                             Manifold("m004").fundamental_group()
        #                     a  a  a  b  A  B  B  A  b
        t = SimsTree(2, 6, [[ 1, 1, 1, 2,-1,-2,-2,-1, 2],
                            [ 2, 1, 1, 1, 2,-1,-2,-2,-1],
                            [-1, 2, 1, 1, 1, 2,-1,-2,-2],
                            [-2,-1, 2, 1, 1, 1, 2,-1,-2],
                            [-2,-2,-1, 2, 1, 1, 1, 2,-1],
                            [-1,-2,-2,-1, 2, 1, 1, 1, 2],
                            [ 2,-1,-2,-2,-1, 2, 1, 1, 1],
                            [ 1, 2,-1,-2,-2,-1, 2, 1, 1],
                            [ 1, 1, 2,-1,-2,-2,-1, 2, 1],
                            [ 1, 1, 1, 2,-1,-2,-2,-1, 2]
                            ],
                     [])

        degrees = Counter([cover.degree for cover in t.list()])
        self.assertEqual(degrees[0], 0)
        self.assertEqual(degrees[1], 1)
        self.assertEqual(degrees[2], 1)
        self.assertEqual(degrees[3], 1)
        self.assertEqual(degrees[4], 2)
        self.assertEqual(degrees[5], 4)
        self.assertEqual(degrees[6], 11)
        self.assertEqual(degrees[7], 0) # Beyond what we counted

class TestPermutationRep(unittest.TestCase):
    def _test_K11n34_7(self, num_threads):
        reps = permutation_reps(
            3,
            ["aaBcbbcAc"],
            ["aacAbCBBaCAAbbcBc"],
            7,
            num_threads = num_threads)

        degrees = Counter([len(rep[0]) for rep in reps])
        self.assertEqual(degrees[1], 1)
        self.assertEqual(degrees[2], 1)
        self.assertEqual(degrees[3], 1)
        self.assertEqual(degrees[4], 1)
        self.assertEqual(degrees[5], 2)
        self.assertEqual(degrees[6], 16)
        self.assertEqual(degrees[7], 30)

        self.assertIn([[0], [0], [0]], reps)
        self.assertIn([[0, 3, 5, 4, 1, 2], [1, 0, 5, 2, 4, 3], [1, 4, 0, 3, 5, 2]], reps)

    def test_K11n34_7_single_threaded(self):
        self._test_K11n34_7(num_threads = 1)

    def test_K11n34_7_multi_threaded(self):
        self._test_K11n34_7(num_threads = 0)

    def test_K11n34_7_fixed_multi_threads(self):
        self._test_K11n34_7(num_threads = 48)

    def test_K15n12345_7(self):
        reps = permutation_reps(
            3,
            ["aBcACAcb"],
            ["aBaCacBAcAbaBabaCAcAbaBaCacBAcAbaBabCAcAbABaCabABAbABaCabCAcAb"],
            7)

        degrees = Counter([len(rep[0]) for rep in reps])
        self.assertEqual(degrees[1], 1)
        self.assertEqual(degrees[2], 1)
        self.assertEqual(degrees[3], 1)
        self.assertEqual(degrees[4], 1)
        self.assertEqual(degrees[5], 3)
        self.assertEqual(degrees[6], 11)
        self.assertEqual(degrees[7], 22)

        self.assertIn([[0], [0], [0]], reps)
        self.assertIn([[0, 2, 1, 3, 4, 6, 5], [0, 3, 4, 5, 1, 2, 6], [1, 2, 0, 5, 6, 3, 4]], reps)
        
    def test_o9_15405_9(self):
        reps = permutation_reps(
            2,
            [],
            ["aaaaabbbaabbbaaaaabbbaabbbaaaaaBBBBBBBB"],
            9)

        degrees = Counter([len(rep[0]) for rep in reps])
        self.assertEqual(degrees[1], 1)
        self.assertEqual(degrees[2], 1)
        self.assertEqual(degrees[3], 1)
        self.assertEqual(degrees[4], 1)
        self.assertEqual(degrees[5], 3)
        self.assertEqual(degrees[6], 3)
        self.assertEqual(degrees[7], 9)
        self.assertEqual(degrees[8], 5)
        self.assertEqual(degrees[9], 14)

        self.assertIn([[0], [0]], reps)
        self.assertIn([[0, 2, 4, 1, 5, 3], [1, 0, 5, 4, 2, 3]], reps)

    def _test_o9_03127_9(self, use_long):
        # Relator is longer than 255 and thus tests that we do not accidentally
        # use a 8-bit integer to index into the relator.
        
        relator = "aabbbaabbaabbbaabbbaaabbbaabbbaabbaabbbaabbbaabbaabbbaabbbaaabbbaabbbaabbaabbbaabbbaaabbbaabbbaabbaabbbaabbbaabbaabbbaabbbaaabbbaabbbaabbaabbbaabbbaaabbbaabbbaabbaabbbaabbbaabbaabbbaabbbaaabbbaabbbaabbaabbbaabbbaabbaabbbaabbbaaabbbaabbbaabbaabbbaabbbaaabbbaabbbaabbaabbb"
        
        reps = permutation_reps(
            2,
            [] if use_long else [ relator ],
            [ relator ] if use_long else [],
            4)

        degrees = Counter([len(rep[0]) for rep in reps])
        
        self.assertEqual(degrees[1], 1)
        self.assertEqual(degrees[2], 3)
        self.assertEqual(degrees[3], 2)
        self.assertEqual(degrees[4], 8)

        self.assertIn([[0, 1, 3, 2], [1, 3, 0, 2]], reps)
        
    def test_o9_03127_9_short(self):
        self._test_o9_03127_9(False)

    def test_o9_03127_9_long(self):
        self._test_o9_03127_9(True)

class TestSeededPermutationRep(unittest.TestCase):
    def test_seeded_basic(self):
        reps = permutation_reps_from_seed(
            1, [], [], 3, [[1, 0, 1]], num_threads=1)
        self.assertEqual(reps, [[[1, 0]], [[1, 2, 0]]])

    def test_seeded_preserves_seed(self):
        reps = permutation_reps_from_seed(
            2, [], [], 4, [[1, 0, 1]], num_threads=1)
        self.assertTrue(reps)
        self.assertTrue(all(rep[0][0] == 1 for rep in reps))

    def test_seeded_threading(self):
        args = dict(
            rank=2, short_relators=[], long_relators=[], max_degree=4,
            initial_edges=[[1, 0, 1]])
        single = permutation_reps_from_seed(**args, num_threads=1)
        threaded = permutation_reps_from_seed(**args, num_threads=4)
        self.assertEqual(single, threaded)

    def test_empty_seed_with_minimal_pruning_matches_unseeded(self):
        ordinary = permutation_reps(
            2, [], [], 4, num_threads=1)
        seeded = permutation_reps_from_seed(
            2, [], [], 4, [], num_threads=1, minimal_pruning=True)
        self.assertEqual(ordinary, seeded)

    def test_invalid_seed(self):
        self.assertEqual(
            permutation_reps_from_seed(
                2, [], [], 4, [[1, 0, 1], [1, 0, 0]], num_threads=1),
            [])
        self.assertEqual(
            permutation_reps_from_seed(
                2, [], [], 4, [[1, 2, 2]], num_threads=1),
            [])
        with self.assertRaises(ValueError):
            permutation_reps_from_seed(
                2, [], [], 4, [[1, 0]], num_threads=1)


    @staticmethod
    def _seed_symmetry_rep_key(rep):
        return tuple(tuple(generator) for generator in rep)

    @staticmethod
    def _seed_symmetry_conjugate_rep(rep, relabeling):
        inverse = [0] * len(relabeling)
        for old, new in enumerate(relabeling):
            inverse[new] = old

        return [
            [
                relabeling[generator[inverse[new]]]
                for new in range(len(generator))
            ]
            for generator in rep
        ]

    @classmethod
    def _three_cycle_seed_canonical(cls, rep):
        # Automorphisms of the prescribed 3-cycle seed.
        relabelings = (
            [0, 1, 2],
            [1, 2, 0],
            [2, 0, 1],
        )

        return min(
            cls._seed_symmetry_rep_key(
                cls._seed_symmetry_conjugate_rep(rep, relabeling))
            for relabeling in relabelings
        )

    @staticmethod
    def _three_cycle_seed_args():
        return dict(
            rank=2,
            short_relators=[],
            long_relators=[],
            max_degree=3,
            initial_edges=[
                [1, 0, 1],
                [1, 1, 2],
                [1, 2, 0],
            ],
            minimal_pruning=False,
        )

    def test_seed_symmetry_default_is_off(self):
        args = self._three_cycle_seed_args()

        omitted = permutation_reps_from_seed(
            **args, num_threads=1)

        explicit_off = permutation_reps_from_seed(
            **args,
            num_threads=1,
            seed_symmetry_pruning=False)

        self.assertEqual(omitted, explicit_off)

    def test_seed_symmetry_prunes_seed_automorphism_duplicates(self):
        args = self._three_cycle_seed_args()

        off = permutation_reps_from_seed(
            **args,
            num_threads=1,
            seed_symmetry_pruning=False)

        on = permutation_reps_from_seed(
            **args,
            num_threads=1,
            seed_symmetry_pruning=True)

        # The optimization must actually remove labeled duplicates.
        self.assertLess(len(on), len(off))

        # But OFF and ON must represent exactly the same quotient by
        # automorphisms of the fixed seed.
        off_classes = {
            self._three_cycle_seed_canonical(rep)
            for rep in off
        }
        on_classes = {
            self._three_cycle_seed_canonical(rep)
            for rep in on
        }

        self.assertEqual(off_classes, on_classes)
        self.assertEqual(len(on_classes), len(on))

    def test_seed_symmetry_threading(self):
        args = self._three_cycle_seed_args()

        single = permutation_reps_from_seed(
            **args,
            num_threads=1,
            seed_symmetry_pruning=True)

        threaded = permutation_reps_from_seed(
            **args,
            num_threads=4,
            seed_symmetry_pruning=True)

        self.assertEqual(
            {
                self._seed_symmetry_rep_key(rep)
                for rep in single
            },
            {
                self._seed_symmetry_rep_key(rep)
                for rep in threaded
            })

    def test_seed_symmetry_disconnected_seed_falls_back(self):
        # Degree is already 3, but the prescribed seed graph has
        # components {0,1} and {2}.  Seed-symmetry pruning is therefore
        # intentionally inapplicable and must leave behavior unchanged.
        args = dict(
            rank=2,
            short_relators=[],
            long_relators=[],
            max_degree=3,
            initial_edges=[
                [1, 0, 1],
                [1, 1, 0],
                [1, 2, 2],
            ],
            num_threads=1,
            minimal_pruning=False,
        )

        off = permutation_reps_from_seed(
            **args,
            seed_symmetry_pruning=False)

        on = permutation_reps_from_seed(
            **args,
            seed_symmetry_pruning=True)

        self.assertEqual(off, on)

    def test_seed_symmetry_invalid_seed_stays_invalid(self):
        self.assertEqual(
            permutation_reps_from_seed(
                2, [], [], 4,
                [[1, 0, 1], [1, 0, 0]],
                num_threads=1,
                seed_symmetry_pruning=True),
            [])

    def test_seed_negative_letter_rejected(self):
        with self.assertRaises(TypeError):
            permutation_reps_from_seed(
                1, [], [], 2,
                [[-1, 0, 1]],
                num_threads=1,
                seed_symmetry_pruning=False)

    def test_seed_symmetry_public_api_documented(self):
        self.assertIn(
            "seed_symmetry_pruning",
            permutation_reps_from_seed.__doc__)


class TestFullCycleConstraints(unittest.TestCase):

    @staticmethod
    def _rank1_args(max_degree=4):
        return dict(
            rank=1,
            short_relators=[],
            long_relators=[],
            max_degree=max_degree,
            initial_edges=[[1, 0, 1]],
            minimal_pruning=False,
        )

    @staticmethod
    def _rep_key(rep):
        return tuple(
            tuple(generator)
            for generator in rep
        )

    def test_full_cycle_default_empty_matches_omitted(self):
        args = self._rank1_args(max_degree=3)

        omitted = permutation_reps_from_seed(
            **args,
            num_threads=1)

        explicit_empty = permutation_reps_from_seed(
            **args,
            num_threads=1,
            full_cycle_constraints=[])

        self.assertEqual(
            omitted,
            explicit_empty)

    def test_full_cycle_order2(self):
        reps = permutation_reps_from_seed(
            **self._rank1_args(max_degree=3),
            num_threads=1,
            full_cycle_constraints=[
                ([1], 2),
            ])

        self.assertEqual(
            reps,
            [[[1, 0]]])

    def test_full_cycle_order3(self):
        reps = permutation_reps_from_seed(
            **self._rank1_args(max_degree=3),
            num_threads=1,
            full_cycle_constraints=[
                ([1], 3),
            ])

        self.assertEqual(
            reps,
            [[[1, 2, 0]]])

    def test_full_cycle_order4_open_path_survives(self):
        # The seed specifies only 0 -> 1, so at the root the
        # four-cycle is still open.  The search must not prune it.
        reps = permutation_reps_from_seed(
            **self._rank1_args(max_degree=4),
            num_threads=1,
            full_cycle_constraints=[
                ([1], 4),
            ])

        self.assertTrue(reps)

        for rep in reps:
            permutation = rep[0]

            self.assertEqual(
                len(permutation),
                4)

            v = 0
            seen = []

            for _ in range(4):
                seen.append(v)
                v = permutation[v]

            self.assertEqual(v, 0)
            self.assertEqual(
                len(set(seen)),
                4)

    def test_full_cycle_root_fixed_point_rejected(self):
        reps = permutation_reps_from_seed(
            rank=1,
            short_relators=[],
            long_relators=[],
            max_degree=1,
            initial_edges=[
                [1, 0, 0],
            ],
            num_threads=1,
            minimal_pruning=False,
            full_cycle_constraints=[
                ([1], 2),
            ])

        self.assertEqual(reps, [])

    def test_full_cycle_shorter_closed_cycle_rejected(self):
        # A completely specified two-cycle cannot later become
        # a four-cycle.
        reps = permutation_reps_from_seed(
            rank=1,
            short_relators=[],
            long_relators=[],
            max_degree=2,
            initial_edges=[
                [1, 0, 1],
                [1, 1, 0],
            ],
            num_threads=1,
            minimal_pruning=False,
            full_cycle_constraints=[
                ([1], 4),
            ])

        self.assertEqual(reps, [])

    def test_full_cycle_signed_word(self):
        # x = (0 1), y = identity.
        # Therefore y^-1 x = x, which is a full two-cycle.
        accepted = permutation_reps_from_seed(
            rank=2,
            short_relators=[],
            long_relators=[],
            max_degree=2,
            initial_edges=[
                [1, 0, 1],
                [1, 1, 0],
                [2, 0, 0],
                [2, 1, 1],
            ],
            num_threads=1,
            minimal_pruning=False,
            full_cycle_constraints=[
                ([-2, 1], 2),
            ])

        self.assertTrue(accepted)

        # x = y = (0 1), hence y^-1 x is the identity.
        # It must fail the order-two full-cycle condition.
        rejected = permutation_reps_from_seed(
            rank=2,
            short_relators=[],
            long_relators=[],
            max_degree=2,
            initial_edges=[
                [1, 0, 1],
                [1, 1, 0],
                [2, 0, 1],
                [2, 1, 0],
            ],
            num_threads=1,
            minimal_pruning=False,
            full_cycle_constraints=[
                ([-2, 1], 2),
            ])

        self.assertEqual(rejected, [])

    def test_full_cycle_single_and_multithread_agree(self):
        args = self._rank1_args(max_degree=4)

        single = permutation_reps_from_seed(
            **args,
            num_threads=1,
            full_cycle_constraints=[
                ([1], 4),
            ])

        threaded = permutation_reps_from_seed(
            **args,
            num_threads=4,
            full_cycle_constraints=[
                ([1], 4),
            ])

        self.assertEqual(
            {
                self._rep_key(rep)
                for rep in single
            },
            {
                self._rep_key(rep)
                for rep in threaded
            })

    def test_full_cycle_invalid_constraints(self):
        args = self._rank1_args(max_degree=3)

        invalid = [
            [([1], 0)],
            [([], 2)],
            [([2], 2)],
        ]

        for constraints in invalid:
            with self.assertRaises(ValueError):
                permutation_reps_from_seed(
                    **args,
                    num_threads=1,
                    full_cycle_constraints=constraints)

    def test_full_cycle_public_api_documented(self):
        self.assertIn(
            "full_cycle_constraints",
            permutation_reps_from_seed.__doc__)

if __name__ == '__main__':
    print("Number of cores reported by the operating system:",
          hardware_concurrency())
    unittest.main()


