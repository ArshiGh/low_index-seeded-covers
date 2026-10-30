#ifndef LOW_INDEX_SIMS_TREE_BASE_H
#define LOW_INDEX_SIMS_TREE_BASE_H

#include "fullCycleConstraints.h"
#include "simsNode.h"
#include <algorithm>
#include <vector>

namespace low_index {

/// A base class for algorithms to find all covering graphs of a given
/// finitely presented group G up to a given degree.
///
/// The constructor of an implementation of this class will at least take
/// the number of generators of G (rank), the max_degree and the relators of
/// the group partitioned into short_relators and long_relators.
///
/// As explained for permutation_reps, short_relators are
/// checked when recursing the incomplete covering subgraphs and long_relators
/// are only checked at the end.
///
class SimsTreeBase
{
public:
    /// Find all complete covering subgraphs for the given group G.
    ///
    /// Call this only once - some implementations store state so
    /// the second call will give the wrong result.
    ///
    std::vector<SimsNode> list();

    virtual ~SimsTreeBase();
    
protected:
    /// Create root SimsNode with no edges.
    SimsTreeBase(
        RankType rank,
        DegreeType max_degree,
        const std::vector<Relator> &short_relators,
        const std::vector<Relator> &long_relators,
        const std::vector<std::vector<DegreeType>> &initial_edges = {},
        bool minimal_pruning = true,
        bool seed_symmetry_pruning = false,
        const FullCycleConstraints &full_cycle_constraints = {});

    // Implements list()
    virtual std::vector<SimsNode> _list() = 0;

    bool _apply_initial_edges(
        const std::vector<std::vector<DegreeType>> &initial_edges);

    void _initialize_seed_relabelings(
        const std::vector<std::vector<DegreeType>> &initial_edges);

    SimsNode _root;
    const std::vector<Relator> _short_relators;
    const std::vector<Relator> _long_relators;
    FullCycleConstraints _full_cycle_constraints;
    const bool _minimal_pruning;
    const bool _seed_symmetry_pruning;

    // Non-identity automorphisms of a connected exact-degree seed,
    // together with their inverses. Empty when no such symmetry
    // is available.
    std::vector<std::vector<DegreeType>> _seed_relabelings;
    std::vector<std::vector<DegreeType>> _seed_inverse_relabelings;

    bool _root_valid;
};

}

#endif
