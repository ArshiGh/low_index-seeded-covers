#ifndef LOW_INDEX_SIMS_TREE_H
#define LOW_INDEX_SIMS_TREE_H

#include "simsTreeBase.h"

namespace low_index {

class RelatorWakeIndex;

/// A single-threaded implementation of SimsTreeBase.

/// It simply recursively adds edges to a SimsNode (starting with the
/// SimsNode with no edges) to find all covering graphs up to the given
/// degree.
///
class SimsTree : public SimsTreeBase
{
public:
    /// See SimsTreeBase for arguments.
    SimsTree(
        RankType rank,
        DegreeType max_degree,
        const std::vector<Relator> &short_relators,
        const std::vector<Relator> &long_relators,
        const std::vector<std::vector<DegreeType>> &initial_edges = {},
        bool minimal_pruning = true,
        bool seed_symmetry_pruning = false,
        const FullCycleConstraints &full_cycle_constraints = {});

protected:
    std::vector<SimsNode> _list() override;

private:
    void _recurse(
        const class StackedSimsNode &n,
        RelatorWakeIndex &wake_index);

    std::vector<SimsNode> _complete_nodes;
};

} // Namespace low_index

#endif
