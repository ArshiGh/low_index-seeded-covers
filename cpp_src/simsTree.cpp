#include "simsTree.h"
#include "stackedSimsNode.h"

namespace low_index {

SimsTree::SimsTree(
    const RankType rank,
    const DegreeType max_degree,
    const std::vector<Relator> &short_relators,
    const std::vector<Relator> &long_relators,
    const std::vector<std::vector<DegreeType>> &initial_edges,
    const bool minimal_pruning,
    const bool seed_symmetry_pruning,
    const FullCycleConstraints &full_cycle_constraints)
  : SimsTreeBase(
        rank, max_degree, short_relators, long_relators,
        initial_edges, minimal_pruning, seed_symmetry_pruning,
        full_cycle_constraints)
{
}

std::vector<SimsNode>
SimsTree::_list()
{
    if (!_root_valid) {
        return {};
    }
    // Allocate all memory needed to recurse up front.
    SimsNodeStack stack(_root);

    RelatorWakeIndex wake_index(
        stack.get_node(),
        _short_relators);

    _recurse(
        stack.get_node(),
        wake_index);
    return std::move(_complete_nodes);
}

void
SimsTree::_recurse(
    const StackedSimsNode &n,
    RelatorWakeIndex &wake_index)
{
    if (n.is_complete()) {
        if (!n.relators_lift(
                _long_relators)) {
            return;
        }

        SimsNode copy(n);

        if (!copy.relators_may_lift(
                _short_relators)) {
            return;
        }

        _complete_nodes.push_back(
            std::move(copy));

        return;
    }

    const std::pair<
        LetterType,
        DegreeType
    > slot =
        n.first_empty_slot();

    const DegreeType m =
        std::min<DegreeType>(
            n.degree() + 1,
            n.max_degree());

    for (DegreeType v = 1;
         v <= m;
         v++) {

        if (n.act_by(
                -slot.first,
                v) != 0) {

            continue;
        }

        const RelatorWakeIndex::Checkpoint
            checkpoint =
                wake_index.checkpoint();

        StackedSimsNode new_subgraph(n);

        new_subgraph.add_edge(
            slot.first,
            slot.second,
            v);

        if (!new_subgraph.relators_may_lift_after_edge(
                _short_relators,
                slot.first,
                slot.second,
                v,
                wake_index)) {

            wake_index.rollback(
                checkpoint);

            continue;
        }

        if (!_full_cycle_constraints.empty() &&
            !new_subgraph.full_cycle_constraints_may_hold(
                _full_cycle_constraints)) {

            wake_index.rollback(
                checkpoint);

            continue;
        }

        if (_minimal_pruning &&
            !new_subgraph.may_be_minimal()) {

            wake_index.rollback(
                checkpoint);

            continue;
        }

        if (!_seed_relabelings.empty() &&
            !new_subgraph.may_be_minimal_under_relabelings(
                _seed_relabelings,
                _seed_inverse_relabelings)) {

            wake_index.rollback(
                checkpoint);

            continue;
        }

        _recurse(
            new_subgraph,
            wake_index);

        wake_index.rollback(
            checkpoint);
    }
}


} // Namespace low_index
