#include "simsTreeBase.h"

#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace low_index {

SimsTreeBase::SimsTreeBase(
    const RankType rank,
    const DegreeType max_degree,
    const std::vector<Relator> &short_relators,
    const std::vector<Relator> &long_relators,
    const std::vector<std::vector<DegreeType>> &initial_edges,
    const bool minimal_pruning,
    const bool seed_symmetry_pruning,
    const FullCycleConstraints &full_cycle_constraints)
  : _root(rank, max_degree, short_relators.size())
  , _short_relators(short_relators)
  , _long_relators(long_relators)
  , _full_cycle_constraints(full_cycle_constraints)
  , _minimal_pruning(minimal_pruning)
  , _seed_symmetry_pruning(seed_symmetry_pruning)
  , _root_valid(true)
{
    for (const Relator &relator : short_relators) {
        if (!(relator.size() < std::numeric_limits<RelatorLengthType>::max())) {
            throw std::domain_error(
                "Length of a relator can be at most " +
                std::to_string(
                    static_cast<int>(
                        std::numeric_limits<RelatorLengthType>::max())));
        }
    }
    for (const FullCycleConstraint &constraint :
         _full_cycle_constraints) {

        const Relator &word = constraint.first;
        const unsigned int order = constraint.second;

        if (word.empty()) {
            throw std::domain_error(
                "A full-cycle constraint word cannot be empty.");
        }

        if (order == 0) {
            throw std::domain_error(
                "A full-cycle constraint order must be positive.");
        }

        for (const LetterType letter : word) {

            const int value =
                static_cast<int>(letter);

            if (value == 0 ||
                value > static_cast<int>(rank) ||
                value < -static_cast<int>(rank)) {

                throw std::domain_error(
                    "A full-cycle constraint contains an invalid generator.");
            }
        }
    }

    _root_valid = _apply_initial_edges(initial_edges);

    if (_root_valid &&
        !_full_cycle_constraints.empty() &&
        !_root.full_cycle_constraints_may_hold(
            _full_cycle_constraints)) {

        _root_valid = false;
    }

    // If the root already has the maximal degree, no new sheets can
    // appear below it.  A full-cycle constraint whose complete action
    // is already known at the root therefore cannot change later.
    //
    // The root has already passed all constraints above, so such
    // permanently decided constraints may be omitted from descendant
    // checks without changing the search result.
    if (_root_valid &&
        _root.degree() == _root.max_degree() &&
        !_full_cycle_constraints.empty()) {

        _full_cycle_constraints.erase(
            std::remove_if(
                _full_cycle_constraints.begin(),
                _full_cycle_constraints.end(),
                [this](
                    const FullCycleConstraint &constraint)
                {
                    return
                        _root.
                        full_cycle_constraint_is_fully_determined(
                            constraint);
                }),
            _full_cycle_constraints.end());
    }

    if (_root_valid && _seed_symmetry_pruning) {
        _initialize_seed_relabelings(initial_edges);
    }
}

void
SimsTreeBase::_initialize_seed_relabelings(
    const std::vector<std::vector<DegreeType>> &initial_edges)
{
    _seed_relabelings.clear();
    _seed_inverse_relabelings.clear();

    if (initial_edges.empty()) {
        return;
    }

    DegreeType seed_degree = 0;

    for (const std::vector<DegreeType> &edge : initial_edges) {
        if (edge.size() != 3) {
            return;
        }

        seed_degree = std::max<DegreeType>(
            seed_degree,
            std::max<DegreeType>(
                edge[1] + 1,
                edge[2] + 1));
    }

    // For now this optimization is intentionally restricted to
    // exact-degree seeded completion searches. This matches our
    // current targeted-cover workflow.
    if (seed_degree == 0 ||
        seed_degree != _root.max_degree()) {
        return;
    }

    const RankType rank = _root.rank();
    const size_t width =
        static_cast<size_t>(seed_degree) + 1;

    // Store only the explicitly prescribed seed graph here.
    // Do not use relator-deduced edges when computing Aut(seed).
    std::vector<std::vector<DegreeType>> outgoing(
        static_cast<size_t>(rank) + 1,
        std::vector<DegreeType>(width, 0));

    std::vector<std::vector<DegreeType>> incoming(
        static_cast<size_t>(rank) + 1,
        std::vector<DegreeType>(width, 0));

    for (const std::vector<DegreeType> &edge : initial_edges) {
        const LetterType letter =
            static_cast<LetterType>(edge[0]);

        const DegreeType from = edge[1] + 1;
        const DegreeType to = edge[2] + 1;

        if (letter < 1 ||
            letter > rank ||
            from < 1 ||
            from > seed_degree ||
            to < 1 ||
            to > seed_degree) {
            return;
        }

        DegreeType &out = outgoing[letter][from];
        DegreeType &in = incoming[letter][to];

        if ((out != 0 && out != to) ||
            (in != 0 && in != from)) {
            return;
        }

        out = to;
        in = from;
    }

    // Require the prescribed seed graph to be connected when
    // directions are ignored. Then an automorphism is determined
    // completely by the image of vertex 1.
    std::vector<bool> seen(width, false);
    std::vector<DegreeType> queue;

    queue.push_back(1);
    seen[1] = true;

    for (size_t q = 0; q < queue.size(); q++) {
        const DegreeType v = queue[q];

        for (RankType l = 1; l <= rank; l++) {
            for (const DegreeType w : {
                     outgoing[l][v],
                     incoming[l][v] }) {

                if (w != 0 && !seen[w]) {
                    seen[w] = true;
                    queue.push_back(w);
                }
            }
        }
    }

    if (queue.size() !=
        static_cast<size_t>(seed_degree)) {
        return;
    }

    // Try every possible image of vertex 1.
    for (DegreeType image1 = 1;
         image1 <= seed_degree;
         image1++) {

        std::vector<DegreeType> pi(width, 0);
        std::vector<DegreeType> work;

        pi[1] = image1;
        work.push_back(1);

        bool ok = true;

        for (size_t q = 0;
             q < work.size() && ok;
             q++) {

            const DegreeType v = work[q];
            const DegreeType image_v = pi[v];

            for (RankType l = 1;
                 l <= rank && ok;
                 l++) {

                for (const bool forward : {
                         true, false }) {

                    const DegreeType w =
                        forward
                        ? outgoing[l][v]
                        : incoming[l][v];

                    if (w == 0) {
                        continue;
                    }

                    const DegreeType image_w =
                        forward
                        ? outgoing[l][image_v]
                        : incoming[l][image_v];

                    if (image_w == 0) {
                        ok = false;
                        break;
                    }

                    if (pi[w] == 0) {
                        pi[w] = image_w;
                        work.push_back(w);
                    } else if (pi[w] != image_w) {
                        ok = false;
                        break;
                    }
                }
            }
        }

        if (!ok ||
            work.size() !=
                static_cast<size_t>(seed_degree)) {
            continue;
        }

        // Verify pi is actually a permutation and compute pi^-1.
        std::vector<DegreeType> pi_inv(width, 0);

        for (DegreeType v = 1;
             v <= seed_degree;
             v++) {

            if (pi[v] == 0 ||
                pi[v] > seed_degree ||
                pi_inv[pi[v]] != 0) {

                ok = false;
                break;
            }

            pi_inv[pi[v]] = v;
        }

        if (!ok) {
            continue;
        }

        // Explicitly verify every prescribed seed edge is preserved.
        for (const std::vector<DegreeType> &edge :
             initial_edges) {

            const LetterType letter =
                static_cast<LetterType>(edge[0]);

            const DegreeType from = edge[1] + 1;
            const DegreeType to = edge[2] + 1;

            if (outgoing[letter][pi[from]] != pi[to]) {
                ok = false;
                break;
            }
        }

        if (!ok) {
            continue;
        }

        // Identity gives no pruning benefit.
        bool identity = true;

        for (DegreeType v = 1;
             v <= seed_degree;
             v++) {

            if (pi[v] != v) {
                identity = false;
                break;
            }
        }

        if (identity) {
            continue;
        }

        _seed_relabelings.push_back(
            std::move(pi));

        _seed_inverse_relabelings.push_back(
            std::move(pi_inv));
    }
}

bool
SimsTreeBase::_apply_initial_edges(
    const std::vector<std::vector<DegreeType>> &initial_edges)
{
    if (initial_edges.empty()) {
        return true;
    }

    std::vector<bool> added(initial_edges.size(), false);
    size_t num_added = 0;

    while (num_added < initial_edges.size()) {
        bool progressed = false;

        for (size_t i = 0; i < initial_edges.size(); i++) {
            if (added[i]) {
                continue;
            }

            const std::vector<DegreeType> &edge = initial_edges[i];
            if (edge.size() != 3) {
                throw std::domain_error(
                    "Each initial edge must be [letter, from_vertex, to_vertex].");
            }

            const LetterType letter = static_cast<LetterType>(edge[0]);
            const DegreeType from_vertex = edge[1] + 1;
            const DegreeType to_vertex = edge[2] + 1;

            if (letter < 1 || letter > _root.rank()) {
                throw std::domain_error("Initial edge letter is out of range.");
            }
            if (from_vertex < 1 || to_vertex < 1 ||
                from_vertex > _root.max_degree() ||
                to_vertex > _root.max_degree()) {
                throw std::domain_error("Initial edge vertex is out of range.");
            }
            if (from_vertex > _root.degree() + 1 ||
                to_vertex > _root.degree() + 1) {
                continue;
            }

            const DegreeType existing_out = _root.act_by(letter, from_vertex);
            const DegreeType existing_in = _root.act_by(-letter, to_vertex);
            if (existing_out != 0 || existing_in != 0) {
                if (existing_out == to_vertex && existing_in == from_vertex) {
                    added[i] = true;
                    num_added++;
                    progressed = true;
                    continue;
                }
                return false;
            }

            if (!_root.verified_add_edge(letter, from_vertex, to_vertex)) {
                return false;
            }
            added[i] = true;
            num_added++;
            progressed = true;
        }

        if (!progressed) {
            return false;
        }
    }

    return _root.relators_may_lift(_short_relators);
}

SimsTreeBase::~SimsTreeBase() = default;

std::vector<SimsNode>
SimsTreeBase::list() {
    return _list();
}

}
