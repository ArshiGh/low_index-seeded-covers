#include "abstractSimsNode.h"

#include <stdexcept>

#include <limits>
#include <stdexcept>
#include <cstring>

namespace low_index {

AbstractSimsNode::AbstractSimsNode(
        const RankType rank,
        const DegreeType max_degree,
        const unsigned int num_relators)
 : CoveringSubgraph(rank, max_degree)
 , _num_relators(num_relators)
{
}

AbstractSimsNode::AbstractSimsNode(
    const AbstractSimsNode &other)
 : CoveringSubgraph(other)
 , _num_relators(other._num_relators)
{
}

// Make sure n is a multiple of the alignment of T.
template<typename T>
static
size_t
_align(const size_t n)
{
    return ((n + alignof(T) - 1 ) / alignof(T)) * alignof(T);
}

AbstractSimsNode::_MemoryLayout::_MemoryLayout(
    const AbstractSimsNode &node)
{
    // DegreeType *outgoing;
    size_t t = outgoing_offset;
    t += node.rank() * node.max_degree() * sizeof(DegreeType);

    // DegreeType *incoming;
    t = _align<DegreeType>(t);
    incoming_offset = t;
    t += node.rank() * node.max_degree() * sizeof(DegreeType);

    // RelatorLengthType *lift_indices;
    t = _align<RelatorLengthType>(t);
    lift_indices_offset = t;
    t += node.num_relators() * node.max_degree() * sizeof(RelatorLengthType);

    // DegreeType *lift_vertices;
    t = _align<DegreeType>(t);
    lift_vertices_offset = t;
    t += node.num_relators() * node.max_degree() * sizeof(DegreeType);

    size = _align<uint64_t>(t);
}

void
AbstractSimsNode::_apply_memory_layout(
    const _MemoryLayout &layout,
    uint8_t * const memory)
{
    _outgoing =
        reinterpret_cast<DegreeType*>(
            memory + layout.outgoing_offset);
    _incoming =
        reinterpret_cast<DegreeType*>(
            memory + layout.incoming_offset);
    _lift_indices =
        reinterpret_cast<RelatorLengthType*>(
            memory + layout.lift_indices_offset);
    _lift_vertices =
        reinterpret_cast<DegreeType*>(
            memory + layout.lift_vertices_offset);
    _memory_size = layout.size;
}

void
AbstractSimsNode::_initialize_memory()
{
    std::memset(_memory_start(), 0, _memory_size);

    // The lift of vertex v + 1 by the empty word is just v + 1.
    // Use this to initialize _lift_vertices.
    for (size_t n = 0; n < _num_relators; n++) {
        for (DegreeType v = 0; v < max_degree(); v++) {
            const size_t j = n * max_degree() + v;
            _lift_vertices[j] = v + 1;
        }
    }
}

void
AbstractSimsNode::_copy_memory(const AbstractSimsNode &other)
{
    std::memcpy(_memory_start(), other._memory_start(), _memory_size);
}


RelatorWakeIndex::RelatorWakeIndex(
    const AbstractSimsNode &node,
    const std::vector<Relator> &relators)
  : _relators(&relators)
  , _rank(node.rank())
  , _max_degree(node.max_degree())
  , _active_degree(node.degree())
  , _num_states(
        static_cast<size_t>(node.num_relators()) *
        static_cast<size_t>(node.max_degree()))
  , _num_slots(
        2 *
        static_cast<size_t>(node.rank()) *
        static_cast<size_t>(node.max_degree()))
  , _words_per_bucket(
        (_num_states + 63) / 64)
  , _no_slot(
        std::numeric_limits<size_t>::max())
  , _buckets(
        _num_slots * _words_per_bucket,
        uint64_t(0))
  , _state_slots(
        _num_states,
        _no_slot)
{
    if (relators.size() != node.num_relators()) {
        throw std::invalid_argument(
            "RelatorWakeIndex: relator count mismatch.");
    }

    _undo.reserve(_num_states);

    _events.reserve(
        2 * (
            static_cast<size_t>(_rank) *
            static_cast<size_t>(_max_degree)
            + 1));

    constexpr DegreeType finished =
        std::numeric_limits<DegreeType>::max();

    for (size_t n = 0;
         n < relators.size();
         n++) {

        for (DegreeType v = 0;
             v < _active_degree;
             v++) {

            const size_t j =
                n * static_cast<size_t>(_max_degree)
                + static_cast<size_t>(v);

            const DegreeType end_vertex =
                node._lift_vertices[j];

            if (end_vertex == finished) {
                continue;
            }

            const RelatorLengthType lift_index =
                node._lift_indices[j];

            const size_t slot =
                _slot_key(
                    relators[n][lift_index],
                    end_vertex);

            _state_slots[j] = slot;

            _set_bit(
                j,
                slot);
        }
    }
}


size_t
RelatorWakeIndex::_slot_key(
    const LetterType letter,
    const DegreeType vertex) const
{
    const int signed_letter =
        static_cast<int>(letter);

    const int abs_letter =
        signed_letter > 0
            ? signed_letter
            : -signed_letter;

    if (signed_letter == 0 ||
        abs_letter > static_cast<int>(_rank) ||
        vertex < 1 ||
        vertex > _max_degree) {

        throw std::domain_error(
            "RelatorWakeIndex: invalid directed slot.");
    }

    const size_t letter_index =
        2 * static_cast<size_t>(abs_letter - 1)
        + static_cast<size_t>(
            signed_letter < 0);

    return
        letter_index *
        static_cast<size_t>(_max_degree)
        + static_cast<size_t>(vertex - 1);
}


void
RelatorWakeIndex::_set_bit(
    const size_t state,
    const size_t slot)
{
    const size_t word =
        state / 64;

    const unsigned int bit =
        static_cast<unsigned int>(
            state % 64);

    _buckets[
        slot * _words_per_bucket
        + word
    ] |= (uint64_t(1) << bit);
}


void
RelatorWakeIndex::_clear_bit(
    const size_t state,
    const size_t slot)
{
    const size_t word =
        state / 64;

    const unsigned int bit =
        static_cast<unsigned int>(
            state % 64);

    _buckets[
        slot * _words_per_bucket
        + word
    ] &= ~(uint64_t(1) << bit);
}


void
RelatorWakeIndex::_move_state(
    const size_t state,
    const size_t new_slot)
{
    const size_t old_slot =
        _state_slots[state];

    if (old_slot == new_slot) {
        return;
    }

    _undo.push_back(
        _Undo{
            state,
            old_slot,
            new_slot
        });

    if (old_slot != _no_slot) {
        _clear_bit(
            state,
            old_slot);
    }

    if (new_slot != _no_slot) {
        _set_bit(
            state,
            new_slot);
    }

    _state_slots[state] =
        new_slot;
}


void
RelatorWakeIndex::_activate_degree(
    const AbstractSimsNode &node)
{
    if (node.degree() < _active_degree) {
        throw std::logic_error(
            "RelatorWakeIndex: degree moved backwards.");
    }

    constexpr DegreeType finished =
        std::numeric_limits<DegreeType>::max();

    for (DegreeType v = _active_degree;
         v < node.degree();
         v++) {

        for (size_t n = 0;
             n < _relators->size();
             n++) {

            const size_t j =
                n * static_cast<size_t>(_max_degree)
                + static_cast<size_t>(v);

            const DegreeType end_vertex =
                node._lift_vertices[j];

            if (end_vertex == finished) {
                continue;
            }

            const RelatorLengthType lift_index =
                node._lift_indices[j];

            const size_t slot =
                _slot_key(
                    (*_relators)[n][lift_index],
                    end_vertex);

            _move_state(
                j,
                slot);
        }
    }

    _active_degree =
        node.degree();
}


void
RelatorWakeIndex::_push_edge_events(
    const LetterType letter,
    const DegreeType source,
    const DegreeType target)
{
    _events.push_back(
        _Event{
            letter,
            source
        });

    _events.push_back(
        _Event{
            static_cast<LetterType>(
                -static_cast<int>(letter)),
            target
        });
}


void
RelatorWakeIndex::rollback(
    const Checkpoint &checkpoint)
{
    if (checkpoint.undo_size > _undo.size()) {
        throw std::logic_error(
            "RelatorWakeIndex: invalid checkpoint.");
    }

    while (_undo.size() >
           checkpoint.undo_size) {

        const _Undo change =
            _undo.back();

        _undo.pop_back();

        if (change.new_slot != _no_slot) {
            _clear_bit(
                change.state,
                change.new_slot);
        }

        if (change.old_slot != _no_slot) {
            _set_bit(
                change.state,
                change.old_slot);
        }

        _state_slots[change.state] =
            change.old_slot;
    }

    _active_degree =
        checkpoint.active_degree;
}



bool
AbstractSimsNode::full_cycle_constraints_may_hold(
    const FullCycleConstraints &constraints) const
{
    for (const FullCycleConstraint &constraint : constraints) {

        const Relator &word = constraint.first;
        const unsigned int order = constraint.second;

        if (word.empty() || order == 0) {
            return false;
        }

        for (DegreeType start = 1;
             start <= degree();
             start++) {

            DegreeType vertex = start;

            for (unsigned int repetition = 1;
                 repetition <= order;
                 repetition++) {

                bool known = true;

                for (const LetterType letter : word) {

                    vertex = act_by(
                        letter,
                        vertex);

                    if (vertex == 0) {
                        known = false;
                        break;
                    }
                }

                // Open path: later edges may still complete a valid cycle.
                if (!known) {
                    break;
                }

                if (vertex == start) {

                    // Closing before the required order is irreversible.
                    if (repetition != order) {
                        return false;
                    }

                    break;
                }

                // The entire word^order path is already known but failed
                // to return to its starting sheet.
                if (repetition == order) {
                    return false;
                }
            }
        }
    }

    return true;
}


bool
AbstractSimsNode::full_cycle_constraint_is_fully_determined(
    const FullCycleConstraint &constraint) const
{
    const Relator &word = constraint.first;
    const unsigned int order = constraint.second;

    if (word.empty() || order == 0) {
        return false;
    }

    for (DegreeType start = 1;
         start <= degree();
         start++) {

        DegreeType vertex = start;

        for (unsigned int repetition = 1;
             repetition <= order;
             repetition++) {

            for (const LetterType letter : word) {

                vertex = act_by(
                    letter,
                    vertex);

                if (vertex == 0) {
                    return false;
                }
            }

            // Once the orbit returns to its starting sheet, the
            // full-cycle checker needs no further information for
            // this start.  The caller has already verified that any
            // such closure occurs at the required order.
            if (vertex == start) {
                break;
            }
        }
    }

    return true;
}



bool
AbstractSimsNode::relators_may_lift(const std::vector<Relator> &relators)
{
    for (size_t n = 0; n < relators.size(); n++) {
        for (DegreeType v = 0; v < degree(); v++) {
	    if (!_relator_may_lift(relators[n], n, v)) {
                return false;
	    }
        }
    }
    return true;
}

bool
AbstractSimsNode::relators_may_lift_after_edge(
    const std::vector<Relator> &relators,
    const LetterType letter,
    const DegreeType source,
    const DegreeType target,
    RelatorWakeIndex &wake_index)
{
    if (relators.size() != num_relators()) {
        throw std::invalid_argument(
            "relators_may_lift_after_edge: "
            "relator count mismatch.");
    }

    // An explicit DFS edge can activate one newly introduced vertex.
    // Register its cached lift states before delivering this edge.
    wake_index._activate_degree(*this);

    wake_index._events.clear();

    wake_index._push_edge_events(
        letter,
        source,
        target);

    constexpr DegreeType finished =
        std::numeric_limits<DegreeType>::max();

    for (size_t event_pos = 0;
         event_pos < wake_index._events.size();
         event_pos++) {

        // _relator_may_lift can append forced-edge events.
        const RelatorWakeIndex::_Event event =
            wake_index._events[event_pos];

        const size_t slot =
            wake_index._slot_key(
                event.letter,
                event.vertex);

        const size_t bucket_base =
            slot *
            wake_index._words_per_bucket;

        for (size_t word = 0;
             word < wake_index._words_per_bucket;
             word++) {

            // Snapshot this word for this event.
            uint64_t pending =
                wake_index._buckets[
                    bucket_base + word];

            while (pending != 0) {

                const unsigned int bit =
                    static_cast<unsigned int>(
                        __builtin_ctzll(
                            pending));

                pending &=
                    pending - 1;

                const size_t j =
                    word * 64
                    + static_cast<size_t>(
                        bit);

                if (j >= wake_index._num_states) {
                    continue;
                }

                // The state may already have moved because of an
                // earlier event in this deduction cascade.
                if (wake_index._state_slots[j] != slot) {
                    continue;
                }

                const size_t n =
                    j /
                    static_cast<size_t>(
                        max_degree());

                const DegreeType v =
                    static_cast<DegreeType>(
                        j %
                        static_cast<size_t>(
                            max_degree()));

                if (v >= degree()) {
                    continue;
                }

                if (!_relator_may_lift(
                        relators[n],
                        n,
                        v,
                        &wake_index)) {

                    return false;
                }

                size_t new_slot =
                    wake_index._no_slot;

                const DegreeType new_end_vertex =
                    _lift_vertices[j];

                if (new_end_vertex != finished) {

                    const RelatorLengthType
                        new_lift_index =
                            _lift_indices[j];

                    new_slot =
                        wake_index._slot_key(
                            relators[n][new_lift_index],
                            new_end_vertex);
                }

                wake_index._move_state(
                    j,
                    new_slot);
            }
        }
    }

    return true;
}


bool
AbstractSimsNode::_relator_may_lift(
    const Relator &relator,
    const size_t n,
    const DegreeType v,
    RelatorWakeIndex * const wake_index)
{
    const size_t j = n * max_degree() + v;

    constexpr DegreeType finished =
        std::numeric_limits<DegreeType>::max();

    DegreeType vertex = _lift_vertices[j];
    // We already determined in an earlier run of _relator_may_lift
    // that this relator lifts.
    if (vertex == finished) {
        return true;
    }

    // Continue lifting the relator where we left of.
    DegreeType next_vertex;
    for (RelatorLengthType i = _lift_indices[j]; true; i++) {
        // Result of lifting the edge given by the next letter in
        // the relator.
        next_vertex = act_by(relator[i], vertex);
        if (i == relator.size() - 1) {
            // We are at the last letter of the relator.
            // This case is handled specially below.
            break;
        }
        if (next_vertex == 0) {
            // The is no edge yet corresponding to the next letter by which
            // we lift the vertex. Store how far we were able to lift the
            // relator for the next call to _relator_may_lift.
            _lift_vertices[j] = vertex;
            _lift_indices[j] = i;
            return true;
        }
        // Move on to the next vertex before looking at the next
        // letter in the relator.
        vertex = next_vertex;
    }

    // We are at the last letter in the relator.

    if (next_vertex == v + 1) {
        // We were able to lift the relator to a loop.  Record this fact and
        // return true.
        _lift_vertices[j] = finished;
        return true;
    }

    if (next_vertex == 0) {
        // The relator almost lifted completely, but the last edge that should
        // appear in the lift was missing from the graph.  The only way that the
        // relator could lift to a loop would be if the missing edge joined
        // the last vertex of the lift to the initial vertex v+1. We attempt
        // to add an edge using the empty slot at the end of the lift and the
        // appropriate slot at the initial vertex v + 1.
        //
        // Note that it is only possible to add such an edge if both of those
        // slots are empty.  We know that the slot at the end of the lift is
        // empty, but the slot at the initial vertex might have already been
        // filled by an edge that was added earlier.  So we must call
        // verified_add_edge here to avoid corrupting the structure of the
        // CoveringSubgraph.
        if (verified_add_edge(relator.back(), vertex, v + 1)) {
            if (wake_index != nullptr) {
                wake_index->_push_edge_events(
                    relator.back(),
                    vertex,
                    static_cast<DegreeType>(v + 1));
            }
            // Record that the relator lifts to a loop and return true.
            _lift_vertices[j] = finished;
            return true;
        }
    }

    // In all other cases return false to signal that the relator does not lift
    // to a loop.
    return false;
}

bool
AbstractSimsNode::relators_lift(const std::vector<Relator> &relators) const
{
    for (const Relator &relator : relators) {
        for (DegreeType v = 1; v <= degree(); v++) {
            // Start with vertex v.
            DegreeType vertex = v;
            for (const int letter : relator) {
                // Traverse the edges labeled by the letters in the relator.
                vertex = act_by(letter, vertex);
                if (vertex == 0) {
                    throw std::domain_error(
                        "relators_lift: The graph is not a covering.");
                }
            }
            if (vertex != v) {
                return false;
            }
        }
    }

    return true;
}

bool
AbstractSimsNode::may_be_minimal() const
{
    // Return False if the subgraph can provably not be extended to a cover
    // which is minimal in its conjugacy class, True otherwise.
    //
    // The ordering of based coverings is given as follows.  Given a covering
    // with rank R and degree D, consider all pairs (N, L) where N is a vertex
    // index in [1, D] and L is a signed label in [-R, -1] or [1, R].  Order
    // the signed labels as 1, -1, 2, -2, ..., R, -R.  Use that to lex order
    // the 2*D*R pairs (N, L).  For each pair (N, L) assign an edge E(N, L) and
    // a number s(N, L) as follows:
    //
    //     * if L > 0, E(N, L) is the edge with initial vertex N and
    //       label L and s(N, L) is index of its terminal vertex.
    //
    //     * if L < 0, E(N, L) is the edge with terminal vertex N and
    //       label -L and s(N, L) is the index of its initial vertex.
    //
    // The list [s(N, L) for all (N, L)] is the complexity of the covering;
    // complexities are ordered lexicographically.
    //
    // Any conjugate of the subgroup corresponding to a covering graph can be
    // obtained by changing the basepoint.  A choice of basepoint determines an
    // indexing of the vertices where the index of a vertex v is the next
    // available index at the moment when the first E(N, L) having v as its
    // terminal (initial) vertex if L > 0 (L < 0) is added.  We only enumerate
    // based coverings for which the basepoint is chosen to minimize
    // complexity.  The subgroups corresponding to these coverings will then by
    // unique up to conjugacy.
    //
    // Even when working with a proper subgraph of a covering graph (containing
    // the basepoint) it may be possible to deduce that the subgraph cannot be
    // extended to a minimal covering.  In such a case the subgraph should not
    // be added to the Sims tree.  If this test is done at each stage of the
    // construction of the tree, the complete graphs in the final tree will
    // automatically be minimal.
    //
    // This method iterates through all vertices of the subgraph, constructing
    // as much as possible of the vertex ordering determined by taking that
    // vertex as a basepoint.  It returns False as soon as it encounters an
    // edge which would result in a higher complexity. If no such edge is found
    // for any choice of basepoint it returns True.

    for (DegreeType basepoint = 2; basepoint <= degree(); basepoint++) {
        if (!_may_be_minimal(basepoint)) {
            return false;
        }
    }
    return true;
}

bool
AbstractSimsNode::may_be_minimal_under_relabelings(
    const std::vector<std::vector<DegreeType>> &relabelings,
    const std::vector<std::vector<DegreeType>> &inverse_relabelings) const
{
    if (relabelings.size() != inverse_relabelings.size()) {
        throw std::domain_error(
            "Seed relabelings and inverse relabelings have different sizes.");
    }

    for (size_t i = 0; i < relabelings.size(); i++) {
        const std::vector<DegreeType> &pi = relabelings[i];
        const std::vector<DegreeType> &pi_inv = inverse_relabelings[i];

        // Compare the current partial graph with
        // q = pi * p * pi^-1
        // using the same lexicographic slot order as may_be_minimal().
        //
        // If an earlier slot is unresolved, no safe pruning decision
        // can yet be made for this relabeling.
        bool comparison_finished = false;

        for (DegreeType slot_vertex = 1;
             slot_vertex <= degree() && !comparison_finished;
             slot_vertex++) {

            const DegreeType source = pi_inv[slot_vertex];

            for (RankType l = 1;
                 l <= rank() && !comparison_finished;
                 l++) {

                for (const LetterType letter : {
                         static_cast<LetterType>(l),
                         static_cast<LetterType>(
                             -static_cast<LetterType>(l)) }) {

                    const DegreeType a =
                        act_by(letter, slot_vertex);

                    const DegreeType old_b =
                        act_by(letter, source);

                    if (a == 0 || old_b == 0) {
                        comparison_finished = true;
                        break;
                    }

                    const DegreeType b = pi[old_b];

                    if (b < a) {
                        // A seed-equivalent image is already
                        // provably lexicographically smaller.
                        return false;
                    }

                    if (b > a) {
                        // The current graph is already smaller
                        // than this particular seed image.
                        comparison_finished = true;
                        break;
                    }
                }
            }
        }
    }

    return true;
}

bool
AbstractSimsNode::_may_be_minimal(const DegreeType basepoint) const
{
    // Unfortunately, MS Visual Studio does not support
    // C99-style Variable Length Arrays, so we allocate a fixed length array.
    constexpr size_t m = std::numeric_limits<DegreeType>::max() + 1;

    static_assert(m < 1000,
                  "Allocating too large array on the stack.");
    
    // We are working with the standard indexing (determined by putting
    // the basepoint at vertex 1) and an alternate indexing determined by
    // a different basepoint.  We construct mappings between the two
    // indexings and store them in the arrays std_to_alt and
    // alt_to_std. (For convenience when dealing with 1-based indices,
    // just ignore the 0 entry).
    DegreeType std_to_alt[m];
    std::memset(std_to_alt, 0, sizeof(DegreeType) * (degree() + 1));

    DegreeType alt_to_std[m];
    // It is not necessary to clear alt_to_std
    std::memset(alt_to_std, 0, sizeof(DegreeType) * (degree() + 1));

    // Initial state.
    std_to_alt[basepoint] = 1;
    alt_to_std[1] = basepoint;

    DegreeType max_index = 1;

    // Iterate over all possible slots
    for (DegreeType slot_vertex = 1; slot_vertex <= degree(); slot_vertex++) {
        for (RankType l = 0; l < rank(); l++) {
            for (const DegreeType * const edges : { _outgoing, _incoming }) {
                // Check that the slot is filled with repect to both indexings.
                const DegreeType a = edges[
                    (slot_vertex - 1) * rank() + l];
                const DegreeType b = edges[
                    (alt_to_std[slot_vertex] - 1) * rank() + l];
                if (a == 0 || b == 0) {
                    // The slot was empty in one indexing, so we cannot decide.
                    return true;
                }
                // Update the mappings.
                DegreeType &c = std_to_alt[b];
                if (c == 0) {
                    // This edge is the first, with respect to the alternate
                    // indexing, that is incident to the vertex with standard
                    // index b.  We now know its alternate index.
                    max_index++;
                    c = max_index;
                    alt_to_std[max_index] = b;
                }
                // Compare the old and new indices of the other end of the edge.
                if (c < a) {
                    // The new basepoint is better - discard this graph.
                    return false;
                }
                if (c > a) {
                    // The old basepoint is better.
                    return true;
                }
            }
        }
    }

    return true;
}

} // Namespace low_index
