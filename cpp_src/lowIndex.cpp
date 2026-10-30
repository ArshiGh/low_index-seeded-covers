#include "lowIndex.h"

#include "words.h"
#include "simsTree.h"
#include "simsTreeMultiThreaded.h"

#include <thread>
#include <memory>

namespace low_index {

const std::string spin_short_strategy = "spin_short";

static
std::vector<std::vector<std::vector<DegreeType>>>
permutation_reps_impl(
    const RankType rank,
    const std::vector<Relator> &short_relators,
    const std::vector<Relator> &long_relators,
    const DegreeType max_degree,
    const std::string &strategy,
    const unsigned int num_threads,
    const std::vector<std::vector<DegreeType>> &initial_edges,
    const bool minimal_pruning,
    const bool seed_symmetry_pruning,
    const FullCycleConstraints &full_cycle_constraints)
{
    // Apply strategy to short relators.
    const std::vector<Relator> all_short_relators =
        (strategy == spin_short_strategy)
            ? spin_short(short_relators, max_degree)
            : short_relators;

    // Determine number of threads to use
    const unsigned int resolved_num_threads =
        (num_threads > 0)
            ? num_threads
            : std::thread::hardware_concurrency();

    // Instantiate appropriate SimsTree implementation
    std::unique_ptr<SimsTreeBase> t;
    if (resolved_num_threads > 1) {
        t = std::make_unique<SimsTreeMultiThreaded>(
                rank, max_degree, all_short_relators, long_relators,
                resolved_num_threads, initial_edges, minimal_pruning,
                seed_symmetry_pruning, full_cycle_constraints);
    } else {
        t = std::make_unique<SimsTree>(
                rank, max_degree, all_short_relators, long_relators,
                initial_edges, minimal_pruning,
                seed_symmetry_pruning, full_cycle_constraints);
    }

    // Convert SimsNode's to permutation representations.
    std::vector<std::vector<std::vector<DegreeType>>> result;
    for (const SimsNode &n : t->list()) {
        result.push_back(n.permutation_rep());
    }
    return result;
}

std::vector<std::vector<std::vector<DegreeType>>>
permutation_reps(
    const RankType rank,
    const std::vector<Relator> &short_relators,
    const std::vector<Relator> &long_relators,
    const DegreeType max_degree,
    const std::string &strategy,
    const unsigned int num_threads)
{
    return permutation_reps_impl(
        rank, short_relators, long_relators, max_degree, strategy,
        num_threads, {}, true, false, {});
}

std::vector<std::vector<std::vector<DegreeType>>>
permutation_reps_from_seed(
    const RankType rank,
    const std::vector<Relator> &short_relators,
    const std::vector<Relator> &long_relators,
    const DegreeType max_degree,
    const std::vector<std::vector<DegreeType>> &initial_edges,
    const std::string &strategy,
    const unsigned int num_threads,
    const bool minimal_pruning,
    const bool seed_symmetry_pruning,
    const FullCycleConstraints &full_cycle_constraints)
{
    return permutation_reps_impl(
        rank, short_relators, long_relators, max_degree, strategy,
        num_threads, initial_edges, minimal_pruning,
        seed_symmetry_pruning, full_cycle_constraints);
}

// Parse a list of SnapPy-words
static
std::vector<Relator>
parse_words(
    const RankType rank,
    const std::vector<std::string> &words)
{
    std::vector<Relator> result;
    result.reserve(words.size());
    for (const std::string &word : words) {
        result.push_back(parse_word(rank, word));
    }
    return result;
}

std::vector<std::vector<std::vector<DegreeType>>>
permutation_reps(
    const RankType rank,
    const std::vector<std::string> &short_relators,
    const std::vector<std::string> &long_relators,
    const DegreeType max_degree,
    const std::string &strategy,
    const unsigned int num_threads)
{
    return permutation_reps(
        rank,
        parse_words(rank, short_relators),
        parse_words(rank, long_relators),
        max_degree,
        strategy,
        num_threads);
}

std::vector<std::vector<std::vector<DegreeType>>>
permutation_reps_from_seed(
    const RankType rank,
    const std::vector<std::string> &short_relators,
    const std::vector<std::string> &long_relators,
    const DegreeType max_degree,
    const std::vector<std::vector<DegreeType>> &initial_edges,
    const std::string &strategy,
    const unsigned int num_threads,
    const bool minimal_pruning,
    const bool seed_symmetry_pruning,
    const FullCycleConstraints &full_cycle_constraints)
{
    return permutation_reps_from_seed(
        rank,
        parse_words(rank, short_relators),
        parse_words(rank, long_relators),
        max_degree,
        initial_edges,
        strategy,
        num_threads,
        minimal_pruning,
        seed_symmetry_pruning,
        full_cycle_constraints);
}

}
