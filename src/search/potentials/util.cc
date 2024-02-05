#include "util.h"

#include "potential_function.h"
#include "potential_optimizer.h"

#include "../heuristic.h"

#include "../plugins/plugin.h"
#include "../task_utils/sampling.h"
#include "../utils/markup.h"

#include <limits>

using namespace std;

namespace potentials {
vector<State> sample_without_dead_end_detection(
    PotentialOptimizer &optimizer,
    int num_samples,
    utils::RandomNumberGenerator &rng) {
    const shared_ptr<AbstractTask> task = optimizer.get_task();
    const TaskProxy task_proxy(*task);
    State initial_state = task_proxy.get_initial_state();
    optimizer.optimize_for_state(initial_state);
    int init_h = optimizer.get_potential_function()->get_value(initial_state);
    sampling::RandomWalkSampler sampler(task_proxy, rng);
    vector<State> samples;
    samples.reserve(num_samples);
    for (int i = 0; i < num_samples; ++i) {
        samples.push_back(sampler.sample_state(init_h));
    }
    return samples;
}

string get_admissible_potentials_reference() {
    return "The algorithm is based on" + utils::format_conference_reference(
        {"Jendrik Seipp", "Florian Pommerening", "Malte Helmert"},
        "New Optimization Functions for Potential Heuristics",
        "https://ai.dmi.unibas.ch/papers/seipp-et-al-icaps2015.pdf",
        "Proceedings of the 25th International Conference on"
        " Automated Planning and Scheduling (ICAPS 2015)",
        "193-201",
        "AAAI Press",
        "2015");
}

void prepare_parser_for_admissible_potentials(plugins::Feature &feature) {
    feature.document_language_support("action costs", "supported");
    feature.document_language_support("conditional effects", "not supported");
    feature.document_language_support("axioms", "not supported");
    feature.document_property("admissible", "yes");
    feature.document_property("consistent", "yes");
    feature.document_property("safe", "yes");
    feature.document_property("preferred operators", "no");
    feature.add_option<double>(
        "max_potential",
        "Bound potentials by this number. Using the bound {{{infinity}}} "
        "disables the bounds. In some domains this makes the computation of "
        "weights unbounded in which case no weights can be extracted. Using "
        "very high weights can cause numerical instability in the LP solver, "
        "while using very low weights limits the choice of potential "
        "heuristics. For details, see the ICAPS paper cited above.",
        "1e8",
        plugins::Bounds("0.0", "infinity"));
    
    feature.add_option<bool>(
        "use_presolve",
        "turn presolving on or off. Using presolve creates an overhead "
        "so turning presolve off might decrease runtime.",
        "true");

    feature.add_option<int>(
        "lp_solve_method",
        "determine which method the LP solver should use to solve the "
        "LPs given when using standard CPLEX (not twophase).",
        "0"
    );

    feature.add_option<bool>(
        "use_presolve",
        "turn presolving on or off. Using presolve creates an overhead "
        "so turning presolve off might decrease runtime.",
        "true");

    feature.add_option<bool>(
        "crossover",
        "when set to false, turns off crossover when using the barrier optimizer. "
        "Crossover finds a basis, which is costly, but required for a warm start. "
        "Because of the warm start, it is turned on by default.",
        "true");
    
    feature.add_option<int>(
        "folding_level",
        "turn folding_level off or on. -1 is automatic, 0 is off, "
        "1 to 5 are the different levels of folding, where "
        "1 is a moderate and 5 is an extremely aggressive level of folding.",
        "-1");
    
    feature.add_option<int>(
        "solve_dual",
        "sets whether this the solver should solve the dual formulation "
        "of the LP instead of the Primal LP. Per default this is set to "
        "0 (automatic), which will usually solve the primal LP. Setting "
        "it to 1 will turn it on, and setting it to -1 will turn it off "
        "entirely.",
        "0"
        );

    lp::add_lp_solver_option_to_feature(feature);
    Heuristic::add_options_to_feature(feature);
}
}
