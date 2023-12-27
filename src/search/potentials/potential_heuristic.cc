#include "potential_heuristic.h"

#include "potential_function.h"

#include "../plugins/plugin.h"

using namespace std;

namespace potentials {
PotentialHeuristic::PotentialHeuristic(
    const plugins::Options &opts, unique_ptr<PotentialFunction> function)
    : Heuristic(opts),
    lp_solver(opts.get<lp::LPSolverType>("lpsolver")),
      function(move(function)) {
}

PotentialHeuristic::~PotentialHeuristic() {
}

void PotentialHeuristic::print_statistics() const
{
    utils::g_log << "LP statistics for potential heuristic:" << endl;
    lp_solver.print_statistics();
}

int PotentialHeuristic::compute_heuristic(const State &ancestor_state) {
    State state = convert_ancestor_state(ancestor_state);
    return max(0, function->get_value(state));
}
}
