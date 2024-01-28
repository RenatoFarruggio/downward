#include "operator_counting_heuristic.h"

#include "constraint_generator.h"

#include "../plugins/plugin.h"
#include "../utils/markup.h"
#include "../utils/logging.h"

#include <cmath>

using namespace std;

namespace operator_counting {
OperatorCountingHeuristic::OperatorCountingHeuristic(const plugins::Options &opts)
    : Heuristic(opts),
      constraint_generators(
          opts.get_list<shared_ptr<ConstraintGenerator>>("constraint_generators")),
      lp_solver(opts.get<lp::LPSolverType>("lpsolver")),
      use_integer_operator_counts(opts.get<bool>("use_integer_operator_counts")),
      use_presolve(opts.get<bool>("use_presolve")),
      use_crossover(opts.get<bool>("crossover")),
      symmetry_breaking_level(opts.get<int>("symmetry_breaking_level")),
      folding_level(opts.get<int>("folding_level")),
      save_presolved_lp(opts.get<bool>("save_presolved_lp")),
      use_warm_starts(opts.get<bool>("use_warm_starts")),
      lp_solve_method_id(opts.get<int>("lp_solve_method")),
      is_first(true) {
    lp_solver.set_mip_gap(0);
    named_vector::NamedVector<lp::LPVariable> variables;
    double infinity = lp_solver.get_infinity();
    for (OperatorProxy op : task_proxy.get_operators()) {
        int op_cost = op.get_cost();
        variables.push_back(lp::LPVariable(0, infinity, op_cost, use_integer_operator_counts));
    }
    lp::LinearProgram lp(lp::LPObjectiveSense::MINIMIZE, move(variables), {}, infinity);
    for (const auto &generator : constraint_generators) {
        generator->initialize_constraints(task, lp);
    }
    lp_solver.set_use_presolve(use_presolve);
    lp_solver.set_crossover(use_crossover);
    lp_solver.set_symmetry_breaking(symmetry_breaking_level);
    lp_solver.set_folding_level(folding_level);
    lp_solver.set_save_presolved_lp(save_presolved_lp);
    lp_solver.set_use_warm_starts(use_warm_starts);
    lp_solver.lp_solve_method(lp_solve_method_id);
    lp_solver.load_problem(lp);

}

OperatorCountingHeuristic::~OperatorCountingHeuristic() {
}

void OperatorCountingHeuristic::print_statistics() const
{
    utils::g_log << "LP statistics for operator-counting heuristic:" << endl;
    lp_solver.print_statistics();
}

int OperatorCountingHeuristic::compute_heuristic(const State &ancestor_state) {
    State state = convert_ancestor_state(ancestor_state);
    assert(!lp_solver.has_temporary_constraints());
    for (const auto &generator : constraint_generators) {
        bool dead_end = generator->update_constraints(state, lp_solver);
        if (dead_end) {
            lp_solver.clear_temporary_constraints();
            return DEAD_END;
        }
    }
    int result;
    if (is_first) {
        lp_solver.solve_with_statistics();
        is_first = false;
    } else {
        lp_solver.solve();
    }
    if (lp_solver.has_optimal_solution()) {
        double epsilon = 0.01;
        double objective_value = lp_solver.get_objective_value();
        result = static_cast<int>(ceil(objective_value - epsilon));
    } else {
        result = DEAD_END;
    }
    lp_solver.clear_temporary_constraints();
    return result;
}

class OperatorCountingHeuristicFeature : public plugins::TypedFeature<Evaluator, OperatorCountingHeuristic> {
public:
    OperatorCountingHeuristicFeature() : TypedFeature("operatorcounting") {
        document_title("Operator-counting heuristic");
        document_synopsis(
            "An operator-counting heuristic computes a linear program (LP) in each "
            "state. The LP has one variable Count_o for each operator o that "
            "represents how often the operator is used in a plan. Operator-"
            "counting constraints are linear constraints over these varaibles that "
            "are guaranteed to have a solution with Count_o = occurrences(o, pi) "
            "for every plan pi. Minimizing the total cost of operators subject to "
            "some operator-counting constraints is an admissible heuristic. "
            "For details, see" + utils::format_conference_reference(
                {"Florian Pommerening", "Gabriele Roeger", "Malte Helmert",
                 "Blai Bonet"},
                "LP-based Heuristics for Cost-optimal Planning",
                "http://www.aaai.org/ocs/index.php/ICAPS/ICAPS14/paper/view/7892/8031",
                "Proceedings of the Twenty-Fourth International Conference"
                " on Automated Planning and Scheduling (ICAPS 2014)",
                "226-234",
                "AAAI Press",
                "2014"));

        add_list_option<shared_ptr<ConstraintGenerator>>(
            "constraint_generators",
            "methods that generate constraints over operator-counting variables");
        add_option<bool>(
            "use_integer_operator_counts",
            "restrict operator-counting variables to integer values. Computing the "
            "heuristic with integer variables can produce higher values but "
            "requires solving a MIP instead of an LP which is generally more "
            "computationally expensive. Turning this option on can thus drastically "
            "increase the runtime.",
            "false");
        
        add_option<bool>(
            "use_presolve",
            "turn presolving on or off. Using presolve creates an overhead "
            "so turning presolve off might decrease runtime.",
            "true");

        add_option<bool>(
            "crossover",
            "when set to false, turns off crossover when using the barrier optimizer. "
            "Crossover finds a basis, which is costly, but required for a warm start. "
            "Because of the warm start, it is turned on by default.",
            "true");

        add_option<int>(
            "symmetry_breaking_level",
            "turn symmetry_breaking off or on. -1 is automatic, 0 is off, "
            "1 to 5 are the different levels of symmetry breaking on, where "
            "1 is minimal and 5 is very aggressive.",
            "-1");

        add_option<int>(
            "folding_level",
            "turn folding_level off or on. -1 is automatic, 0 is off, "
            "1 to 5 are the different levels of folding, where "
            "1 is a moderate and 5 is an extremely aggressive level of folding.",
            "-1");

        add_option<bool>(
            "save_presolved_lp",
            "turn saving of the presolved lp to the disk off or on. "
            "When turned on, the presolved LP will be saved to the disk. "
            "When turned on, only the initial state will be evaluated, "
            "and then the planner will exit. This setting is not intended "
            "to be used in competitions and actual planning, but for "
            "research and understanding the presolving step only.",
            "false"
        );

        add_option<bool>(
            "use_warm_starts",
            "turn warm starts off or on. When turned on, the solution "
            "vector of the last LP is saved, and the next LP will "
            "be solved in two steps, taking the last solution vector as "
            "starting point. First, the constraints that relax "
            "the problem are added, then the new LP is solved using "
            "the primal simplex algorithm. In the second step, the tightening "
            "constraints are added, and this new LP is solved using "
            "the dual simplex algorithm.",
            "false");

        add_option<int>(
            "lp_solve_method",
            "determine which method the LP solver should use to solve the "
            "LPs given when using standard CPLEX (not twophase).",
            "0"
        );

        lp::add_lp_solver_option_to_feature(*this);
        Heuristic::add_options_to_feature(*this);

        document_language_support("action costs", "supported");
        document_language_support(
            "conditional effects",
            "not supported (the heuristic supports them in theory, but none of "
            "the currently implemented constraint generators do)");
        document_language_support(
            "axioms",
            "not supported (the heuristic supports them in theory, but none of "
            "the currently implemented constraint generators do)");

        document_property("admissible", "yes");
        document_property(
            "consistent",
            "yes, if all constraint generators represent consistent heuristics");
        document_property("safe", "yes");
        // TODO: prefer operators that are non-zero in the solution.
        document_property("preferred operators", "no");
    }

    virtual shared_ptr<OperatorCountingHeuristic> create_component(const plugins::Options &options, const utils::Context &context) const override {
        plugins::verify_list_non_empty<shared_ptr<ConstraintGenerator>>(
            context, options, "constraint_generators");
        return make_shared<OperatorCountingHeuristic>(options);
    }
};

static plugins::FeaturePlugin<OperatorCountingHeuristicFeature> _plugin;
}
