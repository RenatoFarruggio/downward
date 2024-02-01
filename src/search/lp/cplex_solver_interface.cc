#include "cplex_solver_interface.h"
#ifdef HAS_CPLEX

#include "lp_solver.h"

#include "../utils/language.h"
#include "../utils/logging.h"
#include "../utils/system.h"

#include <algorithm>
#include <cstring>
#include <numeric>

using namespace std;

namespace lp {
NO_RETURN
static void handle_cplex_error(CPXENVptr env, int error_code) {
    /*
      We checked for additional warnings when using the OSI framework, and
      either ignored them so they would not output anything or treated them a
      out-of-memory errors. In particular, we treated the following messages as
      running out of memory.

          "Insufficient memory for presolve."
          "Not enough memory for devex."
          "MIP starts not constructed because of out-of-memory status."

      It is unclear which channel these messages are sent trough, and so far we
      have not seen them in our tests. If they reoccur, we can probably capture
      them by adding a callback function to the correct message handler.
      https://www.ibm.com/docs/en/icos/22.1.1?topic=output-controlling-message-channels
      It could be that these were just warnings that would lead to the error
      below eventually. In that case there is nothing else we have to do.
    */
    if (error_code == CPXERR_NO_MEMORY) {
        utils::exit_with(utils::ExitCode::SEARCH_OUT_OF_MEMORY);
    }
    char buffer[CPXMESSAGEBUFSIZE];
    const char *error_string = CPXgeterrorstring(env, error_code, buffer);
    if (error_string) {
        cerr << error_string << endl;
    } else {
        cerr << "CPLEX error: Unknown error code " << error_code << endl;
    }
    utils::exit_with(utils::ExitCode::SEARCH_CRITICAL_ERROR);
}

/* Make a call to a CPLEX API function checking its return status. */
template<typename Func, typename ... Args>
static void CPX_CALL(Func cpxfunc, CPXENVptr env, Args && ... args) {
    int status = cpxfunc(env, forward<Args>(args) ...);
    if (status) {
        handle_cplex_error(env, status);
    }
}

static CPXLPptr createProblem(CPXENVptr env, const string &name) {
    int status = 0;
    CPXLPptr problem = CPXcreateprob(env, &status, name.c_str());
    if (status) {
        handle_cplex_error(env, status);
    }
    return problem;
}

static void freeProblem(CPXENVptr env, CPXLPptr *problem) {
    CPX_CALL(CPXfreeprob, env, problem);
}

static tuple<char, double, double> bounds_to_sense_rhs_range(double lb, double ub) {
    if (lb <= -CPX_INFBOUND && ub >= CPX_INFBOUND) {
        // CPLEX does not support <= or >= constraints without bounds.
        return {'R', -CPX_INFBOUND, 2 * CPX_INFBOUND};
    } else if (lb <= -CPX_INFBOUND) {
        return {'L', ub, 0};
    } else if (ub >= CPX_INFBOUND) {
        return {'G', lb, 0};
    } else if (lb == ub) {
        return {'E', lb, 0};
    } else {
        return {'R', lb, ub - lb};
    }
}

static int sense_to_cplex_sense(LPObjectiveSense sense) {
    if (sense == LPObjectiveSense::MINIMIZE) {
        return CPX_MIN;
    } else {
        return CPX_MAX;
    }
}

void CplexSolverInterface::CplexMatrix::assign_column_by_column(
    const named_vector::NamedVector<LPConstraint> &constraints,
    int num_cols) {
    coefficients.clear();
    indices.clear();
    starts.clear();
    counts.clear();

    // Set column starts and number of entries in each column.
    int num_rows = constraints.size();
    starts.resize(num_cols, 0);
    counts.resize(num_cols, 0);
    int num_nonzeros = 0;
    for (int row_index = 0; row_index < num_rows; ++row_index) {
        const vector<int> &vars = constraints[row_index].get_variables();
        num_nonzeros += vars.size();
        for (int var : vars) {
            ++counts[var];
        }
    }
    for (int var = 1; var < num_cols; ++var) {
        starts[var] = starts[var - 1] + counts[var - 1];
    }
    assert(num_nonzeros == starts[num_cols - 1] + counts[num_cols - 1]);

    /*
      Set non-zero coefficients in the correct places of 'coefficients'
      according to their column, and store the corresponding row in 'indices'.
     */
    indices.resize(num_nonzeros);
    coefficients.resize(num_nonzeros);
    for (int row_index = 0; row_index < num_rows; ++row_index) {
        const LPConstraint &constraint = constraints[row_index];
        const vector<int> &vars = constraint.get_variables();
        const vector<double> &coeffs = constraint.get_coefficients();
        int num_entries = vars.size();
        for (int i = 0; i < num_entries; ++i) {
            int var = vars[i];
            double coeff = coeffs[i];
            assert(starts[var] < num_nonzeros);
            indices[starts[var]] = row_index;
            coefficients[starts[var]] = coeff;
            /*
              We shift the start of the column to the right here and fix this
              after the loop. This way, we avoid an extra allocation and copy
              of the vector.
             */
            ++starts[var];
        }
    }
    /*
      We shifted col_start[i] to the right once for every entry. Now we have to
      shift it back to recover its original value.
     */
    for (int var = 0; var < num_cols; ++var) {
        starts[var] -= counts[var];
    }
    assert(static_cast<int>(starts.size()) == num_cols);
    assert(static_cast<int>(counts.size()) == num_cols);
    assert(indices.size() == coefficients.size());
}

void CplexSolverInterface::CplexMatrix::assign_row_by_row(
    const named_vector::NamedVector<LPConstraint> &constraints) {
    coefficients.clear();
    indices.clear();
    starts.clear();
    counts.clear();

    // Collect non-zero matrix entries.
    for (const LPConstraint &constraint : constraints) {
        const vector<int> &vars = constraint.get_variables();
        const vector<double> &coeffs = constraint.get_coefficients();
        assert(vars.size() == coeffs.size());
        starts.push_back(coefficients.size());
        indices.insert(indices.end(), vars.begin(), vars.end());
        coefficients.insert(coefficients.end(), coeffs.begin(), coeffs.end());
    }
    starts.push_back(coefficients.size());
    assert(static_cast<int>(starts.size()) == constraints.size() + 1);
    assert(indices.size() == coefficients.size());
}

void CplexSolverInterface::CplexColumnsInfo::assign(const named_vector::NamedVector<LPVariable> &variables) {
    lb.clear();
    ub.clear();
    type.clear();
    objective.clear();

    int num_cols = variables.size();
    lb.reserve(num_cols);
    ub.reserve(num_cols);
    type.reserve(num_cols);
    objective.reserve(num_cols);
    for (const LPVariable &var : variables) {
        lb.push_back(var.lower_bound);
        ub.push_back(var.upper_bound);
        if (var.is_integer) {
            type.push_back(CPX_INTEGER);
        } else {
            type.push_back(CPX_CONTINUOUS);
        }
        objective.push_back(var.objective_coefficient);
    }
    assert(static_cast<int>(lb.size()) == variables.size());
    assert(static_cast<int>(ub.size()) == variables.size());
    assert(static_cast<int>(type.size()) == variables.size());
    assert(static_cast<int>(objective.size()) == variables.size());
}

void CplexSolverInterface::CplexRowsInfo::assign(const named_vector::NamedVector<LPConstraint> &constraints, int offset, bool dense_range_values) {
    rhs.clear();
    sense.clear();
    range_values.clear();
    range_indices.clear();

    int num_rows = constraints.size();
    sense.resize(num_rows);
    rhs.resize(num_rows);
    if (dense_range_values) {
        range_values.resize(num_rows, 0);
    }
    for (int row_index = 0; row_index < num_rows; ++row_index) {
        const LPConstraint &constraint = constraints[row_index];
        double lb = constraint.get_lower_bound();
        double ub = constraint.get_upper_bound();
        const auto &[sense_value, rhs_value, range_value] = bounds_to_sense_rhs_range(lb, ub);
        sense[row_index] = sense_value;
        rhs[row_index] = rhs_value;
        if (sense_value == 'R') {
            if (dense_range_values) {
                range_values[row_index] = range_value;
            } else {
                range_values.push_back(range_value);
                range_indices.push_back(offset + row_index);
            }
        }
    }
    assert(static_cast<int>(rhs.size()) == constraints.size());
    assert(static_cast<int>(sense.size()) == constraints.size());
    assert(static_cast<int>(range_values.size()) <= constraints.size());
    assert((dense_range_values && (static_cast<int>(range_values.size()) == constraints.size()) && (range_indices.size() == 0)) ||
           (!dense_range_values && (range_values.size() == range_indices.size())));
}

CplexSolverInterface::CplexSolverInterface()
    : env(nullptr), problem(nullptr), 
      save_presolved_problem_to_file_and_exit(false), is_mip(false),
      num_permanent_constraints(0), 
      start_time(0), end_time(0),
      ticks_sum(0), iterations_sum_phase_1(0), iterations_sum_total(0),
      no_iterations_counter_phase_1(0), non_zero_iterations_counter_phase_1(0),
      no_iterations_counter_total(0), non_zero_iterations_counter_total(0),
      init_phase(true),
      cpx_lp_solve_method(CPXlpopt),
      num_unsatisfiable_constraints(0),
      num_unsatisfiable_temp_constraints(0) {
    int status = 0;
    env = CPXopenCPLEX(&status);
    if (status) {
        cerr << "Could not construct CPLEX interface (error_code: " << status
             << ")." << endl;
        utils::exit_with(utils::ExitCode::SEARCH_CRITICAL_ERROR);
    }
    CPX_CALL(CPXsetintparam, env, CPX_PARAM_THREADS, 1);
    // TODO handle output, catch oom
}

CplexSolverInterface::~CplexSolverInterface() {
    int status = CPXcloseCPLEX(&env);
    if (status) {
        cerr << "Failed to close CPLEX interface (error_code: " << status
             << ")." << endl;
        utils::exit_with(utils::ExitCode::SEARCH_CRITICAL_ERROR);
    }
}

bool CplexSolverInterface::is_trivially_unsolvable() const {
    return num_unsatisfiable_constraints + num_unsatisfiable_temp_constraints > 0;
}

void CplexSolverInterface::change_constraint_bounds(int index, double lb, double ub) {
    double current_lb = constraint_lower_bounds[index];
    double current_ub = constraint_upper_bounds[index];
    if (current_lb == lb && current_ub == ub) {
        return;
    }
    const auto &[sense, rhs, range] = bounds_to_sense_rhs_range(lb, ub);

    CPX_CALL(CPXchgsense, env, problem, 1, &index, &sense);
    CPX_CALL(CPXchgrhs, env, problem, 1, &index, &rhs);
    CPX_CALL(CPXchgrngval, env, problem, 1, &index, &range);

    if (current_lb > current_ub && lb <= ub) {
        if (index < num_permanent_constraints) {
            --num_unsatisfiable_constraints;
        } else {
            --num_unsatisfiable_temp_constraints;
        }
    } else if (current_lb <= current_ub && lb > ub) {
        if (index < num_permanent_constraints) {
            ++num_unsatisfiable_constraints;
        } else {
            ++num_unsatisfiable_temp_constraints;
        }
    }
    constraint_lower_bounds[index] = lb;
    constraint_upper_bounds[index] = ub;
}

void CplexSolverInterface::load_problem(const LinearProgram &lp) {
    if (problem) {
        freeProblem(env, &problem);
    }
    problem = createProblem(env, "");

    const named_vector::NamedVector<LPVariable> &variables = lp.get_variables();
    is_mip = any_of(variables.begin(), variables.end(), [](const LPVariable &v) {
                        return v.is_integer;
                    });

    const named_vector::NamedVector<LPConstraint> &constraints = lp.get_constraints();
    num_permanent_constraints = constraints.size();
    num_unsatisfiable_constraints = 0;
    for (const LPConstraint &constraint : constraints) {
        if (constraint.get_lower_bound() > constraint.get_upper_bound()) {
            ++num_unsatisfiable_constraints;
        }
    }

    matrix.assign_column_by_column(constraints, variables.size());
    columns.assign(variables);
    rows.assign(constraints);
    CPX_CALL(CPXcopylp, env, problem, variables.size(), constraints.size(),
             sense_to_cplex_sense(lp.get_sense()),
             columns.get_objective(),
             rows.get_rhs(),
             rows.get_sense(),
             matrix.get_starts(),
             matrix.get_counts(),
             matrix.get_indices(),
             matrix.get_coefficients(),
             columns.get_lb(),
             columns.get_ub(),
             rows.get_range_values());

    if (is_mip) {
        CPX_CALL(CPXcopyctype, env, problem, columns.get_type());
    } else {
        assert(CPXgetprobtype(env, problem) == CPXPROB_LP);
    }

    constraint_lower_bounds.clear();
    constraint_upper_bounds.clear();
    for (const LPConstraint &constraint : constraints) {
        constraint_lower_bounds.push_back(constraint.get_lower_bound());
        constraint_upper_bounds.push_back(constraint.get_upper_bound());
    }

    // Optionally set names.
    if (!lp.get_objective_name().empty()) {
        CPX_CALL(CPXchgprobname, env, problem, lp.get_objective_name().c_str());
    }
    if (variables.has_names()) {
        CplexNameData col_names(variables);
        CPX_CALL(CPXchgcolname, env, problem,
                 col_names.size(),
                 col_names.get_indices(),
                 col_names.get_names());
    }
    if (constraints.has_names()) {
        CplexNameData row_names(constraints);
        CPX_CALL(CPXchgrowname, env, problem,
                 row_names.size(),
                 row_names.get_indices(),
                 row_names.get_names());
    }
}

void CplexSolverInterface::add_temporary_constraints(
    const named_vector::NamedVector<LPConstraint> &constraints) {
    for (const LPConstraint &constraint : constraints) {
        if (constraint.get_lower_bound() > constraint.get_upper_bound()) {
            ++num_unsatisfiable_temp_constraints;
        }
    }

    matrix.assign_row_by_row(constraints);
    rows.assign(constraints, get_num_constraints(), false);
    CplexNameData row_names(constraints);
    // CPXaddrows can add new variables as well, but we do not want any.
    static const int num_extra_columns = 0;
    char **extra_column_names = nullptr;
    CPX_CALL(CPXaddrows, env, problem, num_extra_columns,
             constraints.size(),
             matrix.get_num_nonzeros(),
             rows.get_rhs(),
             rows.get_sense(),
             matrix.get_starts(),
             matrix.get_indices(),
             matrix.get_coefficients(),
             extra_column_names,
             row_names.get_names());

    /*
      If there are any ranged rows, we have to set up their ranges with a
      separate call.
    */
    if (rows.get_num_ranged_rows() > 0) {
        CPX_CALL(CPXchgrngval, env, problem,
                 rows.get_num_ranged_rows(),
                 rows.get_range_indices(),
                 rows.get_range_values());
    }

    for (const LPConstraint &constraint : constraints) {
        constraint_lower_bounds.push_back(constraint.get_lower_bound());
        constraint_upper_bounds.push_back(constraint.get_upper_bound());
    }
}

void CplexSolverInterface::clear_temporary_constraints() {
    int start = num_permanent_constraints;
    int end = get_num_constraints() - 1;
    if (start <= end) {
        CPX_CALL(CPXdelrows, env, problem, start, end);
        num_unsatisfiable_temp_constraints = 0;

        constraint_lower_bounds.resize(num_permanent_constraints);
        constraint_upper_bounds.resize(num_permanent_constraints);
    }
}

double CplexSolverInterface::get_infinity() const {
    return CPX_INFBOUND;
}

void CplexSolverInterface::set_objective_coefficients(const vector<double> &coefficients) {
    objective_indices.clear();
    objective_indices.resize(coefficients.size());
    iota(objective_indices.begin(), objective_indices.end(), 0);
    CPX_CALL(CPXchgobj, env, problem, coefficients.size(), objective_indices.data(), coefficients.data());
}

void CplexSolverInterface::set_objective_coefficient(int index, double coefficient) {
    CPX_CALL(CPXchgobj, env, problem, 1, &index, &coefficient);
}

void CplexSolverInterface::set_constraint_lower_bound(int index, double bound) {
    change_constraint_bounds(index, bound, constraint_upper_bounds[index]);
}

void CplexSolverInterface::set_constraint_upper_bound(int index, double bound) {
    change_constraint_bounds(index, constraint_lower_bounds[index], bound);
}

void CplexSolverInterface::set_variable_lower_bound(int index, double bound) {
    static const char bound_type = 'L';
    CPX_CALL(CPXchgbds, env, problem, 1, &index, &bound_type, &bound);
}

void CplexSolverInterface::set_variable_upper_bound(int index, double bound) {
    static const char bound_type = 'U';
    CPX_CALL(CPXchgbds, env, problem, 1, &index, &bound_type, &bound);
}

void CplexSolverInterface::set_mip_gap(double gap) {
    CPX_CALL(CPXsetdblparam, env, CPXPARAM_MIP_Tolerances_MIPGap, gap);
}

void CplexSolverInterface::solve() {
//    static int counter = 0;
//    cout << "Counter at " << counter << endl;
//    
//    if (counter < 2) {
//        cout << "Saving current problem as p" << counter << ".lp" << endl;
//        write_lp("p" + to_string(counter) + ".lp");
//    }

    if (is_trivially_unsolvable()) {
        return;
    } else if (is_mip) {
        CPX_CALL(CPXmipopt, env, problem);
    } else {
        if (init_phase) {
            CPX_CALL(CPXgetdettime, env, &start_time);
            CPX_CALL(CPXlpopt, env, problem);
            CPX_CALL(CPXgetdettime, env, &end_time);

            int iterations_phase_1 = CPXgetphase1cnt(env, problem);
            int iterations_total = CPXgetitcnt(env, problem);
            double delta_time = (end_time - start_time);

            utils::g_log << "First LP solve phase 1 iterations: " << iterations_phase_1 << endl;
            utils::g_log << "First LP solve phase 2 iterations: " << iterations_total - iterations_phase_1 << endl;
            utils::g_log << "First LP solve iterations total: " << iterations_total << endl;
            utils::g_log << "First LP solve ticks: " << delta_time << endl;
            init_phase = false;
            return;
        }

        CPX_CALL(CPXgetdettime, env, &start_time);
        CPX_CALL(cpx_lp_solve_method, env, problem);
        CPX_CALL(CPXgetdettime, env, &end_time);

        int iterations_phase_1 = CPXgetphase1cnt(env, problem);
        int iterations_total = CPXgetitcnt(env, problem);
        double delta_time = (end_time - start_time);


        iterations_sum_phase_1 += iterations_phase_1;
        iterations_sum_total += iterations_total;
        ticks_sum += delta_time;

        if (iterations_phase_1 == 0) {
            no_iterations_counter_phase_1++;
        } else {
            non_zero_iterations_counter_phase_1++;
        } 

        if (iterations_total == 0) {
            no_iterations_counter_total++;
        } else {
            non_zero_iterations_counter_total++;
        }
        
//        int method_used = CPXgetmethod(env, problem);
//        if (method_used == 1) {
//            cout << "Using Primal Simplex" << endl;
//        } else if (method_used == 2) {
//            cout << "Using Dual Simplex" << endl;
//        } else {
//            cout << "Using method " << method_used << endl;
//        }
//        cout << "Simplex in step " << counter << " took " << CPXgetphase1cnt(env, problem) << " iterations in phase 1" << endl;
//        cout << "Simplex in step " << counter << " took " << CPXgetitcnt(env, problem) << " iterations in total" << endl;
//        cout << "Simplex in step " << counter << " took " << end_time - start_time << " ticks" << endl;
//        cout << endl;
    }
//    counter++;
}

void CplexSolverInterface::solve_with_statistics() {
    print_statistics();
    // For log file: CPXsetlogfilename
    CPX_CALL(CPXsetintparam, env, CPXPARAM_ScreenOutput, CPX_ON);
    solve();
    CPX_CALL(CPXsetintparam, env, CPXPARAM_ScreenOutput, CPX_OFF);
}

void CplexSolverInterface::write_lp(const string &filename) const {
    if (is_trivially_unsolvable()) {
        cerr << "The LP has trivially unsatisfiable constraints that are not "
             << "accurately represented in CPLEX. Writing it to a file would "
             << "misrepresent the LP." << endl;
        utils::exit_with(utils::ExitCode::SEARCH_CRITICAL_ERROR);
    }
    // By not passing in a filetype, we let CPLEX infer it from the filename.
    static const char *filetype = nullptr;
    CPX_CALL(CPXwriteprob, env, problem, filename.c_str(), filetype);
}

void CplexSolverInterface::print_failure_analysis() const {
    if (is_trivially_unsolvable()) {
        cout << "LP/MIP is infeasible because of a trivially unsatisfiable "
             << "constraint" << endl;
        return;
    }
    int status = CPXgetstat(env, problem);
    switch (status) {
    case CPX_STAT_OPTIMAL:
        cout << "LP has an optimal solution." << endl;
        break;
    case CPXMIP_OPTIMAL:
        cout << "MIP has an optimal solution." << endl;
        break;
    case CPXMIP_OPTIMAL_TOL:
        cout << "MIP has an optimal solution within the tolerances of the "
             << "absolute/relative MIP gap." << endl;
        break;
    case CPX_STAT_OPTIMAL_INFEAS:
        cout << "LP has an optimal solution, but with infeasibilities after "
             << "unscaling." << endl;
        break;
    case CPX_STAT_UNBOUNDED:
        cout << "LP/MIP is unbounded" << endl;
        break;
    case CPX_STAT_INFEASIBLE:
        cout << "LP is infeasible" << endl;
        break;
    case CPXMIP_INFEASIBLE:
        cout << "MIP is infeasible" << endl;
        break;
    default:
        cout << "Unexpected status after solving LP/MIP: " << status << endl;
    }
}

bool CplexSolverInterface::is_infeasible() const {
    if (is_trivially_unsolvable()) {
        return true;
    }
    int status = CPXgetstat(env, problem);
    return status == CPX_STAT_INFEASIBLE || status == CPXMIP_INFEASIBLE;
}

bool CplexSolverInterface::is_unbounded() const {
    if (is_trivially_unsolvable()) {
        return false;
    }
    int status = CPXgetstat(env, problem);
    return status == CPX_STAT_UNBOUNDED;
}

bool CplexSolverInterface::has_optimal_solution() const {
    if (is_trivially_unsolvable()) {
        return false;
    }
    int status = CPXgetstat(env, problem);
    switch (status) {
    case CPX_STAT_OPTIMAL:
    case CPXMIP_OPTIMAL:
    /*
      The following status was returned in some cases, for example when
      computing h^+ for childsnack-opt14-strips/child-snack_pfile03-2.pddl.
      It means that the solution is optimal within the tolerances defined by
      the relative or absolute MIP gap.
      TODO: is it safe to treat this as an optimal solution (OSI did)?
    */
    case CPXMIP_OPTIMAL_TOL:
    /*
      The following status was returned in some cases, for example when
      computing diverse potential heuristics for Airport/p29.
      It means that there is an optimal solution of the LP but there are
      "infeasibilities after unscaling".
      TODO: is it safe to treat this as an optimal solution (OSI did)?
    */
    case CPX_STAT_OPTIMAL_INFEAS:
        return true;
    case CPX_STAT_UNBOUNDED:
    case CPX_STAT_INFEASIBLE:
    case CPXMIP_INFEASIBLE:
        return false;
    /*
      The following status was returned in some cases, for example when
      computing the State Equation Heuristic using the Barrier method for 
      organic-synthesis-split-opt18-strips/p03.pddl
      It means that a solution is available, but not proved optimal, due 
      to numeric difficulties during optimization. Note that the solution 
      may not be feasible. For example, it may be only dual feasible.
      Returning true yields a negative heuristic value, so we return false.
      TODO: maybe change something to prevent this error
    */
    case CPX_STAT_NUM_BEST:
        return false;
    /*
      The following status was returned in some cases, for example when
      computing the State Equation Heuristic using the Barrier method for 
      organic-synthesis-split-opt18-strips/p02.pddl
      It means that the solver stopped due to a limit on the dual objective.
      This is an indication that the primal is infeasible.
    */
    case CPX_STAT_ABORT_DUAL_OBJ_LIM:
        return false;
    default:
        cerr << "Unexpected status after solving LP/MIP: " << status << endl;
        utils::exit_with(utils::ExitCode::SEARCH_CRITICAL_ERROR);
    }
}

double CplexSolverInterface::get_objective_value() const {
    assert(has_optimal_solution());
    double value;
    CPX_CALL(CPXgetobjval, env, problem, &value);
    return value;
}

vector<double> CplexSolverInterface::extract_solution() const {
    assert(has_optimal_solution());
    int num_variables = get_num_variables();
    vector<double> solution(num_variables);
    CPX_CALL(CPXgetx, env, problem, solution.data(), 0, num_variables - 1);
    return solution;
}

int CplexSolverInterface::get_num_variables() const {
    return CPXgetnumcols(env, problem);
}

int CplexSolverInterface::get_num_constraints() const {
    return CPXgetnumrows(env, problem);
}

bool CplexSolverInterface::has_temporary_constraints() const {
    return num_permanent_constraints < get_num_constraints();
}

void CplexSolverInterface::print_statistics() const {
    utils::g_log << "LP variables: " << get_num_variables() << endl;
    utils::g_log << "LP constraints: " << get_num_constraints() << endl;
    utils::g_log << "LP non-zero entries: " << CPXgetnumnz(env, problem) << endl;

    if (!init_phase) {
        utils::g_log << "LP solve phase 1 iterations: " << iterations_sum_phase_1 << endl;
        utils::g_log << "LP solve phase 2 iterations: " << iterations_sum_total - iterations_sum_phase_1 << endl;
        utils::g_log << "LP solve iterations total: " << iterations_sum_total << endl;
        utils::g_log << "LP solve ticks: " << ticks_sum << endl;

        utils::g_log << "LP solve phase 1 zero iterations count: " << no_iterations_counter_phase_1 << endl;
        utils::g_log << "LP solve phase 1 nonzero iterations count: " << non_zero_iterations_counter_phase_1 << endl;
        utils::g_log << "LP solve total zero iterations count: " << no_iterations_counter_total << endl;
        utils::g_log << "LP solve total nonzero iterations count: " << non_zero_iterations_counter_total << endl;
    }
}

void CplexSolverInterface::set_use_presolve(bool use_presolve) {
    CPXINT value = use_presolve? CPX_ON : CPX_OFF;
    if (use_presolve) {
        cout << "switching presolve on" << endl;
    } else {
        cout << "switching presolve off" << endl;
    }
    CPX_CALL(CPXsetintparam, env, CPXPARAM_Preprocessing_Presolve, value);
}

void CplexSolverInterface::set_symmetry_breaking(int symmetry_breaking_level) {
    cout << "Symmetry breaking is set to " << symmetry_breaking_level << endl;
    CPXINT value = symmetry_breaking_level;
    CPX_CALL(CPXsetintparam, env, CPXPARAM_Preprocessing_Symmetry, value);
}

void CplexSolverInterface::set_crossover(bool use_crossover) {
    // See 
    // https://www.ibm.com/docs/en/icos/22.1.1?topic=parameters-solution-type-lp-qp
    if (use_crossover) {
        CPX_CALL(CPXsetintparam, env, CPXPARAM_SolutionType, CPX_AUTO_SOLN);
        cout << "Crossover enabled" << endl;
    } else {
        CPX_CALL(CPXsetintparam, env, CPXPARAM_SolutionType, CPX_NONBASIC_SOLN);
        cout << "Crossover disabled" << endl;
    }
}

void CplexSolverInterface::set_folding_level(int folding_level) {
    cout << "Folding level is set to " << folding_level << endl;
    CPXINT value = folding_level;
    CPX_CALL(CPXsetintparam, env, CPXPARAM_Preprocessing_Folding, value);
}

void CplexSolverInterface::set_save_presolved_lp(bool save_presolved_lp) {
    //CPXINT value = save_presolved_lp? CPX_ON: CPX_OFF;
    if (save_presolved_lp) {
        save_presolved_problem_to_file_and_exit = true;
        cout << "Saving presolved lp is turned on" << endl;
    } else {
        save_presolved_problem_to_file_and_exit = false;
        cout << "Saving presolved lp is turned off" << endl;
    }
}

void CplexSolverInterface::save_presolved_problem_to_file(std::string filename) {
    cout << "Saving presolved problem to file " << filename << " ..." << endl;
    double objOffset = 0.0;
    CPX_CALL(CPXpreslvwrite, env, problem, filename.c_str(), &objOffset);
    cout << "Offset is: " << objOffset << endl;
}

void CplexSolverInterface::set_use_warm_starts(bool use_warm_starts) {
    //CPXINT value = use_warm_starts? CPX_ON : CPX_OFF;
    if (use_warm_starts) {
        cout << "Using warm starts is turned on" << endl;
        cout << "[TODO] Warm starting LPs are not yet implemented! (this message comes from cplex_solver_interface.cc)" << endl;
    } else {
        cout << "Using warm starts is turned off" << endl;
    }
}

void CplexSolverInterface::lp_solve_method(int method_id) {
    // These ids are consistent (but not including network simplex and concurrent) with the CPLEX ids at:
    //  https://www.ibm.com/docs/en/icos/22.1.1?topic=parameters-algorithm-continuous-linear-problems


    cout << "Preferred method for solving LPs: ";
    if (method_id == 0) {
        cpx_lp_solve_method = CPXlpopt;
        cout << "None: Let CPLEX choose";
    } else if (method_id == 1) {
        cpx_lp_solve_method = CPXprimopt;
        cout << "Primal simplex";
    } else if (method_id == 2) {
        cpx_lp_solve_method = CPXdualopt;
        cout << "Dual simplex";
    } else if (method_id == 4) {
        cpx_lp_solve_method = CPXbaropt;
        cout << "Barrier optimizer";
    } else if (method_id == 5) {
        cpx_lp_solve_method = CPXsiftopt;
        cout << "Sifting";
    } else {
        cout << "WARNING: Admissible ids for lp_solve_method are 0 (auto), 1 (primal simplex), 2 (dual simplex), 4 (barrier), 5 (sifting)." << endl;
        cerr << "Unknown method name (" << method_id << ") provided for solving LPs!" << endl;
        exit(-1);
    }
    cout << endl;
}

void CplexSolverInterface::set_solve_dual(int solve_dual) {
    // The according documentation can be found at:
    //  https://www.ibm.com/docs/en/cofz/12.10.0?topic=performance-preprocessing


    cout << "Solving dual formulation is ";
    CPXINT value = solve_dual;
    CPX_CALL(CPXsetintparam, env, CPX_PARAM_PREDUAL, value);
    if (solve_dual == 0) {
        cout << "set to automatic" << endl;
    } else if (solve_dual == -1) {
        cout << "turned OFF" << endl;
    } else if (solve_dual == 1) {
        cout << "turned ON" << endl;
    } else {
        cerr << "provided with wrong input for solve_dual! Valid options are -1, 0, 1, but provided was " << solve_dual << endl;
    }
}

}
#endif
