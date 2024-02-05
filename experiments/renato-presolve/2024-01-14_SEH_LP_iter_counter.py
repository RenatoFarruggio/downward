#! /usr/bin/env python

import os
import shutil

import project


REPO = project.get_repo_base()
BENCHMARKS_DIR = os.environ["DOWNWARD_BENCHMARKS"]
SCP_LOGIN = "farren00@login-infai.scicore.unibas.ch"
REMOTE_REPOS_DIR = "/infai/farren00/workspace/"
# If REVISION_CACHE is None, the default "./data/revision-cache/" is used.
REVISION_CACHE = os.environ.get("DOWNWARD_REVISION_CACHE")
if project.REMOTE:
    print("Running on REMOTE.")
    SUITE = project.SUITE_OPTIMAL_STRIPS
    ENV = project.BaselSlurmEnvironment(email="renato.farruggio@unibas.ch")
else:
    print("Running LOCALLY.")
    SUITE = ["depot:p01.pddl", "grid:prob01.pddl", "gripper:prob01.pddl"
             , "miconic:s1-0.pddl"
             , "miconic:s2-0.pddl"
    ]
    ENV = project.LocalEnvironment(processes=2)

cost_partitioning_uniform_optimal = "optimal"

lpsolver="cplex_twophase"

lm_factory_for_disjunctive_landmarks = "lm_rhw(disjunctive_landmarks=true, verbosity=normal, use_orders=true, only_causal_landmarks=false)"
landmark_disjunctive_landmarks = f"landmark_cost_partitioning({lm_factory_for_disjunctive_landmarks}, pref=false, prog_goal=true, prog_gn=true, prog_r=true, verbosity=normal, transform=no_transform(), cache_estimates=true, cost_partitioning={cost_partitioning_uniform_optimal}, alm=true, lpsolver={lpsolver})"

lm_factory_for_optimal_cost_partitioning = "lm_hm(m=2, conjunctive_landmarks=true, verbosity=normal, use_orders=true)"
landmark_optimal_cost_partitioning = f"landmark_cost_partitioning({lm_factory_for_optimal_cost_partitioning}, pref=false, prog_goal=true, prog_gn=true, prog_r=true, verbosity=normal, transform=no_transform(), cache_estimates=true, cost_partitioning=optimal, alm=true, lpsolver={lpsolver})"

CONFIGS = [
    (f"config_{lpsolver}_state_equation_heuristic", ["--search", f"astar(operatorcounting([state_equation_constraints()],lpsolver={lpsolver}))"]),  # state equation heuristic
    #(f"config_{lpsolver}_delete_relaxation_heuristic", ["--search", f"astar(operatorcounting([delete_relaxation_constraints()],lpsolver={lpsolver}))"]),  # delete relaxation heuristic
    #(f"config_{lpsolver}_disjunctive_action_landmark_heuristic", ["--search", f"astar({landmark_disjunctive_landmarks})"]),  # disjunctive action landmark heuristic for abstractions 
    #(f"config_{lpsolver}_optimal_cost_partitioning_heuristic", ["--search", f"astar({landmark_optimal_cost_partitioning})"]),  # optimal cost partitioning 
    #(f"config_{lpsolver}_post_hoc_optimization_heuristic", ["--search", f"astar(operatorcounting([pho_constraints()],lpsolver={lpsolver}))"]),  # post hoc optimization
    #(f"config_{lpsolver}_initial_state_potential_heuristic_presolve_off", ["--search", f"astar(initial_state_potential(use_presolve=false, lpsolver={lpsolver}))"]),   # potential heuristic
]


BUILD_OPTIONS = []
DRIVER_OPTIONS = ["--overall-time-limit", "5m"]
REV_NICKS = [
    ("symmetrybreaking", "state_equation_heuristic_LP_iter_counter_first_step"),
]
ATTRIBUTES = [
    "error",
    "run_dir",
    #"search_start_time",
    #"search_start_memory",
    #"total_time",
    #"h_values",
    "coverage",
    "expansions",
    "evaluations",
    #"memory",
    #"presolve_time",
    #"search_time",
    #"lp_solve_time_sum",
    "lp_variables",
    "lp_constraints",
    "lp_nonzero_entries",
    "iterations_for_problem_0",
    "iterations_for_problem_1",
    project.EVALUATIONS_PER_TIME,
    project.SPARSITY,
]

exp = project.FastDownwardExperiment(environment=ENV, revision_cache=REVISION_CACHE)
for config_nick, config in CONFIGS:
    for rev, rev_nick in REV_NICKS:
        algo_name = f"{rev_nick}:{config_nick}" if rev_nick else config_nick
        exp.add_algorithm(
            algo_name,
            REPO,
            rev,
            config,
            build_options=BUILD_OPTIONS,
            driver_options=DRIVER_OPTIONS,
        )
exp.add_suite(BENCHMARKS_DIR, SUITE)

exp.add_parser(exp.EXITCODE_PARSER)
exp.add_parser(exp.TRANSLATOR_PARSER)
exp.add_parser(exp.SINGLE_SEARCH_PARSER)
exp.add_parser(project.DIR / "parser.py")
exp.add_parser(exp.PLANNER_PARSER)

exp.add_step("build", exp.build)
exp.add_step("start", exp.start_runs)
exp.add_fetcher(name="fetch")

def add_presolve_time_per_total_time(run):
    exit()
    total_time = run.get("total_time")
    lp_solve_time = run.get("lp_solve_time")
    # TODO
    run["presolve_time_per_total_time"] = lp_solve_time / total_time
    return run

def add_total_time_minus_search_time(run):
    total_time = run.get("total_time")
    actual_search_time = run.get("search_time")
    if total_time is not None and actual_search_time is not None:
        run["total_time_minus_search_time"] = total_time - actual_search_time
    return run

def add_warm_starts_plus_cold_starts_minus_evaluations(run):
    num_warm_starts = run.get("num_warm_starts")
    num_cold_starts = run.get("num_cold_starts")
    evaluations = run.get("evaluations")
    if num_warm_starts is not None and num_cold_starts is not None and evaluations is not None:
        run["warm_starts_plus_cold_starts_minus_evaluations"] = num_warm_starts + num_cold_starts - evaluations
    return run


def add_repairs_per_warm_starts(run):
    num_tried_possible_repairs = run.get("num_tried_possible_repairs")
    num_warm_starts = run.get("num_warm_starts")
    if num_tried_possible_repairs is not None and num_warm_starts is not None and num_warm_starts != 0:
        run["repairs_per_warm_starts"] = num_tried_possible_repairs / num_warm_starts
    return run

def remove_short_running_lps(run):
    total_time = run.get("total_time")
    if total_time is not None:
        return total_time > 1
    return run

def remove_failed_runs(run):
    error_type = run.get("error")
    if error_type is not None:
        return error_type == "success"
    return run

def remove_trivial_runs_initial(run):
    iterations_0 = run.get("iterations_for_problem_0")
    if iterations_0 is not None:
        return iterations_0 > 0
    return run

def remove_trivial_runs_after_first_step(run):
    iterations_1 = run.get("iterations_for_problem_1")
    if iterations_1 is not None:
        return iterations_1 > 0
    return run

project.add_absolute_report(
    exp,
    attributes=ATTRIBUTES, 
    #filter=[add_total_time_minus_search_time, project.add_evaluations_per_time], 
    #filter=[remove_short_running_lps, add_total_time_minus_search_time, project.add_evaluations_per_time], 
    filter=[remove_trivial_runs_initial,
            remove_trivial_runs_after_first_step,
            project.add_variables_per_constraint,
            project.add_sparsity,
            project.add_evaluations_per_time]
)


# Plotting
#attributes = ["total_time"]
#pairs = [
#    ("presolvetiming:config_initial_state_potential_heuristic_presolve_on",
#     "presolvetiming:config_initial_state_potential_heuristic_presolve_off"),
#]
#suffix = "-rel" if project.RELATIVE else ""
#for algo1, algo2 in pairs:
#    for attr in attributes:
#        exp.add_report(
#            project.ScatterPlotReport(
#                relative=project.RELATIVE,
#                get_category=None if project.TEX else lambda run1, run2: run1["domain"],
#                attributes=[attr],
#                filter_algorithm=[algo1, algo2],
#                format="tex" if project.TEX else "png",
#            ),
#            name=f"{exp.name}-{algo1}-vs-{algo2}-{attr}{suffix}",
#        )

project.add_scp_step(exp, "grid", "workspace/downward-projects/")

exp.add_parse_again_step()

exp.run_steps()