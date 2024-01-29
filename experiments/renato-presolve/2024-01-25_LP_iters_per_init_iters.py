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
             , "miconic:s5-0.pddl"
             #, "miconic:s20-0.pddl"
             #, "airport:p22-airport4halfMUC-p3.pddl"
    ]
    ENV = project.LocalEnvironment(processes=2)


CONFIGS = [
    (f"config_state_equation_heuristic", ["--search", "astar(operatorcounting([state_equation_constraints()]))"]),
    (f"config_delete_relaxation_heuristic", ["--search", "astar(operatorcounting([delete_relaxation_constraints()]))"]),
    (f"config_optimal_cost_partitioning_heuristic", ["--search", "astar(landmark_cost_partitioning(lm_rhw(),cost_partitioning=optimal))"]),
    (f"config_pho_heuristic", ["--search", "astar(operatorcounting([pho_constraints()]))"]),
    (f"config_initial_potential_heuristic", ["--search", "astar(initial_state_potential(use_presolve=false))"]),
    (f"config_diverse_potential_heuristic", ["--search", "astar(diverse_potentials(use_presolve=false))"]),
]

BUILD_OPTIONS = []
DRIVER_OPTIONS = ["--overall-time-limit", "5m"]
REV_NICKS = [
    ("symmetrybreaking", "iterations_per_initial_iterations"),
]
ATTRIBUTES = ["error",
    "run_dir",
    #"search_start_time",
    #"search_start_memory",
    #"total_time",
    #"h_values",
    "coverage",
    "expansions",
    "evaluations",
    #"memory",
    "presolve_time",
    "presolve_ticks",
    "search_time",
    "lp_solve_time_sum",
    "lp_variables",
    "lp_constraints",
    "lp_nonzero_entries",
    "expansions_until_last_jump",
    "iterations_phase_1",
    "iterations_phase_2",
    "iterations_total",
    "lp_solve_ticks",
    "initial_iterations_phase_1",
    "initial_iterations_phase_2",
    "initial_iterations_total",
    "initial_lp_solve_ticks",
    "phase_1_zero_iterations_count",
    "phase_1_nonzero_iterations_count",
    "total_zero_iterations_count",
    "total_nonzero_iterations_count",
    "lp_count",
    "average_iterations_after_initial",
    project.EVALUATIONS_PER_TIME,
    project.SPARSITY,
    project.AVERAGE_ITERATIONS_AFTER_INITIAL_PER_INITIAL_ITERATIONS
    
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

def replace_None_presolve_time_with_zero(run):
    presolve_time = run.get("presolve_time")
    if presolve_time is None:
        run["presolve_time"] = 0.0
    return run

def replace_None_presolve_ticks_with_zero(run):
    presolve_ticks = run.get("presolve_ticks")
    if presolve_ticks is None:
        run["presolve_ticks"] = 0.0
    return run

project.add_absolute_report(
    exp,
    attributes=ATTRIBUTES,
    #filter=[add_total_time_minus_search_time, project.add_evaluations_per_time], 
    #filter=[remove_short_running_lps, add_total_time_minus_search_time, project.add_evaluations_per_time], 
    filter=[project.add_lp_count,
            project.add_average_iterations_after_initial,
            project.add_average_iterations_after_initial_per_initial_iterations,
            project.add_variables_per_constraint,
            project.add_sparsity,
            project.add_evaluations_per_time]
)



#for i in range(int(len(CONFIGS)/2)):
#    algo1 = CONFIGS[2*i][0]
#    algo2 = CONFIGS[2*i+1][0]
#
#    project.add_comparative_report(
#        exp,
#        algorithm_pairs = [(f"iter_counter:{algo1}", 
#                            f"iter_counter:{algo2}")],
#        attributes=["coverage", 
#                    "lp_solve_time_sum", 
#                    "total_time", 
#                    "search_time",
#                    "expansions_until_last_jump",
#                    "iterations_phase_1",
#                    "iterations_phase_2",
#                    "iterations_total",
#                    "lp_solve_ticks",
#                    "initial_iterations_phase_1",
#                    "initial_iterations_phase_2",
#                    "initial_iterations_total",
#                    "initial_lp_solve_ticks",    
#                    "phase_1_zero_iterations_count",
#                    "phase_1_nonzero_iterations_count",
#                    "total_zero_iterations_count",
#                    "total_nonzero_iterations_count"
#                    ],
#        name_postfix=f"comparison-{algo1[22:]}"
#    )

algorithms = [REV_NICKS[0][1] + ":" + config[0] for config in CONFIGS]

print()
for algorithm_name in algorithms:
    algorithm_name_short = algorithm_name.split(':')[1][7:-10]

    project.add_comparative_report(
        exp,
        algorithm_pairs = [(algorithm_name)],
        attributes=["average_iterations_after_initial_per_initial_iterations"],
        name_postfix=f"{algorithm_name_short}"
)
print()


# TODO: Add custom filter or use pandas

# # Plotting
# def domain_as_category(run1, run2):
#     # run2['domain'] has the same value, because we always
#     # compare two runs of the same problem.
#     return run1["domain"]
# 
# exp.add_report(
#     project.ScatterPlotReport(
#         relative=True,
#         get_category=domain_as_category,
#         attributes=["coverage"],
#         filter_algorithm=["presolve_timing:config_initial_potential_heuristic_presolve_ON", 
#                           "presolve_timing:config_initial_potential_heuristic_presolve_OFF"
#                           ],
#         format="png",
#     ),
#     name=f"{exp.name}-coverage-rel"
# )



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


#for config in CONFIGS:
#    algorithm_name = REV_NICKS[0][1] + ":" + config[0]
#    exp.add_report(
#        project.ScatterPlotReport(
#            relative=False,
#            attributes=["initial_iterations_total", "iterations_total"],
#            filter_algorithm=[algorithm_name],
#            title="Initial vs Total Iterations for Algorithm-1",
#            xlabel="Initial Iterations Total",
#            ylabel="Iterations Total",
#            format="png",
#            ),
#        name=f"{config[0][20:]}-iterations_per_initial",
#    )

project.add_scp_step(exp, "grid", "workspace/downward-projects/")

exp.add_parse_again_step()

exp.run_steps()
