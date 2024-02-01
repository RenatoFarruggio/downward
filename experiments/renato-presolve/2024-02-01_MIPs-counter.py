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
             , "elevators-opt08-strips:p04.pddl"    # both finish, and take the same time
             , "elevators-opt08-strips:p05.pddl"    # both finish, and take the same time
             , "elevators-opt08-strips:p06.pddl"    # both finish, and take the same time
             #, "woodworking-opt08-strips:p04.pddl"  # ON runs out of memory
             #, "woodworking-opt08-strips:p13.pddl"  # ON runs out of time
             #, "woodworking-opt08-strips:p14.pddl"  # ON runs out of memory
             #, "visitall-opt14-strips:p-1-7.pddl"  # OFF runs out of memory
    ]
    ENV = project.LocalEnvironment(processes=2)


CONFIGS = [
    (f"config_SEH_auto", ["--search", "astar(operatorcounting([state_equation_constraints()]),bound=0)"]),
    
    (f"config_DEL_auto", ["--search", "astar(operatorcounting([delete_relaxation_constraints()]),bound=0)"]),
    
    (f"config_OCP_auto", ["--search", "astar(landmark_cost_partitioning(lm_rhw(),cost_partitioning=optimal),bound=0)"]),
    
    (f"config_PHO_auto", ["--search", "astar(operatorcounting([pho_constraints()]),bound=0)"]),
    
    (f"config_IPOT_auto", ["--search", "astar(initial_state_potential(use_presolve=false),bound=0)"]),
    
    (f"config_DPOT_auto", ["--search", "astar(diverse_potentials(use_presolve=false),bound=0)"]),

]


BUILD_OPTIONS = []
DRIVER_OPTIONS = ["--overall-time-limit", "5m"]
REV_NICKS = [
    #("symmetrybreaking", "presolvetiming"),
    ("symmetrybreaking", "symmetrybreaking"),
]
ATTRIBUTES = [
    "error",
    "run_dir",
    #"search_start_time",
    #"search_start_memory",
    #"total_time",
    #"h_values",
    #"coverage",
    #"expansions",
    #"memory",
    #"presolve_time",
    #"search_time",
    #"lp_solve_time_sum",
    "lps_detected",
    "mips_detected",
    #"warm_start_counter",
    #"cold_start_counter",
    #project.EVALUATIONS_PER_TIME,
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

project.add_absolute_report(
    exp,
    attributes=ATTRIBUTES, 
    #filter=[add_total_time_minus_search_time, project.add_evaluations_per_time], 
    #filter=[remove_short_running_lps, add_total_time_minus_search_time, project.add_evaluations_per_time], 
    #filter=[project.add_evaluations_per_time]
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

# Comparative report
comparative_attributes_list = [
    ["coverage"],
    ["error"]
]
for comparative_attributes in comparative_attributes_list:
    project.add_comparative_report(
        exp, 
        attributes=comparative_attributes,
        algorithm_pairs=[("presolvetiming:config_initial_state_potential_heuristic_presolve_on",
                          "presolvetiming:config_initial_state_potential_heuristic_presolve_off")],
        name_postfix=f"{comparative_attributes[0]}",
        filter=[remove_failed_runs, project.add_evaluations_per_time]
    )


project.add_scp_step(exp, "grid", "workspace/downward-projects/")

exp.add_parse_again_step()

exp.run_steps()
