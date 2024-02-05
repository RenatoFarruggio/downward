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
    #SUITE = project.SUITE_SATISFICING
    SUITE = ["depot:p01.pddl", "grid:prob01.pddl", "gripper:prob01.pddl"
             , "organic-synthesis-split-opt18-strips:p17.pddl"]
    ENV = project.LocalEnvironment(processes=2)

bound = 0
CONFIGS = [
    # timing with presolve vs timing without presolve
    (f"config_timing_with_presolve", ["--search", "astar(operatorcounting([state_equation_constraints()],use_presolve=True))"]),
    (f"config_timing_without_presolve", ["--search", "astar(operatorcounting([state_equation_constraints()],use_presolve=False))"]),
]
BUILD_OPTIONS = []
DRIVER_OPTIONS = ["--overall-time-limit", "5m"]
REV_NICKS = [
    ("symmetrybreaking", "presolvetiming"),
]
ATTRIBUTES = [
    "error",
    "run_dir",
    "search_start_time",
    "search_start_memory",
    "total_time",
    "h_values",
    "coverage",
    "expansions",
    "memory",
    "presolve_time",
    "actual_search_time",
    "search_time",
    "my_total_time",
    "lp_solve_time_sum",
    project.EVALUATIONS_PER_TIME,
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
    total_time = run.get("total_time")
    lp_solve_time = run.get("lp_solve_time")
    # TODO
    run["presolve_time_per_total_time"] = lp_solve_time / total_time
    return run

def add_total_time_minus_my_total_time(run):
    total_time = run.get("total_time")
    my_total_time = run.get("my_total_time")
    if total_time is not None and my_total_time is not None:
        run["total_time_minus_my_total_time"] = total_time-my_total_time
    return run

def add_total_time_minus_actual_search_time(run):
    total_time = run.get("total_time")
    actual_search_time = run.get("actual_search_time")
    if total_time is not None and actual_search_time is not None:
        run["total_time_minus_actual_search_time"] = total_time - actual_search_time
    return run

def remove_short_running_lps(run): # NOT APPLIED
    lp_solve_time_sum = run.get("lp_solve_time_sum")
    return lp_solve_time_sum > 1

project.add_absolute_report(
    exp,
    attributes=ATTRIBUTES + ["total_time_minus_my_total_time", "total_time_minus_actual_search_time"], 
    filter=[add_total_time_minus_my_total_time, add_total_time_minus_actual_search_time, project.add_evaluations_per_time]
)

#comparative_attributes_list = [
#    ["total_time", "my_total_time"],
#    ["search_time", "total_time"]
#]

#for comparative_attributes in comparative_attributes_list:
#    project.add_comparative_report(
#        exp, 
#        attributes=comparative_attributes,
#        algorithm_pairs=[("presolvetiming:config_timing_with_presolve", "presolvetiming:config_timing_without_presolve")],
#        name_postfix=f"{comparative_attributes[0]}-vs-{comparative_attributes[1]}"
#    )

project.add_scp_step(exp, "grid", "workspace/downward-projects/")

exp.run_steps()