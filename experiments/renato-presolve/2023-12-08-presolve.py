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
    SUITE = ["depot:p01.pddl", "grid:prob01.pddl", "gripper:prob01.pddl"]
             #, "organic-synthesis-split-opt18-strips:p17.pddl"]
    ENV = project.LocalEnvironment(processes=2)

bound = 0
CONFIGS = [
    # save_presolved_lp=true
    # use_presolve=false
    (f"config_bound_{bound}_with_auto_symmetry_breaking", ["--search", f"astar(operatorcounting([state_equation_constraints()],symmetry_breaking_level=-1),bound={bound})"]),
    (f"config_bound_{bound}_with_no_symmetry_breaking", ["--search", f"astar(operatorcounting([state_equation_constraints()],symmetry_breaking_level=0),bound={bound})"]),
    (f"config_bound_{bound}_with_symmetry_breaking_level_1", ["--search", f"astar(operatorcounting([state_equation_constraints()],symmetry_breaking_level=1),bound={bound})"]),
    (f"config_bound_{bound}_with_symmetry_breaking_level_2", ["--search", f"astar(operatorcounting([state_equation_constraints()],symmetry_breaking_level=2),bound={bound})"]),
    (f"config_bound_{bound}_with_symmetry_breaking_level_3", ["--search", f"astar(operatorcounting([state_equation_constraints()],symmetry_breaking_level=3),bound={bound})"]),
    (f"config_bound_{bound}_with_symmetry_breaking_level_4", ["--search", f"astar(operatorcounting([state_equation_constraints()],symmetry_breaking_level=4),bound={bound})"]),
    (f"config_bound_{bound}_with_symmetry_breaking_level_5", ["--search", f"astar(operatorcounting([state_equation_constraints()],symmetry_breaking_level=5),bound={bound})"]),
]
BUILD_OPTIONS = []
DRIVER_OPTIONS = ["--overall-time-limit", "5m"]
REV_NICKS = [
    ("symmetrybreaking", "nonzerocounting"),
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
    "lp_solve_time",
    "symmetry_breaking_level",
    project.EVALUATIONS_PER_TIME,
]

ATTRIBUTES_COMPARATIVE = [
    "lp_solve_time",
    "total_time",
]

ALGORITHM_PAIRS = [
    (
        "nonzerocounting:config_bound_0_with_auto_symmetry_breaking",
        "nonzerocounting:config_bound_0_with_auto_symmetry_breaking",
        "nonzerocounting:config_bound_0_with_symmetry_breaking_level_5"
    )
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
    run["presolve_time_per_total_time"] = lp_solve_time / total_time
    return run

project.add_absolute_report(
    exp, attributes=ATTRIBUTES + ["presolve_time_per_total_time"], filter=[add_presolve_time_per_total_time, project.add_evaluations_per_time]
)

project.add_comparative_report(
    exp, 
    attributes=ATTRIBUTES_COMPARATIVE,
    algorithm_pairs=ALGORITHM_PAIRS
)

project.add_scp_step(exp, "grid", "workspace/downward-projects/")

exp.run_steps()
