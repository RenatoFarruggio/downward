#! /usr/bin/env python

import os
import shutil
import itertools

import project

from lab.reports import Attribute, arithmetic_mean

REPO = project.get_repo_base()
BENCHMARKS_DIR = os.environ["DOWNWARD_BENCHMARKS"]
SCP_LOGIN = "farren00@login-infai.scicore.unibas.ch"
REMOTE_REPOS_DIR = "/infai/farren00/workspace/"
# If REVISION_CACHE is None, the default "./data/revision-cache/" is used.
REVISION_CACHE = os.environ.get("DOWNWARD_REVISION_CACHE")
if project.REMOTE:
    print("Running on REMOTE.")
    SUITE = project.SUITE_OPTIMAL_STRIPS
    ENV = project.BaselSlurmEnvironment(email="renato.farruggio@unibas.ch",
                                        partition="infai_2")
else:
    print("Running LOCALLY.")
    SUITE = ["miconic:s1-0.pddl"
             #, "depot:p01.pddl"
             #, "grid:prob01.pddl", "gripper:prob01.pddl"
             #, "organic-synthesis-split-opt18-strips:p02.pddl"
             #, "elevators-opt08-strips:p04.pddl"    # both finish, and take the same time
             #, "elevators-opt08-strips:p05.pddl"    # both finish, and take the same time
             #, "elevators-opt08-strips:p06.pddl"    # both finish, and take the same time
             #, "woodworking-opt08-strips:p04.pddl"  # ON runs out of memory
             #, "woodworking-opt08-strips:p13.pddl"  # ON runs out of time
             #, "woodworking-opt08-strips:p14.pddl"  # ON runs out of memory
             #, "visitall-opt14-strips:p-1-7.pddl"  # OFF runs out of memory
    ]
    ENV = project.LocalEnvironment(processes=8)

#    (f"SEH", ["--search", "astar(operatorcounting([state_equation_constraints()]))"]),
#    (f"DEL", ["--search", "astar(operatorcounting([delete_relaxation_constraints()]))"]),
#    (f"OCP", ["--search", "astar(landmark_cost_partitioning(lm_rhw(),cost_partitioning=optimal))"]),
#    (f"PHO", ["--search", "astar(operatorcounting([pho_constraints()]))"]),
#    (f"IPOT", ["--search", "astar(initial_state_potential(use_presolve=true))"]),
#    (f"DPOT", ["--search", "astar(diverse_potentials(use_presolve=true))"]),
def heuristic_to_command(heuristic_str: str, settings: str):
    match(heuristic_str):
        case "SEH":
            return f"astar(operatorcounting([state_equation_constraints()]{settings}))"
        case "DEL":
            return f"astar(operatorcounting([delete_relaxation_constraints()]{settings}))"  
        case "OCP":
            return f"astar(landmark_cost_partitioning(lm_rhw(),cost_partitioning=optimal{settings}))"
        case "PHO":
            return f"astar(operatorcounting([pho_constraints()]{settings}))"
        case "IPOT":
            return f"astar(initial_state_potential({settings.strip(',')}))"
        case "DPOT":
            return f"astar(diverse_potentials({settings.strip(',')}))"
        case _:
            print(f"ERROR: Unknown heuristic: {heuristic_str}")

MAX_K = 1
CONFIGS = []
heuristics = ["SEH", 
              "DEL", 
              "OCP",
              "PHO", "IPOT", "DPOT"
              ]
optimizers = ["",  # Default: 0 
              #"lp_solve_method=1", 
              #"lp_solve_method=2", 
              #"lp_solve_method=4"
              ]

cplex_settings_presolve = ["",  # Default: true 
                           #"use_presolve=false"
                           ]
cplex_settings_solve_dual = ["",  # Default: 0
                             #"solve_dual=-1", 
                             #"solve_dual=1"
                             ]
cplex_settings_folding_level = ["",  # Default: -1
                                #"folding_level=0",
                                #"folding_level=1"
                                ]

cplex_settings_aggregator_application_limit = ["aggregator_application_limit=-1",  # Default: -1
                                               #"aggregator_application_limit=0",
                                               ]

cplex_preprocessing_settings = itertools.product(cplex_settings_presolve, 
                                                 cplex_settings_solve_dual, 
                                                 cplex_settings_folding_level, 
                                                 cplex_settings_aggregator_application_limit)
cplex_preprocessing_settings = [[a for a in setting if len(a)>0] for setting in cplex_preprocessing_settings]

def cplex_preprocessing_settings_to_string(cplex_preprocessing_setting) -> str:
    if len(cplex_preprocessing_setting) > 0:
        return ",".join(cplex_preprocessing_setting)
    else:
        return ""


#for cplex_preprocessing_setting in cplex_preprocessing_settings:
#    print(f"({len(cplex_preprocessing_setting)}) {cplex_preprocessing_settings_to_string(cplex_preprocessing_setting)}")


settings = [""]
for heuristic in heuristics:
    for optimizer in optimizers:
        for cplex_preprocessing_setting in cplex_preprocessing_settings:
            if len(cplex_preprocessing_setting) >= MAX_K + 1:
                continue

            setting = ","
            setting += "," + optimizer
            setting += "," + cplex_preprocessing_settings_to_string(cplex_preprocessing_setting)

            setting = setting.strip(",")
            if len(setting) > 0:
                setting = "," + setting

            setting_for_name = setting.replace('=', '_').replace(',', '_').strip("_")
            
            name = f"config_{heuristic}_{setting_for_name}".strip("_")
            CONFIGS.append((f"{name}", ["--search", f"{heuristic_to_command(heuristic, setting)}"]))


print(f"Setting up experiment with {len(CONFIGS)} configs with at most {MAX_K} non-default preprocessing settings")
#counter = 0
#for config in CONFIGS:
#    counter += 1
#    print(f"({counter}/{len(CONFIGS)}) {config}")
#    print()



BUILD_OPTIONS = []
DRIVER_OPTIONS = ["--overall-time-limit", "1m"]
REV_NICKS = [
    #("symmetrybreaking", "presolvetiming"),
    ("symmetrybreaking", "bruteforce_mini_short"),
]

DENSITY = Attribute(
    "density", min_wins=False, function=arithmetic_mean, digits=5
)

DENSITY_OF_REDUCED = Attribute(
    "density_of_reduced", min_wins=False, function=arithmetic_mean, digits=5
)

DENSITY_OF_REDUCED_PER_DENSITY = Attribute(
    "density_of_reduced_per_density", min_wins=True, function=arithmetic_mean, digits=2
)
ATTRIBUTES = [
    "error",
    "run_dir",
    #"search_start_time",
    #"search_start_memory",
    project.TOTAL_TIME,
    #"h_values",
    "coverage",
    #"expansions",
    #"memory",
    "presolve_time",
    "presolve_ticks",
    "initial_lp_solve_ticks",
    #"search_time",
    #"lp_solve_time_sum",
    #"lps_detected",
    #"mips_detected",
    #project.EVALUATIONS_PER_TIME,
    "lp_variables",
    "lp_constraints",
    "lp_nonzero_entries",
    "lp_variables_of_reduced",
    "lp_constraints_of_reduced",
    "lp_nonzero_entries_of_reduced",
    "solved_during_presolve",
    DENSITY,
    DENSITY_OF_REDUCED,
    DENSITY_OF_REDUCED_PER_DENSITY,
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


def add_density(run):
    lp_variables = run.get("lp_variables")
    lp_constraints = run.get("lp_constraints")
    lp_nonzero_entries = run.get("lp_nonzero_entries")
    if lp_variables is not None and lp_constraints is not None:
        size = lp_variables * lp_constraints
        if lp_nonzero_entries is not None and size > 0:
            run["density"] = lp_nonzero_entries / size
    
    if "density" not in run:
        run["density"] = None
    return run

def add_density_of_reduced(run):
    lp_variables_of_reduced = run.get("lp_variables_of_reduced")
    lp_constraints_of_reduced = run.get("lp_constraints_of_reduced")
    lp_nonzero_entries_of_reduced = run.get("lp_nonzero_entries_of_reduced")
    if lp_variables_of_reduced is not None and lp_constraints_of_reduced is not None:
        size_of_reduced = lp_variables_of_reduced * lp_constraints_of_reduced
        if lp_nonzero_entries_of_reduced is not None and size_of_reduced > 0:
            run["density_of_reduced"] = lp_nonzero_entries_of_reduced / size_of_reduced
    
    if "density_of_reduced" not in run:
        run["density_of_reduced"] = None
    return run

def add_density_of_reduced_per_density(run):
    density = run.get("density")
    density_of_reduced = run.get("density_of_reduced")
    if density is not None and density_of_reduced is not None and density > 0:
        run["density_of_reduced_per_density"] = density_of_reduced / density
    
    if "density_of_reduced_per_density" not in run:
        run["density_of_reduced_per_density"] = None
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

algorithms = [REV_NICKS[0][1] + ":" + config[0] for config in CONFIGS]

project.add_absolute_report(
    exp,
    attributes=ATTRIBUTES,
    filter=[add_density,
            add_density_of_reduced,
            add_density_of_reduced_per_density],
    #filter=[add_total_time_minus_search_time, project.add_evaluations_per_time], 
    #filter=[remove_short_running_lps, add_total_time_minus_search_time, project.add_evaluations_per_time], 
    #filter=[project.add_evaluations_per_time]
)

#for algorithm in algorithms:
#    project.add_absolute_report(
#        exp,
#        filter_algorithm=[algorithm],
#        attributes=["coverage",
#                    "error",
#                    project.TOTAL_TIME],
#        filter=[remove_failed_runs,
#                add_density,
#                add_density_of_reduced,
#                add_density_of_reduced_per_density],
#        name=f"{exp.name}-{algorithm}-abs-success-only"
#    )

#####

#project.add_absolute_report(
#    exp,
#    attributes=["presolve_time_per_total_time"],
#    format="tex",
#    name="tex-report-presolve_time_per_total_time-all",
#    filter_algorithm=algorithms
#)
#
#for algorithm in algorithms:
#    project.add_absolute_report(
#        exp,
#        attributes=["presolve_time_per_total_time"],
#        format="tex",
#        name=f"tex-report-presolve_time_per_total_time-{algorithm}",
#        filter_algorithm=[algorithm]
#    )


#project.add_absolute_report(
#    exp,
#    attributes=["variables_per_constraint_of_presolved",
#                "constraints_per_variable_of_presolved"],
#    format="tex",
#    name="tex-report-variables_per_constraint_of_presolved-all",
#    filter_algorithm=algorithms
#)
#
#for algorithm in algorithms:
#    project.add_absolute_report(
#        exp,
#        attributes=["variables_per_constraint_of_presolved"],
#        format="tex",
#        name=f"tex-report-variables_per_constraint_of_presolved-{algorithm}",
#        filter_algorithm=[algorithm]
#    )

#####################################################

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
#comparative_attributes_list = [
#    ["coverage"],
#    ["error"]
#]
#for comparative_attributes in comparative_attributes_list:
#    project.add_comparative_report(
#        exp, 
#        attributes=comparative_attributes,
#        algorithm_pairs=[("presolvetiming:config_initial_state_potential_heuristic_presolve_on",
#                          "presolvetiming:config_initial_state_potential_heuristic_presolve_off")],
#        name_postfix=f"{comparative_attributes[0]}",
#        filter=[remove_failed_runs, project.add_evaluations_per_time]
#    )


project.add_scp_step(exp, "grid", "workspace/downward-projects/")

exp.add_parse_again_step()

exp.run_steps()