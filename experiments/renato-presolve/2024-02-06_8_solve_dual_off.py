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
    DRIVER_OPTIONS = ["--overall-time-limit", "5m"]
else:
    print("Running LOCALLY.")
    SUITE = ["depot:p01.pddl", "grid:prob01.pddl", "gripper:prob01.pddl"
             , "miconic:s1-0.pddl"
             , "miconic:s2-0.pddl"
             , "miconic:s5-0.pddl"
             , "movie:prob01.pddl"
             #, "movie:prob02.pddl"
             #, "movie:prob03.pddl"
             #, "movie:prob04.pddl"
             #, "movie:prob05.pddl"
             #, "movie:prob06.pddl"
             #, "movie:prob07.pddl"
             #, "movie:prob08.pddl"
             #, "movie:prob09.pddl"
             #, "movie:prob10.pddl"
             #, "movie:prob11.pddl"
             #, "movie:prob12.pddl"
             #, "movie:prob13.pddl"
             #, "movie:prob14.pddl"
             #, "movie:prob15.pddl"
             #, "movie:prob16.pddl"
             #, "movie:prob17.pddl"
             #, "movie:prob18.pddl"
             #, "movie:prob19.pddl"
             #, "movie:prob20.pddl"
             #, "movie:prob21.pddl"
             #, "movie:prob22.pddl"
             #, "movie:prob23.pddl"
             #, "movie:prob24.pddl"
             #, "movie:prob25.pddl"
             #, "movie:prob26.pddl"
             #, "movie:prob27.pddl"
             #, "movie:prob28.pddl"
             #, "movie:prob29.pddl"
             #, "movie:prob30.pddl"
             , "satellite:p01-pfile1.pddl"
             , "satellite:p02-pfile2.pddl"
             #, "satellite:p03-pfile3.pddl"
             #, "satellite:p04-pfile4.pddl"
             #, "satellite:p05-pfile5.pddl"
             #, "satellite:p06-pfile6.pddl"
             #, "satellite:p07-pfile7.pddl"
             #, "satellite:p08-pfile8.pddl"
             #, "satellite:p09-pfile9.pddl"
             #, "satellite:p10-pfile10.pddl"
             #, "satellite:p11-pfile11.pddl"
             #, "satellite:p12-pfile12.pddl"
             #, "satellite:p13-pfile13.pddl"
             #, "satellite:p14-pfile14.pddl"
             #, "satellite:p15-pfile15.pddl"
             #, "satellite:p16-pfile16.pddl"
             #, "satellite:p17-pfile17.pddl"
             #, "satellite:p18-pfile18.pddl"
             #, "satellite:p19-pfile19.pddl"
             #, "satellite:p20-pfile20.pddl"
             #, "satellite:p21-HC-pfile1.pddl"
             #, "satellite:p22-HC-pfile2.pddl"
             #, "satellite:p23-HC-pfile3.pddl"
             #, "satellite:p24-HC-pfile4.pddl"
             #, "satellite:p25-HC-pfile5.pddl"
             #, "satellite:p26-HC-pfile6.pddl"
             #, "satellite:p27-HC-pfile7.pddl"
             #, "satellite:p28-HC-pfile8.pddl"
             #, "satellite:p29-HC-pfile9.pddl"
             #, "satellite:p30-HC-pfile10.pddl"
             #, "satellite:p31-HC-pfile11.pddl"
             #, "satellite:p32-HC-pfile12.pddl"
             #, "satellite:p33-HC-pfile13.pddl"
             #, "satellite:p34-HC-pfile14.pddl"
             #, "satellite:p35-HC-pfile15.pddl"
             , "satellite:p36-HC-pfile16.pddl"
             #, "miconic:s20-0.pddl"
             #, "airport:p22-airport4halfMUC-p3.pddl"
    ]
    ENV = project.LocalEnvironment(processes=8)
    DRIVER_OPTIONS = ["--overall-time-limit", "1m"]


CONFIGS = [
    (f"config_SEH_default", ["--search", "astar(operatorcounting([state_equation_constraints()]))"]),
    (f"config_SEH_solve_dual_off", ["--search", "astar(operatorcounting([state_equation_constraints()],solve_dual=-1))"]),
    
    (f"config_DEL_default", ["--search", "astar(operatorcounting([delete_relaxation_constraints()]))"]),
    (f"config_DEL_solve_dual_off", ["--search", "astar(operatorcounting([delete_relaxation_constraints()],solve_dual=-1))"]),
    
    (f"config_OCP_default", ["--search", "astar(landmark_cost_partitioning(lm_rhw(),cost_partitioning=optimal))"]),
    (f"config_OCP_solve_dual_off", ["--search", "astar(landmark_cost_partitioning(lm_rhw(),cost_partitioning=optimal,solve_dual=-1))"]),
    
    (f"config_PHO_default", ["--search", "astar(operatorcounting([pho_constraints()]))"]),
    (f"config_PHO_solve_dual_off", ["--search", "astar(operatorcounting([pho_constraints()],solve_dual=-1))"]),
    
    (f"config_IPOT_default", ["--search", "astar(initial_state_potential())"]),
    (f"config_IPOT_solve_dual_off", ["--search", "astar(initial_state_potential(solve_dual=-1))"]),
    
    (f"config_DPOT_default", ["--search", "astar(diverse_potentials())"]),
    (f"config_DPOT_solve_dual_off", ["--search", "astar(diverse_potentials(solve_dual=-1))"]),
]

BUILD_OPTIONS = []
REV_NICKS = [
    ("symmetrybreaking", "solve_dual_off_vs_default"),
]
ATTRIBUTES = ["error",
    "run_dir",
    #"search_start_time",
    #"search_start_memory",
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
    "lp_solve_ticks",
    "initial_iterations_phase_1",
    "initial_iterations_phase_2",
    "initial_iterations_total",
    "initial_iterations_barrier",
    
    "iterations_phase_1",
    "iterations_phase_2",
    "iterations_total",
    "iterations_barrier",

    "initial_lp_solve_ticks",
    "phase_1_zero_iterations_count",
    "phase_1_nonzero_iterations_count",
    "total_zero_iterations_count",
    "total_nonzero_iterations_count",
    "lp_count",
    
    project.TOTAL_TIME,
    project.EVALUATIONS_PER_TIME,
    project.SPARSITY,
    project.AVERAGE_ITERATIONS_AFTER_INITIAL,
    project.AVERAGE_ITERATIONS_AFTER_INITIAL_PER_INITIAL_ITERATIONS,
    project.INITIAL_ITERATIONS_PER_AVERAGE_ITERATIONS_AFTER_INITIAL,
    project.PRESOLVE_TIME_PER_LP_SOLVE_TIME_SUM,
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
            project.add_initial_iterations_per_average_iterations_after_initial,
            project.add_variables_per_constraint,
            project.add_sparsity,
            project.add_evaluations_per_time,
            project.add_presolve_time_per_lp_solve_time_sum,],
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

for i in range(int(len(CONFIGS)/2)):
    algo1 = CONFIGS[2*i][0]
    algo2 = CONFIGS[2*i+1][0]
    project.add_comparative_report(
        exp,
        algorithm_pairs = [(REV_NICKS[0][1] + ":" + algo1,
                            REV_NICKS[0][1] + ":" + algo2)],
        attributes=ATTRIBUTES,
        #attributes=["coverage",
        #            "lp_count",
        #            "iterations_phase_1",
        #            "iterations_phase_2",
        #            "iterations_total",
        #            "iterations_barrier",
        #            "lp_solve_ticks",
        #            "lp_solve_time_sum",
        #            "phase_1_zero_iterations_count",
        #            "total_time",
        #            "average_iterations_after_initial",
        #            project.AVERAGE_ITERATIONS_AFTER_INITIAL_PER_INITIAL_ITERATIONS],
        filter=[project.add_lp_count,
                project.add_average_iterations_after_initial,
                project.add_average_iterations_after_initial_per_initial_iterations,
                project.add_initial_iterations_per_average_iterations_after_initial,
                project.add_variables_per_constraint,
                project.add_presolve_time_per_lp_solve_time_sum,],
        name_postfix=f"{algo1[7:-8]}-comp-to-default"
    )

for i in range(int(len(CONFIGS)/2)):
    algo1 = CONFIGS[2*i][0]
    algo2 = CONFIGS[2*i+1][0]
    project.add_comparative_report(
        exp,
        algorithm_pairs=[(REV_NICKS[0][1] + ":" + algo1,
                            REV_NICKS[0][1] + ":" + algo2)],
        attributes=["coverage",
                    project.TOTAL_TIME,
                    project.AVERAGE_ITERATIONS_AFTER_INITIAL,
                    project.AVERAGE_ITERATIONS_AFTER_INITIAL_PER_INITIAL_ITERATIONS,
                    project.INITIAL_ITERATIONS_PER_AVERAGE_ITERATIONS_AFTER_INITIAL,
                    "presolve_time",
                    "lp_solve_time_sum",
                    "error"],
        format="tex",
        filter=[project.add_lp_count,
                project.add_average_iterations_after_initial,
                project.add_average_iterations_after_initial_per_initial_iterations,
                project.add_initial_iterations_per_average_iterations_after_initial,
                project.add_variables_per_constraint,
                project.add_presolve_time_per_lp_solve_time_sum,],
        name_postfix=f"{algo1[7:-8]}-comp-to-default-tex",
    )


project.add_scp_step(exp, "grid", "workspace/downward-projects/")

exp.add_parse_again_step()

exp.run_steps()
