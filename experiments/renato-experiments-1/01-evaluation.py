#! /usr/bin/env python

from pathlib import Path

from lab.experiment import Experiment

import project


ATTRIBUTES = [
    "error",
    "run_dir",
    "planner_time",
    "initial_h_value",
    "coverage",
    "cost",
    "evaluations",
    "memory",
    project.EVALUATIONS_PER_TIME,
]

exp = Experiment()
exp.add_step(
    "remove-combined-properties", project.remove_file, Path(exp.eval_dir) / "properties"
)

#project.fetch_algorithm(exp, "2023-09-30-presolve", "01-cg", new_algo="cg")
project.fetch_algorithms(exp, "2023-09-30-presolve")

filters = [project.add_evaluations_per_time]

project.add_absolute_report(
    exp, attributes=ATTRIBUTES, filter=filters, name=f"{exp.name}"
)


exp.run_steps()
