#! /usr/bin/env python

import logging
import re

from lab.parser import Parser


class CommonParser(Parser):
    def add_repeated_pattern(
        self, name, regex, file="run.log", required=False, type=int
    ):
        def find_all_occurences(content, props):
            matches = re.findall(regex, content)
            if required and not matches:
                logging.error(f"Pattern {regex} not found in file {file}")
            props[name] = [type(m) for m in matches]

        self.add_function(find_all_occurences, file=file)

    def add_bottom_up_pattern(
        self, name, regex, file="run.log", required=False, type=int
    ):
        def search_from_bottom(content, props):
            reversed_content = "\n".join(reversed(content.splitlines()))
            match = re.search(regex, reversed_content)
            if required and not match:
                logging.error(f"Pattern {regex} not found in file {file}")
            if match:
                props[name] = type(match.group(1))

        self.add_function(search_from_bottom, file=file)


def main():
    parser = CommonParser()
    parser.add_bottom_up_pattern(
        "search_start_time",
        r"\[t=(.+)s, \d+ KB\] g=0, 1 evaluated, 0 expanded",
        type=float,
    )
    parser.add_bottom_up_pattern(
        "search_start_memory",
        r"\[t=.+s, (\d+) KB\] g=0, 1 evaluated, 0 expanded",
        type=int,
    )
    parser.add_pattern(
        "initial_h_value",
        r"f = (\d+) \[1 evaluated, 0 expanded, t=.+s, \d+ KB\]",
        type=int,
    )
    parser.add_repeated_pattern(
        "h_values",
        r"New best heuristic value for .+: (\d+)\n",
        type=int,
    )
    # "presolve time"
    # Presolve time = 0.00 sec.
    parser.add_pattern(
        "presolve_time",
        r"Presolve time = (.+) sec.",
        type=float,
    )

    parser.add_pattern(
        "presolve_ticks",
        r"Presolve time = .+ sec. \((.+) ticks\)",
        type=float,
    )

    # "Actual search time"
    # Actual search time: 0.00147077s
    parser.add_pattern(
        "actual_search_time",
        r"\] Actual search time: (.+)s",
        type=float,
    )

    # "search_time"
    # Search time: 0.00241889s
    parser.add_pattern(
        "search_time",
        r"\] Search time: (.+)s",
        type=float,
    )


    parser.add_pattern(
        "lp_solve_time_sum",
        r"\] LP solve time sum: (.+)s",
        type=float,
    )
    #def sum_up_lp_solve_time(content, props):
    #    matches = re.findall(r"LP solve time sum: (.+)s", content)
    #    props["lp_solve_time_sum"] = sum([float(t) for t in matches])
    #parser.add_function(sum_up_lp_solve_time)

    parser.add_pattern(
        "symmetry_breaking_level",
        r"Symmetry breaking is set to (.+)",
        type=int,
    )

    parser.add_pattern(
        "lps_detected",
        r"LP detected: (.+)",
        type=int,
    )

    parser.add_pattern(
        "mips_detected",
        r"MIP detected: (.+)",
        type=int,
    )

    parser.add_pattern(
        "lp_variables",
        r"LP variables: (.+)",
        type=int,
    )

    parser.add_pattern(
        "lp_constraints",
        r"LP constraints: (.+)",
        type=int,
    )

    parser.add_pattern(
        "lp_nonzero_entries",
        r"LP non-zero entries: (.+)",
        type=int,
    )

    parser.add_pattern(
        "lp_variables_of_reduced",
        r"Reduced LP has .+ rows, (.+) columns, and .+ nonzeros.",
        type=int,
    )

    parser.add_pattern(
        "lp_constraints_of_reduced",
        r"Reduced LP has (.+) rows, .+ columns, and .+ nonzeros.",
        type=int,
    )

    parser.add_pattern(
        "lp_nonzero_entries_of_reduced",
        r"Reduced LP has .+ rows, .+ columns, and (.+) nonzeros.",
        type=int,
    )

    parser.add_pattern(
        "mip_starts_used",
        r"^(\d+) of \d+ MIP starts provided solutions.",
        type=int,
    )

    parser.add_bottom_up_pattern(
        "num_warm_starts",
        r"\] Warm starts: (.+)",
        type=int,
    )

    parser.add_bottom_up_pattern(
        "num_cold_starts",
        r"\] Cold starts: (.+)",
        type=int,
    )

    parser.add_bottom_up_pattern(
        "num_tried_possible_repairs",
        r"\] Attempted repairs: (.+)",
        type=int,
    )

    parser.add_pattern(
        "iterations_for_problem_0",
        r"Optimization in step 0 took (.+) iterations",
        type=int,
    )
    parser.add_pattern(
        "iterations_for_problem_1",
        r"Optimization in step 1 took (.+) iterations",
        type=int,
    )

    parser.add_pattern(
        "expansions_until_last_jump",
        r"\] Expanded until last jump: (.+) state(s).",
        type=int,
    )

    parser.add_bottom_up_pattern(
        "iterations_phase_1",
        r"\] LP solve phase 1 iterations: (.+)",
        type=int,
    )

    parser.add_bottom_up_pattern(
        "iterations_phase_2",
        r"\] LP solve phase 2 iterations: (.+)",
        type=int,
    )

    parser.add_bottom_up_pattern(
        "iterations_total",
        r"\] LP solve iterations total: (.+)",
        type=int,
    )

    parser.add_bottom_up_pattern(
        "lp_solve_ticks",
        r"\] LP solve ticks: (.+)",
        type=float,
    )

    parser.add_bottom_up_pattern(
        "initial_iterations_phase_1",
        r"\] First LP solve phase 1 iterations: (.+)",
        type=int,
    )

    parser.add_bottom_up_pattern(
        "initial_iterations_phase_2",
        r"\] First LP solve phase 2 iterations: (.+)",
        type=int,
    )

    parser.add_bottom_up_pattern(
        "initial_iterations_total",
        r"\] First LP solve iterations total: (.+)",
        type=int,
    )

    parser.add_bottom_up_pattern(
        "initial_lp_solve_ticks",
        r"\] First LP solve ticks: (.+)",
        type=float,
    )


    parser.add_bottom_up_pattern(
        "phase_1_zero_iterations_count",
        r"LP solve phase 1 zero iterations count: (.+)",
        type=int,
    )

    parser.add_bottom_up_pattern(
        "phase_1_nonzero_iterations_count",
        r"LP solve phase 1 nonzero iterations count: (.+)",
        type=int,
    )

    parser.add_bottom_up_pattern(
        "total_zero_iterations_count",
        r"LP solve total zero iterations count: (.+)",
        type=int,
    )

    parser.add_bottom_up_pattern(
        "total_nonzero_iterations_count",
        r"LP solve total nonzero iterations count: (.+)",
        type=int,
    )

    parser.parse()


if __name__ == "__main__":
    main()
