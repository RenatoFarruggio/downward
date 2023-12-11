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

    # "my_total_time"
    # Total time: 0.00657609s
    parser.add_pattern(
        "my_total_time",
        r"\] Total time: (.+)s",
        type=float,
    )
    parser.add_repeated_pattern(
        "lp_solve_time",
        r"LP solve time: (.+)s",
        type=float,
    )
    parser.add_pattern(
        "symmetry_breaking_level",
        r"Symmetry breaking is set to (.+)",
        type=int,
    )
    parser.parse()


if __name__ == "__main__":
    main()
