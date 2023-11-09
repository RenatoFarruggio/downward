from pathlib import Path
import platform
import re
import subprocess
import sys

from downward.experiment import FastDownwardExperiment
from downward.reports.absolute import AbsoluteReport
from downward.reports.scatter import ScatterPlotReport
from downward.reports.taskwise import TaskwiseReport
from lab.environments import (
    BaselSlurmEnvironment,
    LocalEnvironment,
    TetralithEnvironment,
)
from lab.experiment import ARGPARSER
from lab.reports import Attribute, geometric_mean


# Silence import-unused messages. Experiment scripts may use these imports.
assert (
    BaselSlurmEnvironment
    and FastDownwardExperiment
    and LocalEnvironment
    and ScatterPlotReport
    and TaskwiseReport
    and TetralithEnvironment
)


DIR = Path(__file__).resolve().parent
SCRIPT = Path(sys.argv[0]).resolve()

NODE = platform.node()
# Cover both the Basel and Linköping clusters for simplicity.
REMOTE = NODE.endswith((".scicore.unibas.ch", ".cluster.bc2.ch")) or re.match(
    r"tetralith\d+\.nsc\.liu\.se|n\d+", NODE
)


def parse_args():
    ARGPARSER.add_argument("--tex", action="store_true", help="produce LaTeX output")
    ARGPARSER.add_argument(
        "--relative", action="store_true", help="make relative scatter plots"
    )
    return ARGPARSER.parse_args()


ARGS = parse_args()
TEX = ARGS.tex
RELATIVE = ARGS.relative

EVALUATIONS_PER_TIME = Attribute(
    "evaluations_per_time", min_wins=False, function=geometric_mean, digits=1
)

# Generated by "./suites.py satisficing" in aibasel/downward-benchmarks repo.
# fmt: off
SUITE_SATISFICING = [
    "agricola-sat18-strips", "airport", "assembly", "barman-sat11-strips",
    "barman-sat14-strips", "blocks", "caldera-sat18-adl",
    "caldera-split-sat18-adl", "cavediving-14-adl", "childsnack-sat14-strips",
    "citycar-sat14-adl", "data-network-sat18-strips", "depot", "driverlog",
    "elevators-sat08-strips", "elevators-sat11-strips", "flashfill-sat18-adl",
    "floortile-sat11-strips", "floortile-sat14-strips", "freecell",
    "ged-sat14-strips", "grid", "gripper", "hiking-sat14-strips",
    "logistics00", "logistics98", "maintenance-sat14-adl", "miconic",
    "miconic-fulladl", "miconic-simpleadl", "movie", "mprime", "mystery",
    "nomystery-sat11-strips", "nurikabe-sat18-adl", "openstacks",
    "openstacks-sat08-adl", "openstacks-sat08-strips",
    "openstacks-sat11-strips", "openstacks-sat14-strips", "openstacks-strips",
    "optical-telegraphs", "organic-synthesis-sat18-strips",
    "organic-synthesis-split-sat18-strips", "parcprinter-08-strips",
    "parcprinter-sat11-strips", "parking-sat11-strips", "parking-sat14-strips",
    "pathways", "pegsol-08-strips", "pegsol-sat11-strips", "philosophers",
    "pipesworld-notankage", "pipesworld-tankage", "psr-large", "psr-middle",
    "psr-small", "rovers", "satellite", "scanalyzer-08-strips",
    "scanalyzer-sat11-strips", "schedule", "settlers-sat18-adl",
    "snake-sat18-strips", "sokoban-sat08-strips", "sokoban-sat11-strips",
    "spider-sat18-strips", "storage", "termes-sat18-strips",
    "tetris-sat14-strips", "thoughtful-sat14-strips", "tidybot-sat11-strips",
    "tpp", "transport-sat08-strips", "transport-sat11-strips",
    "transport-sat14-strips", "trucks", "trucks-strips",
    "visitall-sat11-strips", "visitall-sat14-strips",
    "woodworking-sat08-strips", "woodworking-sat11-strips", "zenotravel",
]

SUITE_OPTIMAL_STRIPS = [
    "agricola-opt18-strips", "airport", "barman-opt11-strips",
    "barman-opt14-strips", "blocks", "childsnack-opt14-strips",
    "data-network-opt18-strips", "depot", "driverlog", "elevators-opt08-strips",
    "elevators-opt11-strips", "floortile-opt11-strips", "floortile-opt14-strips",
    "freecell", "ged-opt14-strips", "grid", "gripper", "hiking-opt14-strips",
    "logistics00", "logistics98", "miconic", "movie", "mprime", "mystery",
    "nomystery-opt11-strips", "openstacks-opt08-strips", "openstacks-opt11-strips",
    "openstacks-opt14-strips", "openstacks-strips", "organic-synthesis-opt18-strips",
    "organic-synthesis-split-opt18-strips", "parcprinter-08-strips",
    "parcprinter-opt11-strips", "parking-opt11-strips", "parking-opt14-strips",
    "pathways", "pegsol-08-strips", "pegsol-opt11-strips",
    "petri-net-alignment-opt18-strips", "pipesworld-notankage", "pipesworld-tankage",
    "psr-small", "rovers", "satellite", "scanalyzer-08-strips",
    "scanalyzer-opt11-strips", "snake-opt18-strips", "sokoban-opt08-strips",
    "sokoban-opt11-strips", "spider-opt18-strips", "storage", "termes-opt18-strips",
    "tetris-opt14-strips", "tidybot-opt11-strips", "tidybot-opt14-strips", "tpp",
    "transport-opt08-strips", "transport-opt11-strips", "transport-opt14-strips",
    "trucks-strips", "visitall-opt11-strips", "visitall-opt14-strips",
    "woodworking-opt08-strips", "woodworking-opt11-strips", "zenotravel",
]

DOMAIN_GROUPS = {
    "airport": ["airport"],
    "assembly": ["assembly"],
    "barman": [
        "barman", "barman-opt11-strips", "barman-opt14-strips",
        "barman-sat11-strips", "barman-sat14-strips"],
    "blocksworld": ["blocks", "blocksworld"],
    "cavediving": ["cavediving-14-adl"],
    "childsnack": ["childsnack-opt14-strips", "childsnack-sat14-strips"],
    "citycar": ["citycar-opt14-adl", "citycar-sat14-adl"],
    "depots": ["depot", "depots"],
    "driverlog": ["driverlog"],
    "elevators": [
        "elevators-opt08-strips", "elevators-opt11-strips",
        "elevators-sat08-strips", "elevators-sat11-strips"],
    "floortile": [
        "floortile-opt11-strips", "floortile-opt14-strips",
        "floortile-sat11-strips", "floortile-sat14-strips"],
    "freecell": ["freecell"],
    "ged": ["ged-opt14-strips", "ged-sat14-strips"],
    "grid": ["grid"],
    "gripper": ["gripper"],
    "hiking": ["hiking-opt14-strips", "hiking-sat14-strips"],
    "logistics": ["logistics98", "logistics00"],
    "maintenance": ["maintenance-opt14-adl", "maintenance-sat14-adl"],
    "miconic": ["miconic", "miconic-strips"],
    "miconic-fulladl": ["miconic-fulladl"],
    "miconic-simpleadl": ["miconic-simpleadl"],
    "movie": ["movie"],
    "mprime": ["mprime"],
    "mystery": ["mystery"],
    "nomystery": ["nomystery-opt11-strips", "nomystery-sat11-strips"],
    "openstacks": [
        "openstacks", "openstacks-strips", "openstacks-opt08-strips",
        "openstacks-opt11-strips", "openstacks-opt14-strips",
        "openstacks-sat08-adl", "openstacks-sat08-strips",
        "openstacks-sat11-strips", "openstacks-sat14-strips",
        "openstacks-opt08-adl", "openstacks-sat08-adl"],
    "optical-telegraphs": ["optical-telegraphs"],
    "parcprinter": [
        "parcprinter-08-strips", "parcprinter-opt11-strips", "parcprinter-sat11-strips"],
    "parking": [
        "parking-opt11-strips", "parking-opt14-strips",
        "parking-sat11-strips", "parking-sat14-strips"],
    "pathways": ["pathways"],
    "pathways-noneg": ["pathways-noneg"],
    "pegsol": ["pegsol-08-strips", "pegsol-opt11-strips", "pegsol-sat11-strips"],
    "philosophers": ["philosophers"],
    "pipes-nt": ["pipesworld-notankage"],
    "pipes-t": ["pipesworld-tankage"],
    "psr": ["psr-middle", "psr-large", "psr-small"],
    "rovers": ["rover", "rovers"],
    "satellite": ["satellite"],
    "scanalyzer": [
        "scanalyzer-08-strips", "scanalyzer-opt11-strips", "scanalyzer-sat11-strips"],
    "schedule": ["schedule"],
    "sokoban": [
        "sokoban-opt08-strips", "sokoban-opt11-strips",
        "sokoban-sat08-strips", "sokoban-sat11-strips"],
    "storage": ["storage"],
    "tetris": ["tetris-opt14-strips", "tetris-sat14-strips"],
    "thoughtful": ["thoughtful-sat14-strips"],
    "tidybot": [
        "tidybot-opt11-strips", "tidybot-opt14-strips",
        "tidybot-sat11-strips", "tidybot-sat14-strips"],
    "tpp": ["tpp"],
    "transport": [
        "transport-opt08-strips", "transport-opt11-strips", "transport-opt14-strips",
        "transport-sat08-strips", "transport-sat11-strips", "transport-sat14-strips"],
    "trucks": ["trucks", "trucks-strips"],
    "visitall": [
        "visitall-opt11-strips", "visitall-opt14-strips",
        "visitall-sat11-strips", "visitall-sat14-strips"],
    "woodworking": [
        "woodworking-opt08-strips", "woodworking-opt11-strips",
        "woodworking-sat08-strips", "woodworking-sat11-strips"],
    "zenotravel": ["zenotravel"],
    # IPC 2018:
    "agricola": ["agricola", "agricola-opt18-strips", "agricola-sat18-strips"],
    "caldera": ["caldera-opt18-adl", "caldera-sat18-adl"],
    "caldera-split": ["caldera-split-opt18-adl", "caldera-split-sat18-adl"],
    "data-network": [
        "data-network", "data-network-opt18-strips", "data-network-sat18-strips"],
    "flashfill": ["flashfill-sat18-adl"],
    "nurikabe": ["nurikabe-opt18-adl", "nurikabe-sat18-adl"],
    "organic-split": [
        "organic-synthesis-split", "organic-synthesis-split-opt18-strips",
        "organic-synthesis-split-sat18-strips"],
    "organic" : [
        "organic-synthesis", "organic-synthesis-opt18-strips",
        "organic-synthesis-sat18-strips"],
    "petri-net": [
        "petri-net-alignment", "petri-net-alignment-opt18-strips",
        "petri-net-alignment-sat18-strips"],
    "settlers": ["settlers-opt18-adl", "settlers-sat18-adl"],
    "snake": ["snake", "snake-opt18-strips", "snake-sat18-strips"],
    "spider": ["spider", "spider-opt18-strips", "spider-sat18-strips"],
    "termes": ["termes", "termes-opt18-strips", "termes-sat18-strips"],
}
# fmt: on


DOMAIN_RENAMINGS = {}
for group_name, domains in DOMAIN_GROUPS.items():
    for domain in domains:
        DOMAIN_RENAMINGS[domain] = group_name
for group_name in DOMAIN_GROUPS:
    DOMAIN_RENAMINGS[group_name] = group_name


def group_domains(run):
    old_domain = run["domain"]
    run["domain"] = DOMAIN_RENAMINGS[old_domain]
    run["problem"] = old_domain + "-" + run["problem"]
    run["id"][2] = run["problem"]
    return run


def get_repo_base() -> Path:
    """Get base directory of the repository, as an absolute path.

    Search upwards in the directory tree from the main script until a
    directory with a subdirectory named ".git" is found.

    Abort if the repo base cannot be found."""
    path = Path(SCRIPT)
    while path.parent != path:
        if (path / ".git").is_dir():
            return path
        path = path.parent
    sys.exit("repo base could not be found")


def remove_file(path: Path):
    try:
        path.unlink()
    except FileNotFoundError:
        pass


def add_evaluations_per_time(run):
    evaluations = run.get("evaluations")
    time = run.get("search_time")
    if evaluations is not None and evaluations >= 100 and time:
        run["evaluations_per_time"] = evaluations / time
    return run


def _get_exp_dir_relative_to_repo():
    repo_name = get_repo_base().name
    script = Path(SCRIPT)
    script_dir = script.parent
    rel_script_dir = script_dir.relative_to(get_repo_base())
    expname = script.stem
    return repo_name / rel_script_dir / "data" / expname


def add_scp_step(exp, login, repos_dir):
    remote_exp = Path(repos_dir) / _get_exp_dir_relative_to_repo()
    exp.add_step(
        "scp-eval-dir",
        subprocess.call,
        [
            "scp",
            "-r",  # Copy recursively.
            "-C",  # Compress files.
            f"{login}:{remote_exp}-eval",
            f"{exp.path}-eval",
        ],
    )


def fetch_algorithm(exp, expname, algo, *, new_algo=None):
    """Fetch (and possibly rename) a single algorithm from *expname*."""
    new_algo = new_algo or algo

    def rename_and_filter(run):
        if run["algorithm"] == algo:
            run["algorithm"] = new_algo
            run["id"][0] = new_algo
            return run
        return False

    exp.add_fetcher(
        f"data/{expname}-eval",
        filter=rename_and_filter,
        name=f"fetch-{new_algo}-from-{expname}",
        merge=True,
    )


def fetch_algorithms(exp, expname, *, algos=None, name=None, filters=None):
    """
    Fetch multiple or all algorithms.
    """
    assert not expname.rstrip("/").endswith("-eval")
    algos = set(algos or [])
    filters = filters or []
    if algos:

        def algo_filter(run):
            return run["algorithm"] in algos

        filters.append(algo_filter)

    exp.add_fetcher(
        f"data/{expname}-eval",
        filter=filters,
        name=name or f"fetch-from-{expname}",
        merge=True,
    )


def add_absolute_report(exp, *, name=None, outfile=None, **kwargs):
    report = AbsoluteReport(**kwargs)
    if name and not outfile:
        outfile = f"{name}.{report.output_format}"
    elif outfile and not name:
        name = Path(outfile).name
    elif not name and not outfile:
        name = f"{exp.name}-abs"
        outfile = f"{name}.{report.output_format}"

    if not Path(outfile).is_absolute():
        outfile = Path(exp.eval_dir) / outfile

    exp.add_report(report, name=name, outfile=outfile)
