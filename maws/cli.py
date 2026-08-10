"""
maws.cli
========

The ``maws`` command.

Installing this package puts a ``maws`` command on your path with four
subcommands:

``maws design``
    The main one: design an aptamer against a target molecule.
``maws prepare``
    Work out the parameters for a target molecule and stop there, writing the
    files a later design run will reuse.
``maws inspect``
    Say what a design run would do — how large the molecules are, which
    parameters would be used, how much work it would be — without doing any of
    it.
``maws clean``
    Tidy a downloaded structure file into something the builder will accept.

Settings can also be put in a file and passed with ``--config``, which saves
typing a dozen flags. A flag given on the command line always wins over the
file.

Examples
--------
.. code-block:: console

    $ maws design --target data/pfoa.pdb --length 20 --aptamer DNA --seed 42
    $ maws inspect --target data/pfoa.pdb --aptamer DNA
    $ maws design --config maws.toml --length 25
"""

from __future__ import annotations

import argparse
import logging
import sys
import tomllib
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from maws.errors import MawsError

__all__ = ["main"]

PROGRAM = "maws"
"""Name the command is installed under."""

CONFIG_SECTIONS = ("design", "sampling", "physics", "output")
"""Section names read from a ``--config`` file.

Settings are grouped only to keep the file readable; every section is merged
into one flat set of values, so it does not matter which section a setting
appears in.
"""


def main(argv: Sequence[str] | None = None) -> int:
    """Run the ``maws`` command.

    Parameters
    ----------
    argv : sequence of str, optional
        The arguments, without the program name. Defaults to the ones the
        command was invoked with.

    Returns
    -------
    int
        ``0`` if the command succeeded, ``1`` if it reported a problem, ``2``
        if the arguments were wrong.

    Examples
    --------
    >>> main(["--version"])
    Traceback (most recent call last):
        ...
    SystemExit: 0
    """
    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(level=_log_level(args), format="%(message)s", stream=sys.stderr)
    try:
        return int(args.handler(args))
    except MawsError as error:
        print(f"{PROGRAM}: {error}", file=sys.stderr)
        return 1


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser for the whole command.

    Returns
    -------
    argparse.ArgumentParser
        A parser whose ``handler`` attribute, once arguments are parsed, is the
        function that carries out the chosen subcommand.

    Examples
    --------
    >>> parser = build_parser()
    >>> args = parser.parse_args(["design", "--target", "t.pdb", "--length", "4"])
    >>> args.target, args.length
    ('t.pdb', 4)
    """
    from maws import __version__

    parser = argparse.ArgumentParser(
        prog=PROGRAM,
        description="MAWS - design aptamers against a target molecule.",
    )
    parser.add_argument(
        "--version", action="version", version=f"{PROGRAM} {__version__}"
    )
    subcommands = parser.add_subparsers(dest="command", required=True)

    _add_design(subcommands)
    _add_prepare(subcommands)
    _add_inspect(subcommands)
    _add_clean(subcommands)
    return parser


def _add_verbosity(parser: argparse.ArgumentParser) -> None:
    """Add the mutually exclusive ``-v``/``-q`` flags to a subcommand.

    Parameters
    ----------
    parser : argparse.ArgumentParser
        The subcommand's parser.
    """
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="report each candidate as it is scored",
    )
    group.add_argument(
        "-q", "--quiet", action="store_true", help="report nothing but the result"
    )


def _add_design(subcommands: Any) -> None:
    """Add the ``design`` subcommand.

    Parameters
    ----------
    subcommands : argparse._SubParsersAction
        Where to add it.
    """
    parser = subcommands.add_parser(
        "design",
        help="design an aptamer against a target (the main workflow)",
        description="Design an aptamer against a target molecule.",
    )
    required = parser.add_argument_group("required")
    required.add_argument(
        "--target",
        metavar="PDB",
        help="structure file for the molecule the aptamer should bind to",
    )

    output = parser.add_argument_group("output")
    output.add_argument(
        "-n", "--name", default="MAWS_aptamer", help="stem for the output filenames"
    )
    output.add_argument(
        "-o",
        "--output",
        metavar="DIR",
        default=".",
        help="directory to write all output files to",
    )
    output.add_argument(
        "--format",
        choices=("text", "json"),
        default="text",
        help="how to print the result",
    )

    shape = parser.add_argument_group("design")
    shape.add_argument(
        "--length", type=int, default=15, help="nucleotides in the finished strand"
    )
    shape.add_argument(
        "--aptamer", choices=("RNA", "DNA"), default="RNA", help="nucleic acid to use"
    )
    shape.add_argument(
        "--molecule",
        choices=("protein", "organic", "lipid"),
        default="protein",
        help="what the target is",
    )

    search = parser.add_argument_group("sampling")
    search.add_argument(
        "--samples", type=int, default=5000, help="shapes tried per growth step"
    )
    search.add_argument(
        "--first-samples",
        type=int,
        default=None,
        help="shapes tried on the first step [same as --samples]",
    )
    search.add_argument(
        "--beta", type=float, default=0.01, help="how sharply low energies are favoured"
    )
    search.add_argument(
        "--reach",
        type=float,
        default=10.0,
        help="how far past the target the strand may sit, in angstrom",
    )
    search.add_argument(
        "--probe",
        type=float,
        default=1.4,
        help="radius of the ball used to find the target surface, in angstrom",
    )
    search.add_argument(
        "--sampling",
        choices=("sphere", "surface-following"),
        default="sphere",
        help="which region to draw positions from",
    )
    search.add_argument(
        "--seed",
        type=int,
        default=None,
        help="fix the randomness, for a repeatable run",
    )
    search.add_argument(
        "--relax-iterations",
        type=int,
        default=50,
        help="nudge-and-settle rounds after each nucleotide is joined on",
    )

    physics = parser.add_argument_group("physics")
    physics.add_argument(
        "--salt-conc",
        type=float,
        default=0.15,
        help="dissolved salt concentration in mol/L",
    )

    cleaning = parser.add_argument_group("input cleaning")
    cleaning.add_argument(
        "--clean-pdb", action="store_true", help="tidy the target file before building"
    )
    cleaning.add_argument(
        "--keep-chains", default="all", help="chains to keep: all, one, or e.g. A,B"
    )
    cleaning.add_argument(
        "--remove-h", action="store_true", help="strip hydrogens while cleaning"
    )
    cleaning.add_argument(
        "--drop-hetatm", action="store_true", help="drop water, ions and other extras"
    )

    parser.add_argument(
        "--config", metavar="FILE", help="read settings from a TOML file"
    )
    _add_verbosity(parser)
    # Kept so resolve_settings can tell an option that was typed from one that
    # was left alone: only the latter may be overridden by a config file.
    parser.set_defaults(
        handler=_run_design,
        design_defaults={
            action.dest: action.default
            for action in parser._actions  # noqa: SLF001
        },
    )


def _add_prepare(subcommands: Any) -> None:
    """Add the ``prepare`` subcommand.

    Parameters
    ----------
    subcommands : argparse._SubParsersAction
        Where to add it.
    """
    parser = subcommands.add_parser(
        "prepare",
        help="work out a target's parameters without designing anything",
        description=(
            "Run the parameter-fitting programs for a target molecule and "
            "write the files, so a later design run does not have to."
        ),
    )
    parser.add_argument("--target", required=True, metavar="PDB", help="structure file")
    parser.add_argument(
        "--molecule",
        choices=("protein", "organic", "lipid"),
        default="protein",
        help="what the target is",
    )
    parser.add_argument(
        "--aptamer", choices=("RNA", "DNA"), default="RNA", help="nucleic acid to use"
    )
    parser.add_argument(
        "-o", "--output", metavar="DIR", default="./params", help="where to write them"
    )
    _add_verbosity(parser)
    parser.set_defaults(handler=_run_prepare)


def _add_inspect(subcommands: Any) -> None:
    """Add the ``inspect`` subcommand.

    Parameters
    ----------
    subcommands : argparse._SubParsersAction
        Where to add it.
    """
    parser = subcommands.add_parser(
        "inspect",
        help="show what a run would do, without running it",
        description="Describe a design run without building or scoring anything.",
    )
    parser.add_argument("--target", required=True, metavar="PDB", help="structure file")
    parser.add_argument(
        "--aptamer", choices=("RNA", "DNA"), default="RNA", help="nucleic acid to use"
    )
    parser.add_argument(
        "--molecule",
        choices=("protein", "organic", "lipid"),
        default="protein",
        help="what the target is",
    )
    parser.add_argument(
        "--length", type=int, default=15, help="nucleotides in the finished strand"
    )
    parser.add_argument(
        "--samples", type=int, default=5000, help="shapes tried per growth step"
    )
    parser.add_argument(
        "--salt-conc", type=float, default=0.15, help="dissolved salt in mol/L"
    )
    _add_verbosity(parser)
    parser.set_defaults(handler=_run_inspect)


def _add_clean(subcommands: Any) -> None:
    """Add the ``clean`` subcommand.

    Parameters
    ----------
    subcommands : argparse._SubParsersAction
        Where to add it.
    """
    parser = subcommands.add_parser(
        "clean",
        help="tidy a structure file so the builder will accept it",
        description=(
            "A downloaded structure often holds several copies of the "
            "molecule, bound water, and records the builder cannot use. This "
            "writes a tidied copy alongside the original."
        ),
    )
    parser.add_argument("--target", required=True, metavar="PDB", help="structure file")
    parser.add_argument(
        "--keep-chains", default="all", help="chains to keep: all, one, or e.g. A,B"
    )
    parser.add_argument("--remove-h", action="store_true", help="strip hydrogens")
    parser.add_argument(
        "--drop-hetatm", action="store_true", help="drop water, ions and other extras"
    )
    _add_verbosity(parser)
    parser.set_defaults(handler=_run_clean)


def load_config(path: str | Path) -> dict[str, Any]:
    """Read settings from a TOML file and flatten them into one mapping.

    Parameters
    ----------
    path : str or pathlib.Path
        The file to read.

    Returns
    -------
    dict
        Setting name to value, with any hyphens in names turned into
        underscores so they match the command's own option names.

    Raises
    ------
    maws.errors.ConfigurationError
        If the file cannot be read or is not valid TOML.

    Notes
    -----
    Only the sections named in :data:`CONFIG_SECTIONS` are read. Anything else
    is ignored, so a file may carry notes or settings for other tools.
    """
    from maws.errors import ConfigurationError

    try:
        with open(path, "rb") as handle:
            raw = tomllib.load(handle)
    except OSError as exc:
        raise ConfigurationError(f"cannot read config file {path}: {exc}") from exc
    except tomllib.TOMLDecodeError as exc:
        raise ConfigurationError(f"{path} is not valid TOML: {exc}") from exc

    flat: dict[str, Any] = {}
    for section in CONFIG_SECTIONS:
        for key, value in raw.get(section, {}).items():
            flat[key.replace("-", "_")] = value
    return flat


def resolve_settings(
    args: argparse.Namespace, parser_defaults: dict[str, Any]
) -> dict[str, Any]:
    """Merge command-line flags, a config file, and the built-in defaults.

    Parameters
    ----------
    args : argparse.Namespace
        What was parsed from the command line.
    parser_defaults : dict
        The value each option takes when it is not given.

    Returns
    -------
    dict
        The settled values.

    Notes
    -----
    A flag typed on the command line wins over the config file, which in turn
    wins over the defaults. A flag is taken to have been typed when its value
    differs from the default, which is why the defaults have to be passed in.

    Examples
    --------
    >>> import argparse
    >>> args = argparse.Namespace(length=15, config=None)
    >>> resolve_settings(args, {"length": 15})["length"]
    15
    """
    settings = dict(vars(args))
    from_file = load_config(settings["config"]) if settings.get("config") else {}

    for key, value in from_file.items():
        if key not in settings:
            continue
        if settings[key] == parser_defaults.get(key):
            settings[key] = value
    return settings


def _log_level(args: argparse.Namespace) -> int:
    """Choose a logging level from the verbosity flags.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed arguments, which may carry ``verbose`` and ``quiet``.

    Returns
    -------
    int
        A level from :mod:`logging`.
    """
    if getattr(args, "quiet", False):
        return logging.ERROR
    if getattr(args, "verbose", False):
        return logging.DEBUG
    return logging.INFO


def _run_design(args: argparse.Namespace) -> int:
    """Carry out ``maws design``.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed arguments.

    Returns
    -------
    int
        ``0`` on success, ``2`` if no target was given.
    """
    import json

    from maws.api import design
    from maws.reporting import LogReporter

    settings = resolve_settings(args, getattr(args, "design_defaults", {}))

    if not settings.get("target"):
        print(f"{PROGRAM} design: --target is required", file=sys.stderr)
        return 2

    output = Path(settings["output"])
    output.mkdir(parents=True, exist_ok=True)

    with LogReporter(
        job=settings["name"],
        directory=output,
        logger=logging.getLogger(PROGRAM),
    ) as reporter:
        result = design(
            settings["target"],
            length=settings["length"],
            aptamer=settings["aptamer"],
            molecule=settings["molecule"],
            samples=settings["samples"],
            first_samples=settings["first_samples"],
            beta=settings["beta"],
            salt_conc=settings["salt_conc"],
            reach=settings["reach"],
            probe=settings["probe"],
            sampling=settings["sampling"],
            relax_iterations=settings["relax_iterations"],
            seed=settings["seed"],
            clean_pdb=settings["clean_pdb"],
            keep_chains=settings["keep_chains"],
            remove_h=settings["remove_h"],
            drop_hetatm=settings["drop_hetatm"],
            on_event=reporter,
        )

    if settings["format"] == "json":
        print(
            json.dumps(
                {
                    "sequence": result.sequence,
                    "energy": result.energy,
                    "entropy": result.entropy,
                    "steps": result.steps,
                    "success": result.success,
                }
            )
        )
    else:
        print(result)
    return 0 if result.success else 1


def _run_prepare(args: argparse.Namespace) -> int:
    """Carry out ``maws prepare``.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed arguments.

    Returns
    -------
    int
        ``0`` on success.
    """
    from maws.build import LeapBuilder
    from maws.forcefield import ForceField
    from maws.topology import PdbChain, default_residue_name

    target = Path(args.target)
    forcefield = ForceField.for_target(args.aptamer, args.molecule)
    builder = LeapBuilder(params_dir=args.output)
    chain = builder.prepare(
        PdbChain(
            role="ligand",
            path=target,
            residue_name=default_residue_name(target),
            parameterized=forcefield.parameterized,
        ),
        forcefield,
    )
    name = chain.canonical[0]
    print(f"Wrote {Path(args.output) / f'{name}.lib'} ({chain.n_atoms} atoms)")
    return 0


def _run_inspect(args: argparse.Namespace) -> int:
    """Carry out ``maws inspect``.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed arguments.

    Returns
    -------
    int
        ``0`` on success.

    Notes
    -----
    Counts the target's atoms by reading its file directly rather than by
    building anything, so this works without AmberTools installed.
    """
    from maws.build import FakeBuilder
    from maws.forcefield import ForceField
    from maws.libraries import dna, rna
    from maws.topology import PdbChain, default_residue_name

    target = Path(args.target)
    forcefield = ForceField.for_target(
        args.aptamer, args.molecule, salt_conc=args.salt_conc
    )
    library = rna() if args.aptamer == "RNA" else dna()
    residue = default_residue_name(target)
    counted = FakeBuilder().prepare(
        PdbChain(
            role="ligand",
            path=target,
            residue_name=residue,
            parameterized=forcefield.parameterized,
        ),
        forcefield,
    )

    evaluations = 8 * args.samples * max(args.length - 1, 0) + 4 * args.samples
    print("Assembly")
    print(f"  aptamer   {args.aptamer:<8} {forcefield.aptamer_source:<24} empty")
    print(
        f"  ligand    PDB      {forcefield.ligand_source:<24} "
        f"{target.name} ({counted.n_atoms} atoms)"
    )
    print("Force field")
    print(
        f"  salt {forcefield.salt_conc} mol/L    "
        f"parameters ready: {'yes' if forcefield.parameterized else 'no'}    "
        f"nucleotides: {forcefield.alphabet}"
    )
    print("Residues")
    print(f"  library   {len(library)} residues, {len(library.tokens)} tokens")
    print("Estimated work")
    print(f"  energy evaluations  {evaluations:,}")
    return 0


def _run_clean(args: argparse.Namespace) -> int:
    """Carry out ``maws clean``.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed arguments.

    Returns
    -------
    int
        ``0`` on success.
    """
    from maws.io.pdb_cleaner import resolve_pdb_path

    cleaned, original = resolve_pdb_path(
        args.target,
        "protein",
        clean_pdb=True,
        keep_chains=args.keep_chains,
        remove_h=args.remove_h,
        drop_hetatm=args.drop_hetatm,
        logger=logging.getLogger(PROGRAM),
    )
    if cleaned == original:
        print(f"{original} was left unchanged")
    else:
        print(f"Wrote {cleaned}")
    return 0


if __name__ == "__main__":  # pragma: no cover - module entry point
    raise SystemExit(main())
