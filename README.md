# MAWS

Design an aptamer — a short strand of RNA or DNA — that folds up against a
target molecule and sticks to it.

MAWS grows the strand one nucleotide at a time. At each step it tries every
nucleotide at both ends of the strand it has so far, samples thousands of
shapes for each of those candidates, scores them, and keeps the best. This is
a rewrite of the [MAWS 2023](README_orig.md) code by `dtu-denmark`, for
integration into [pyaptamer](https://github.com/gc-os-ai/pyaptamer/).

> **This is the experimental redesign branch.** It is not compatible with what
> came before: the old `Complex`, `Chain`, `Structure`, `space` and `run`
> modules are gone, and so is the code that imported them. See
> [docs/redesign](docs/redesign/) for the design and
> [docs/redesign/09-audit-response.md](docs/redesign/09-audit-response.md) for
> the scientific defects it fixes and how the answers it gives now differ.

## Installing

AmberTools is not on PyPI, so a working installation comes from conda:

```console
$ conda env create -f environment.yml
$ conda activate maws
$ pip install -e ".[dev]"
```

Without AmberTools and OpenMM the package still imports, and everything except
building real structures and computing real energies still runs — which is
what lets most of the test suite run anywhere.

## Running

```console
$ maws design --target thrombin.pdb --length 15
```

Try a short run first. The defaults perform on the order of a million energy
evaluations and take hours:

```console
$ maws design --target thrombin.pdb --length 3 --samples 50
```

The same thing from Python:

```python
from maws import design

result = design("thrombin.pdb", length=15)
print(result.sequence, result.energy, result.score)
```

`maws design --help` lists every option; each has a matching keyword argument
on `design`.

## How it is put together

Each layer uses only the ones above it.

| Layer | Modules | What lives there |
|---|---|---|
| values | `values`, `errors`, `forcefield` | Sequences, atom ranges, turnable bonds. No behaviour, no I/O. |
| structures | `topology`, `build`, `libraries`, `io` | What a design is made of, and running AmberTools to build it. |
| geometry | `pose`, `geometry`, `regrow` | Atom positions, and every way of moving them. |
| physics | `energy`, `relax` | Turning positions into an energy, and settling them. |
| search | `sampling`, `scoring`, `search` | Proposing shapes, scoring them, and the growth loop. |
| interface | `api`, `cli`, `reporting` | `design()`, `AptamerDesigner`, the command-line program. |

A `Pose` never changes: every method returns a new one. That is what makes the
search safe to write as a plain loop over candidates.

## Testing

```console
$ pytest                  # unit and algorithm tiers, about two seconds
$ pytest -m integration   # needs AmberTools and OpenMM installed
```

The suite is in three tiers, described at the top of
[tests/conftest.py](tests/conftest.py). The middle tier runs the real search
algorithm against a stand-in builder and a stand-in energy, so the whole growth
loop is exercised in milliseconds with nothing installed. The integration tier
is where anything that depends on real chemistry is checked — above all that
turning any bond of a real strand leaves every other bond intact.

Docstring examples are run as tests too, so an example that no longer produces
what it claims fails the build.

## Contributing

Docstrings follow the rules in
[.claude/skills/docstring/SKILL.md](.claude/skills/docstring/SKILL.md): NumPy
format, written for someone meeting the package for the first time.

`pre-commit install` sets up the formatter and linter that run on every commit.
