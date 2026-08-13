# Python API Reference

`maws.run` - Programmatic interface for running MAWS.

## Overview

The `maws.run` module provides a clean Python API for running MAWS without using the CLI. Construct a [`MawsRunner`](#mawsrunner) with the design parameters, then call [`run`](#mawsrunnerrun) once per ligand. It returns a [`MawsResult`](#mawsresult) dataclass.

Both names are re-exported from the package root:

```python
from maws import MawsRunner, MawsResult
```

## Classes

### MawsRunner

Holds the configuration for a MAWS aptamer design run. All constructor arguments are **keyword-only**.

```python
from maws import MawsRunner

runner = MawsRunner(
    num_nucleotides=15,
    aptamer_type="RNA",
    molecule_type="protein",
)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `num_nucleotides` | `int` | *required* | Number of nucleotides to design |
| `aptamer_type` | `"RNA"` or `"DNA"` | *required* | Aptamer type |
| `molecule_type` | `"protein"`, `"organic"`, `"lipid"` | *required* | Ligand type |
| `beta` | `float` | `0.01` | Inverse temperature for the entropy score (see [docs/routines.md](routines.md)) |
| `first_chunk_size` | `int` | `5000` | Samples in first step |
| `second_chunk_size` | `int` | `5000` | Samples in subsequent steps |
| `clean_pdb` | `bool` | `False` | Clean input PDB |
| `keep_chains` | `str` | `"all"` | Chain policy for cleaner |
| `remove_h` | `bool` | `False` | Remove hydrogens |
| `drop_hetatm` | `bool` | `False` | Drop HETATM records |
| `verbose` | `bool` | `False` | Log progress at INFO instead of DEBUG |
| `sampler_mode` | `"sphere"` or `"surface-following"` | `"sphere"` | Shape of the region poses are drawn from |
| `reach` | `float` | `10.0` | Distance the envelope extends past the ligand surface (Å) |
| `d_max` | `float` | `6.0` | How far from the surface a pose may sit (Å), for `surface-following` |
| `site_centre` | `sequence of 3 floats` | `None` | Sample around this point instead of the whole target (Å) |
| `site_radius` | `float` | `15.0` | How far the region reaches from `site_centre` (Å) |
| `probe` | `float` | `1.4` | vdW probe radius for SAS rejection (Å, water-equivalent) |
| `clash_tolerance` | `float` | `1.0` | How far the placed strand may overlap the target's vdW spheres before the pose is redrawn (Å) |
| `salt_conc` | `float` | `0.15` | Monovalent salt conc. (mol/L) for GB Debye–Hückel screening; `0` = unscreened |
| `seed` | `int` | `None` | Seed for every random draw. Omit it and the run draws one, logs it, and returns it on `MawsResult` |

The constructor raises `ValueError` for `num_nucleotides <= 0`, a non-positive chunk size, or a negative `reach`, `probe`, `clash_tolerance`, or `salt_conc`, and `TypeError` for a non-integer `seed`. The sampler arguments are checked when `run()` builds the sampler.

The ligand PDB path is **not** a constructor argument — it is passed to `run()`, so one configured runner can be reused across several ligands.

> **Note (2026): behavior change.** `salt_conc` defaults to `0.15` mol/L
> (~physiological) so the GB implicit solvent now screens the highly charged
> nucleic-acid backbone. Earlier releases ran unscreened (`0.0` mol/L), so
> computed energies — and therefore the sequences MAWS selects — will differ
> from prior versions unless you pass `salt_conc=0.0` (`--salt-conc 0` on the
> CLI). The screening is monovalent only and does not model Mg²⁺.

> **Note (2026):** the initial pose search restricts samples to the region
> just outside the ligand surface, computed automatically from the ligand
> atoms. `reach` and `probe` tune that region; defaults work for most
> targets. Set `sampler_mode="surface-following"` to keep poses within
> `d_max` of the surface instead, or give a `site_centre` to spend the whole
> sample budget on a known binding site. See [docs/space.md](space.md) for
> the underlying API.

> **Note (2026):** a pose whose atoms overlap the target is thrown away and
> redrawn before it costs an energy evaluation. `clash_tolerance` sets how
> much overlap counts as a clash: the default 1.0 Å admits hydrogen bonds,
> which sit about 0.8 Å inside the summed van der Waals radii.

> **Note (2026):** every run is reproducible. Pass `seed` to repeat one; omit
> it and the run draws a seed, logs it as `Random seed: ...`, and returns it
> on `MawsResult`, so a run can be repeated after the fact.

#### MawsRunner.run

```python
run(*, pdb, name="MAWS_aptamer", output_pdb=None) -> MawsResult
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `pdb` | `str \| Path` | *required* | Input ligand PDB file path |
| `name` | `str` | `"MAWS_aptamer"` | Run name, used for log context and artifact naming |
| `output_pdb` | `str \| Path \| None` | `None` | Where to write the result structure. If `None`, no PDB is written |

`output_pdb` is interpreted two ways: an existing **directory** gets `{name}_RESULT.pdb` written inside it; anything else is treated as the exact output file path, with parent directories created as needed.

### MawsResult

Frozen dataclass returned by `MawsRunner.run`.

| Attribute | Type | Description |
|-----------|------|-------------|
| `sequence` | `str` | Best aptamer sequence found |
| `energy` | `float` | Energy of the final best configuration; `nan` if no candidate was scored |
| `entropy` | `float` | Entropy score used for selection (`maws.routines.entropy_score`; ≤ 0, lower is better) |
| `pdb_path` | `str \| None` | Path to the written PDB, or `None` when `output_pdb` was not given |
| `seed` | `int \| None` | The seed the run used. Pass it back as `MawsRunner(seed=...)` to reproduce this result |

## Usage Examples

### Basic Run

```python
from maws import MawsRunner

runner = MawsRunner(
    num_nucleotides=20,
    aptamer_type="RNA",
    molecule_type="protein",
)

result = runner.run(pdb="protein.pdb", name="my_aptamer", output_pdb="results/")

print(f"Best sequence: {result.sequence}")
print(f"Energy: {result.energy}")
print(f"Entropy: {result.entropy}")
print(f"Saved to: {result.pdb_path}")
```

### With PDB Cleaning

```python
runner = MawsRunner(
    num_nucleotides=15,
    aptamer_type="RNA",
    molecule_type="protein",
    clean_pdb=True,
    remove_h=True,
    drop_hetatm=True,
    keep_chains="A,B",
)

result = runner.run(pdb="messy_protein.pdb")
```

### DNA Aptamer for Small Molecule

```python
runner = MawsRunner(
    num_nucleotides=25,
    aptamer_type="DNA",
    molecule_type="organic",
)

result = runner.run(pdb="drug.pdb")
```

### Quick Test Run

```python
runner = MawsRunner(
    num_nucleotides=3,
    aptamer_type="RNA",
    molecule_type="protein",
    first_chunk_size=10,
    second_chunk_size=10,
)

result = runner.run(pdb="ligand.pdb")
```

### Screening Several Ligands

One runner can be reused, since the PDB is an argument to `run()`:

```python
runner = MawsRunner(
    num_nucleotides=15,
    aptamer_type="RNA",
    molecule_type="protein",
)

results = {
    target: runner.run(pdb=f"{target}.pdb", name=target, output_pdb="results/")
    for target in ("1BRQ", "1HAO")
}

for target, result in sorted(results.items(), key=lambda kv: kv[1].entropy):
    print(target, result.sequence, result.entropy)
```

## Logging

`MawsRunner` logs through the standard `logging` module under the `maws.run` logger; it does not configure handlers itself. `verbose=True` promotes a handful of progress messages from DEBUG to INFO. To see output, configure logging in your own code:

```python
import logging

logging.basicConfig(level=logging.INFO)
```

## Output Files

- `{name}_RESULT.pdb` - Final structure PDB, written only when `output_pdb` is given

The CLI (`maws.maws2023`) additionally writes `{name}_output.log`.

## Environment Variables

- `MAWS_OPENMM_PLATFORM` - Override OpenMM platform selection (e.g., "CUDA", "CPU")
