# Python API Reference

`maws.run` - Programmatic interface for running MAWS.

## Overview

The `maws.run` module provides a clean Python API for running MAWS without using the CLI. It uses dataclasses for configuration and results.

## Classes

### MAWSConfig

Configuration for a MAWS aptamer design run.

```python
from maws import MAWSConfig

config = MAWSConfig(
    pdb_path="ligand.pdb",
    num_nucleotides=15,
    name="my_aptamer",
    aptamer_type="RNA",
    molecule_type="protein",
)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `pdb_path` | `str` | *required* | Path to ligand PDB file |
| `num_nucleotides` | `int` | *required* | Number of nucleotides |
| `name` | `str` | `"MAWS_aptamer"` | Job name (for output files) |
| `aptamer_type` | `"RNA"` or `"DNA"` | `"RNA"` | Aptamer type |
| `molecule_type` | `"protein"`, `"organic"`, `"lipid"` | `"protein"` | Ligand type |
| `beta` | `float` | `0.01` | Inverse temperature for sampling |
| `first_chunk_size` | `int` | `5000` | Samples in first step |
| `second_chunk_size` | `int` | `5000` | Samples in subsequent steps |
| `clean_pdb` | `bool` | `False` | Clean input PDB |
| `keep_chains` | `str` | `"all"` | Chain policy for cleaner |
| `remove_h` | `bool` | `False` | Remove hydrogens |
| `drop_hetatm` | `bool` | `False` | Drop HETATM records |
| `verbose` | `bool` | `True` | Log progress to console |

### MAWSResult

Result of a MAWS aptamer design run.

| Attribute | Type | Description |
|-----------|------|-------------|
| `sequence` | `str` | Best aptamer sequence found |
| `energy` | `float` | Final energy |
| `entropy` | `float` | Final entropy score |
| `complex` | `Complex` | OpenMM Complex with final structure |
| `config` | `MAWSConfig` | Configuration used |
| `pdb_path` | `str` | Path to result PDB file |

#### Methods

**`save_pdb(path: str) -> str`**: Save result structure to PDB file.

## Functions

### run_maws

Run the MAWS aptamer design algorithm.

```python
from maws import run_maws, MAWSConfig

config = MAWSConfig(
    pdb_path="data/ligand.pdb",
    num_nucleotides=15,
)

result = run_maws(config)

print(result.sequence)  # Best aptamer sequence
result.save_pdb("final.pdb")
```

## Usage Examples

### Basic Run

```python
from maws import run_maws, MAWSConfig

# Configure
config = MAWSConfig(
    pdb_path="protein.pdb",
    num_nucleotides=20,
    aptamer_type="RNA",
    molecule_type="protein",
)

# Run
result = run_maws(config)

# Results
print(f"Best sequence: {result.sequence}")
print(f"Energy: {result.energy}")
print(f"Saved to: {result.pdb_path}")
```

### With PDB Cleaning

```python
config = MAWSConfig(
    pdb_path="messy_protein.pdb",
    num_nucleotides=15,
    clean_pdb=True,
    remove_h=True,
    drop_hetatm=True,
    keep_chains="A,B",
)

result = run_maws(config)
```

### DNA Aptamer for Small Molecule

```python
config = MAWSConfig(
    pdb_path="drug.pdb",
    num_nucleotides=25,
    aptamer_type="DNA",
    molecule_type="organic",
)

result = run_maws(config)
```

### Quick Test Run

```python
config = MAWSConfig(
    pdb_path="ligand.pdb",
    num_nucleotides=3,
    first_chunk_size=10,
    second_chunk_size=10,
    verbose=False,
)

result = run_maws(config)
```

## Output Files

- `{name}_output.log` - Log file
- `{name}_RESULT.pdb` - Final structure PDB

## Environment Variables

- `MAWS_OPENMM_PLATFORM` - Override OpenMM platform selection (e.g., "CUDA", "CPU")
