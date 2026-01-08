# Chain API Reference

`maws.chain.Chain` - A single polymer chain inside a Complex.

## Overview

Chain is a lightweight view into a Complex's coordinate array. It manages sequence state and provides methods for per-residue rotations. Geometry operations are delegated to the owning Complex.

## Constructor

```python
Chain(
    Complex: Complex,      # Owning complex
    Structure: Structure,  # Residue templates
    sequence: str = None,  # Initial alias sequence
    start: int = 0,        # Global atom index
    ID: int = 0,           # Chain identifier
)
```

## Attributes

| Attribute | Type | Description |
|-----------|------|-------------|
| `id` | `int` | Chain index within Complex.chains |
| `start` | `int` | Global atom index where chain begins |
| `length` | `int` | Total atoms in this chain |
| `sequence` | `str` | Space-separated canonical residue names |
| `alias_sequence` | `str` | Space-separated alias tokens |
| `sequence_array` | `list[str]` | Tokenized canonical sequence |
| `alias_sequence_array` | `list[str]` | Tokenized alias sequence |
| `structure` | `Structure` | Residue template bank |
| `complex` | `Complex` | Back-reference to owner |
| `element` | `list[int]` | Triple `[start, start+1, start+length]` |
| `residues_start` | `list[int]` | Per-residue atom offsets |
| `start_history` | `int` | Pre-edit start index |
| `length_history` | `int` | Pre-edit length |
| `append_history` | `list[str]` | Residues appended in last edit |
| `prepend_history` | `list[str]` | Residues prepended in last edit |

## Methods

### Sequence Manipulation

#### `create_sequence(sequence: str) -> None`

Overwrite the entire chain sequence with alias tokens.

```python
aptamer.create_sequence("G A U C")
# Must call complex.build() or complex.rebuild() afterwards
```

#### `append_sequence(sequence: str) -> None`

Append one alias residue to the 3′ end.

```python
aptamer.append_sequence("A")
# Populates append_history, must rebuild afterwards
```

#### `prepend_sequence(sequence: str) -> None`

Prepend one alias residue to the 5′ end.

```python
aptamer.prepend_sequence("G")
# Populates prepend_history, must rebuild afterwards
```

### Rotation Methods

#### `rotate_in_residue(residue_index, residue_element_index, angle, reverse=False)`

Rotate a template-defined torsion inside a specific residue.

**Parameters:**
- `residue_index` (int): Residue index within chain (negative OK)
- `residue_element_index` (int): Torsion index from Structure
- `angle` (float): Rotation angle in radians
- `reverse` (bool): Whether to reverse rotation direction

```python
# Rotate first torsion in first residue by 0.5 radians
aptamer.rotate_in_residue(0, 0, 0.5)
```

#### `rotate_element(element, angle, reverse=False)`

Rotate a chain-local element (atom triple) by angle radians.

#### `rotate_historic_element(historic_element, angle)`

Rotate using element indices captured before a prepend/append.

#### `rotate_in_historic_residue(historic_index, element_index, angle)`

Rotate using a residue index from before a prepend operation.

### Internal Methods

#### `update_chains() -> None`

Recompute length and offsets, then shift downstream chains.

Called internally after sequence changes.

## Usage Example

```python
from maws.complex import Complex
from maws.rna_structure import load_rna_structure

# Create complex and add chain
cpx = Complex()
rna = load_rna_structure()
cpx.add_chain("", rna)

# Get aptamer chain
aptamer = cpx.aptamer_chain()

# Set sequence
aptamer.create_sequence("G A U")
cpx.build()

# Rotate torsion in first residue
aptamer.rotate_in_residue(0, 0, 0.5)

# Extend sequence
aptamer.append_sequence("C")
cpx.rebuild()
```
