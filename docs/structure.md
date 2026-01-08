# Structure API Reference

`maws.structure.Structure` - Container for residue templates and per-residue topology rules.

## Overview

Structure defines the static chemistry used by MAWS:
- Residue metadata (atom counts)
- Alias mapping for LEaP names
- Polymer connectivity rules
- Torsion definitions

Structure is the "vocabulary" that Chain uses to translate sequences.

## Constructor

```python
Structure(
    residue_names: Sequence[str],
    residue_length: Sequence[int] = None,
    rotating_elements: Sequence[RotationSpec] = None,
    backbone_elements: Sequence[BackboneSpec] = None,
    connect: Sequence[ConnectEntry] = None,
    residue_path: str = None,
    alias: Sequence[Sequence[str]] = None,
)
```

### Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `residue_names` | `list[str]` | Template names in order (e.g., `["A", "C", "G", "U"]`) |
| `residue_length` | `list[int]` | Atom count per residue |
| `rotating_elements` | `list[tuple]` | Torsion definitions: `(residue, start, bond, end)` |
| `backbone_elements` | `list[tuple]` | Connection anchors for polymer |
| `connect` | `list` | Connectivity rules: `[[append], [prepend], len, len]` |
| `residue_path` | `str` | Directory containing `.lib`/`.frcmod` files |
| `alias` | `list[list[str]]` | Alias mapping: `[name, alone, start, middle, end]` |

## Attributes

| Attribute | Type | Description |
|-----------|------|-------------|
| `residue_names` | `list[str]` | Ordered template names |
| `residue_length` | `dict[str, int]` | Map residue → atom count |
| `alias` | `dict[str, list[str]]` | Map residue → `[alone, start, middle, end]` |
| `connect` | `dict[str, ConnectEntry]` | Map residue → connectivity |
| `rotating_elements` | `dict[str, RotListStored]` | Map residue → torsions |
| `backbone_elements` | `dict[str, BackboneEntry]` | Map residue → backbone anchors |
| `init_string` | `str` | LEaP bootstrap commands |

## Methods

### `translate(sequence: str) -> str`

Convert alias sequence to canonical LEaP names.

Rules:
- Single residue → `alone` form
- First residue → `start` form
- Middle residues → `middle` form
- Last residue → `end` form

```python
# If alias["A"] = ["A", "A5", "A", "A3"]
s.translate("A B C")  # → "A5 B C3"
s.translate("A")      # → "A"
```

### `resolve_index(residue: str, i: int) -> int`

Normalize a negative atom index to absolute.

```python
s.resolve_index("G", -1)  # → 32 (if G has 33 atoms)
```

### `torsions(residue: str) -> list[tuple]`

Get normalized rotation triples for a residue.

Returns list of `(start, bond, end_or_None)`.

```python
s.torsions("G")  # → [(0, 1, None), (1, 2, 3), ...]
```

### `append_bond(residue, prev_residue_length=None) -> tuple`

Get connection info for appending to 3′ end.

Returns `(new_atom_idx, old_atom_idx, bond_length_Å)`.

### `prepend_bond(residue, next_residue_length=None) -> tuple`

Get connection info for prepending to 5′ end.

Returns `(new_atom_idx, old_atom_idx, bond_length_Å)`.

### `add_rotation(residue_name, rotations) -> dict`

Add torsion definitions for a residue.

```python
s.add_rotation("G", (0, 1, None))  # Single torsion
s.add_rotation("G", [(0, 1, 2), (3, 4, 5)])  # Multiple
```

## Type Aliases

```python
ResidueName = str
AtomIndex = int  # 0-based, negative = from end

RotationSpec = tuple[ResidueName, AtomIndex, AtomIndex, AtomIndex | None]
BackboneSpec = tuple[ResidueName, AtomIndex, AtomIndex, AtomIndex, AtomIndex, AtomIndex]
BackboneEntry = list[list[int]]  # [[start, mid_pre, bond], [mid_post, end]]
ConnectEntry = list[list[int] | float]  # [[append], [prepend], len, len]
AliasEntry = list[str]  # [alone, start, middle, end]
```

## Pre-built Structures

MAWS provides pre-built structures for RNA and DNA:

```python
from maws.rna_structure import load_rna_structure
from maws.dna_structure import load_dna_structure

rna = load_rna_structure()  # A, C, G, U
dna = load_dna_structure()  # DA, DC, DG, DT
```

## Usage Example

```python
from maws.structure import Structure

# Define a simple structure
s = Structure(
    residue_names=["X", "Y"],
    residue_length=[10, 12],
    rotating_elements=[
        ("X", 0, 1, None),
        ("X", 1, 2, 3),
        ("Y", 0, 1, 2),
    ],
    connect=[
        [[0, -1], [-2, 0], 1.6, 1.6],
        [[0, -1], [-2, 0], 1.6, 1.6],
    ],
    alias=[
        ["X", "X", "X5", "X", "X3"],
        ["Y", "Y", "Y5", "Y", "Y3"],
    ],
    residue_path="/path/to/lib/files",
)

# Translate sequence
canonical = s.translate("X Y X")  # → "X5 Y X3"

# Get torsions
torsions = s.torsions("X")  # → [(0, 1, None), (1, 2, 3)]
```
