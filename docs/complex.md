# Complex API Reference

`maws.complex.Complex` - Container for Chain objects that builds AMBER topology and runs OpenMM simulations.

## Overview

Complex is the high-level container that:
- Owns all Chain objects
- Builds AMBER topology/coordinates using AmberTools LEaP
- Creates OpenMM simulation for energy calculations and MD
- Manages a content-addressed cache for builds

## Constructor

```python
Complex(
    force_field_aptamer: str = "leaprc.RNA.OL3",
    force_field_ligand: str = "leaprc.protein.ff19SB",
)
```

## Attributes

| Attribute | Type | Description |
|-----------|------|-------------|
| `chains` | `list[Chain]` | All chains in insertion order |
| `positions` | `list[Vec3]` | Flat coordinate array (all atoms) |
| `topology` | `Topology` | OpenMM topology |
| `system` | `System` | OpenMM system |
| `integrator` | `Integrator` | OpenMM integrator |
| `simulation` | `Simulation` | OpenMM Simulation object |
| `prmtop` | `AmberPrmtopFile` | AMBER topology file |
| `inpcrd` | `AmberInpcrdFile` | AMBER coordinate file |
| `build_string` | `str` | LEaP force field preamble |

## Methods

### Chain Management

#### `add_chain(sequence: str, structure: Structure) -> None`

Append a Chain with the given sequence and Structure.

```python
cpx.add_chain("G A U", rna_structure)
```

#### `add_chain_from_pdb(pdb_path, force_field_aptamer, force_field_ligand, structure=None, pdb_name="LIG", parameterized=False) -> None`

Add a ligand chain from a PDB file.

```python
cpx.add_chain_from_pdb(
    pdb_path="ligand.pdb",
    force_field_aptamer="leaprc.RNA.OL3",
    force_field_ligand="leaprc.protein.ff19SB",
    parameterized=True,  # For proteins
)
```

#### `get_chain(index: int) -> Chain`

Return chain by index with bounds checking.

#### `aptamer_chain() -> Chain`

Convenience: return chain[0], typically the aptamer.

#### `ligand_chain() -> Chain`

Convenience: return chain[1], typically the ligand.

### Build Methods

#### `build(target_path="", file_name="out") -> None`

Materialize chains into AMBER topology/coordinates using LEaP.

Uses content-addressed cache under `.maws_cache/`.

```python
cpx.build()
# Now cpx.positions, cpx.topology, cpx.simulation are available
```

#### `rebuild(target_path="", file_name="out", exclusion=None) -> None`

Rebuild after sequence edits, preserving coordinates outside modified regions.

```python
aptamer.append_sequence("A")
cpx.rebuild()  # Preserves existing atom positions
```

### Geometry Operations

#### `rotate_element(element, angle, reverse=False) -> None`

Rotate a global element by angle radians.

**Parameters:**
- `element`: `[start, bond, end_exclusive]` in global indices
- `angle`: radians
- `reverse`: direction flag

#### `rotate_global(element, axis, angle, reverse=False, glob=True) -> None`

Core rotation kernel using Rodrigues formula.

**Parameters:**
- `element`: atom triple
- `axis`: rotation axis (Vec3 or array with units)
- `angle`: radians
- `glob`: if True, rotate whole chain; if False, rotate subelement

```python
import numpy as np
from openmm import unit

axis = np.array([0, 0, 1]) * unit.angstrom
cpx.rotate_global(aptamer.element, axis, 0.5)
```

#### `translate_global(element, shift) -> None`

Translate a global element by a displacement vector.

```python
shift = np.array([5, 0, 0]) * unit.angstrom
cpx.translate_global(aptamer.element, shift)
```

### Energy and MD

#### `get_energy() -> tuple[float, list[Vec3]]`

Compute potential energy for current coordinates.

Returns `(energy_kJ_mol, positions)`.

```python
energy, pos = cpx.get_energy()
```

#### `minimize(max_iterations=100) -> float`

Local energy minimization.

```python
energy = cpx.minimize()
```

#### `step(number_of_steps: int) -> tuple[float, list[Vec3]]`

Advance MD by N integrator steps.

```python
energy, pos = cpx.step(1000)
```

#### `rigid_minimize(max_iterations=100, max_step_iterations=100) -> float`

Experimental: random torsions + local minimization.

#### `pert_min(size=0.1, iterations=50) -> float`

Experimental: random kicks + minimization.

## Usage Example

```python
from maws.complex import Complex
from maws.rna_structure import load_rna_structure
import numpy as np
from openmm import unit

# Create complex
cpx = Complex()
rna = load_rna_structure()
cpx.add_chain("", rna)
cpx.add_chain_from_pdb("protein.pdb", ...)

# Build
aptamer = cpx.aptamer_chain()
aptamer.create_sequence("G A U C")
cpx.build()

# Sample conformations
for _ in range(100):
    # Random translation
    shift = np.random.uniform(-10, 10, 3) * unit.angstrom
    cpx.translate_global(aptamer.element, shift)
    
    # Minimize
    cpx.minimize()
    
    # Get energy
    energy, _ = cpx.get_energy()
    
    # Reset for next sample
    cpx.rebuild()
```

## Caching

Builds are cached under `.maws_cache/`:
- Cache key = SHA1 of build_string + sequences
- Files: `<key>.prmtop`, `<key>.inpcrd`

Set `MAWS_OPENMM_PLATFORM` env var to override GPU selection.
