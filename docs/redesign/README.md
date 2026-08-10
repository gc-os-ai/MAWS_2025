# MAWS API Redesign

**Start here:** read [01-values.md](01-values.md). It is the change everything else
depends on.

**Status:** design only. No code written yet.
**Total effort:** about 3 weeks of focused work, split into 10 PRs.
**Smallest useful first step:** fix one bug in `chain.py`. About 2 hours. See
[07-migration.md](07-migration.md).

---

## The problem in one table

People call this codebase "a spaghetti of classes and subclasses." That is not
quite right. There is almost no inheritance here. `Structure`, `Chain`, and
`Complex` are three plain classes. Nothing subclasses anything.

Three real problems make it hard to read:

| # | Problem | Where | Fixed by |
| --- | --- | --- | --- |
| R1 | One list `[start, bond, end]` means 5 different things | everywhere | [01-values.md](01-values.md) |
| R2 | `Chain` and `Complex` point at each other and edit each other | [chain.py:149](../../maws/chain.py#L149) | [03-pose.md](03-pose.md) |
| R3 | Edit history is stored in 6 loose fields, not a type | `*_history` fields | [03-pose.md](03-pose.md) |

Two more problems follow from those three:

| # | Problem | Cost | Fixed by |
| --- | --- | --- | --- |
| C1 | The search algorithm is written out twice | 480 duplicate lines | [05-search.md](05-search.md) |
| C2 | Nothing runs in tests without AmberTools installed | slow CI, no coverage | [04-energy.md](04-energy.md) |

---

## The layers

One idea: keep chemistry, coordinates, and search in separate layers. Right now
`Complex` does all three.

```mermaid
flowchart TB
    L4["Layer 4 — Interface<br/>maws.cli · design()"]
    L3["Layer 3 — Search<br/>grow_aptamer() · Scorer"]
    L2["Layer 2 — State<br/>Pose · ChainView · EnergyModel"]
    L1["Layer 1 — Build<br/>Assembly · BuiltSystem · LeapBuilder"]
    L0["Layer 0 — Values<br/>AtomRange · Torsion · ResidueTemplate"]

    L4 --> L3 --> L2 --> L1 --> L0

    style L0 fill:#e8f4ea,stroke:#4a7c59
    style L1 fill:#eef2f7,stroke:#5a7a9a
    style L2 fill:#fdf3e3,stroke:#b8860b
    style L3 fill:#fae8e8,stroke:#a85454
    style L4 fill:#f0e8f7,stroke:#7a5a9a
```

Each layer may only use layers below it. That rule gives you:

- Layer 0 imports nothing but `dataclasses` and `numpy`.
- Layer 1 is the only layer that reads or writes files.
- Layer 2 is the only layer that imports OpenMM.
- Layer 3 uses protocols, not concrete classes. So it runs in tests in
  milliseconds.

---

## Read in this order

| Doc | What it covers | Read time |
| --- | --- | --- |
| [00-principles.md](00-principles.md) | 6 rules the other docs follow | 8 min |
| [01-values.md](01-values.md) | `AtomRange`, `Torsion`, `ResidueTemplate`, `NucleotideSequence` | 12 min |
| [02-topology.md](02-topology.md) | `ForceField`, `Assembly`, `BuiltSystem`, `LeapBuilder` | 12 min |
| [03-pose.md](03-pose.md) | `Pose`, `ChainView`, `Edit` | 15 min |
| [04-energy.md](04-energy.md) | `EnergyModel` protocol, `StubEnergy` | 8 min |

Then:

| Doc | What it covers | Read time |
| --- | --- | --- |
| [05-search.md](05-search.md) | `grow_aptamer()` and typed events | 10 min |
| [06-cli.md](06-cli.md) | CLI and the Python API | 10 min |
| [07-migration.md](07-migration.md) | 10 PRs in order, with time estimates | 10 min |
| [08-patterns.md](08-patterns.md) | What the design is called; OpenMM/sklearn conventions | 10 min |

---

## Before and after

The clearest summary is what a caller has to type.

### Now

```python
import copy
from openmm import unit
import maws.space as space
from maws.complex import Complex
from maws.rna_structure import load_rna_structure
from maws.routines import entropy_score

cpx = Complex(
    force_field_aptamer="leaprc.RNA.OL3",
    force_field_ligand="leaprc.protein.ff19SB",
    salt_conc=0.15,
)
cpx.add_chain("", load_rna_structure())
cpx.add_chain_from_pdb(                        # force fields passed again
    pdb_path="ligand.pdb",
    force_field_aptamer="leaprc.RNA.OL3",
    force_field_ligand="leaprc.protein.ff19SB",
    parameterized=True,
)

ligand_only = Complex(                         # a second Complex, just to get a centre
    force_field_aptamer="leaprc.RNA.OL3",
    force_field_ligand="leaprc.protein.ff19SB",
    salt_conc=0.15,
)
ligand_only.add_chain_from_pdb(
    pdb_path="ligand.pdb",
    force_field_aptamer="leaprc.RNA.OL3",
    force_field_ligand="leaprc.protein.ff19SB",
    parameterized=True,
)
ligand_only.build()

sampler = space.make_sampler(ligand_only, reach=10.0, probe=1.4)
rotations = space.NAngles(4)

cx = copy.deepcopy(cpx)                        # deepcopy of an object holding a Simulation
aptamer = cx.aptamer_chain()                   # chain 0 is the aptamer, by convention
aptamer.create_sequence("G")
cx.build()

positions0 = cx.positions[:]                   # save coordinates by hand
energies = []
for _ in range(5000):
    pose = sampler.generator()                 # not a generator; returns one sample
    cx.translate_global(aptamer.element, pose.position * unit.angstrom)
    cx.rotate_global(aptamer.element, pose.axis * unit.angstrom, pose.angle)
    for j in range(4):
        aptamer.rotate_in_residue(0, j, rotations.generator()[j])
    energies.append(cx.get_energy()[0])
    cx.positions = positions0[:]               # restore by hand, every loop

score = entropy_score(energies, beta=0.01)
```

### Proposed

```python
from maws import Assembly, ForceField, ResidueLibrary, build, entropy_score
from maws.sampling import SurfaceSampler

ff = ForceField.for_target("RNA", "protein", salt_conc=0.15)   # one value, one place
system = build(
    Assembly()
    .with_aptamer(ResidueLibrary.rna(), sequence="G")
    .with_ligand("ligand.pdb", ff),
    ff,
)

sampler = SurfaceSampler.around(system.chain("ligand"), reach=10.0, probe=1.4)
pose0 = system.pose                                            # a value; never copied

energies = [
    system.energy(
        pose0.place(system.chain("aptamer"), sampler.sample())
            .randomize_torsions(system.chain("aptamer").residue(0))
    )
    for _ in range(5000)
]

score = entropy_score(energies, beta=0.01)
```

Fewer lines is not the goal. Fewer things to know is. Gone: `deepcopy`, the second
`Complex`, manual save/restore, force fields typed twice, `* unit.angstrom` at call
sites, chain-by-index.

---

## Not changing

- **The science.** Bondi radii, SAS rejection, the entropy score, GB salt
  screening, and the grow-by-one-nucleotide strategy stay exactly as they are.
- **`pdb_cleaner.py`** (655 lines). Already clean and tested. It only moves.
- **`prepare.make_lib`.** Same behaviour, plus a caller-supplied residue name.
- **The `.maws_cache/` format.** Same SHA1-keyed `.prmtop`/`.inpcrd` pairs.

---

## 3 questions for you

1. **Is `rotate_in_residue` applying each torsion 3 times on purpose?** See
   [01-values.md](01-values.md#2-the-bug-this-makes-impossible). If it is a bug,
   fixing it changes the numbers every past run produced. Fix it before the
   redesign so the two effects stay separate.
2. **Delete `Complex.rigid_minimize`?** It raises `TypeError` on its first
   iteration and nothing calls it. See [04-energy.md](04-energy.md).
3. **Does anything outside this repo import `maws.complex.Complex`?**
   That sets how long the compatibility shim lives.

**Next:** open [01-values.md](01-values.md).
</content>
