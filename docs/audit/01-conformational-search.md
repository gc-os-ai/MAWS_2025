# Class A — Rotation kernel

Files: `maws/chain.py`, `maws/complex.py`, `maws/dna_structure.py`

This is the code that turns a random angle into a moved atom. It runs about a million times per MAWS run. It has zero test coverage — see [Class G](07-verification-gap.md).

---

<a id="a1"></a>
## A1 — CRITICAL — every torsion is applied 3 times, with broken indices

**Where:** [chain.py:346-372](../../maws/chain.py#L346-L372), `Chain.rotate_in_residue`

### What the code does

```python
for i in range(len(element)):
    if element[i] and element[i] < 0:
        element[i] += self.structure.residue_length[...]

    if element[2] is None:
        revised_element = [...]
    elif element[2] == 0:
        revised_element = [...]
    else:
        revised_element = [...]
        rev = False

    self.rotate_element(revised_element, angle, reverse=rev)   # <-- inside the loop
```

Two problems in one loop:

1. `element` has 3 items, so the rotation runs **3 times** instead of once.
2. On pass `i`, only indices `0..i` have been made positive. The first pass or two build the rotation from indices that are **still negative**.

Negative indices go into `for j in range(element[1], element[2])`. `range(-7, 10)` gives `-7, -6, ... 9`. In Python, index `-7` means "7th from the end". So the code grabs the **last 7 atoms of the chain plus the first 10** and rotates them as one lump.

### Which torsions are affected

2 of the 4 torsions MAWS samples have negative index specs. So this is not an edge case.

| torsion | spec (RNA `G5`) | ranges actually rotated, in order |
|---|---|---|
| j=0 (5'-OH) | `[0, 1, None]` | `[1:66)` `[1:66)` `[1:66)` |
| j=1 (beta) | `[1, 2, None]` | `[2:66)` `[2:66)` `[2:66)` |
| **j=2 (chi)** | `[10, 12, -7]` | **`[-7:10)`** **`[-7:10)`** `[10:25)` |
| **j=3 (epsilon)** | `[-7, -1, None]` | **`[-1:66)`** `[31:66)` `[31:66)` |

For chi, the axis is also wrong on the first two passes: it uses `positions[-1] - positions[25]` instead of the C1'–N9 bond.

### Proof

Built a real RNA `G5–A3` dinucleotide with LEaP. 71 covalent bonds. Applied one 0.7 rad torsion:

| call | max bond error | bonds broken (>0.05 Å) |
|---|---|---|
| `rotate_in_residue(1, 0, θ)` — alpha | 0.000 Å | 0 / 71 |
| `rotate_in_residue(1, 1, θ)` — beta | 0.000 Å | 0 / 71 |
| **`rotate_in_residue(1, 2, θ)` — chi** | **8.484 Å** | **5 / 71** |
| **`rotate_in_residue(0, 3, θ)` — epsilon** | **2.633 Å** | **1 / 71** |

Angle check on a torsion that does not break: requested 17.19°, got −51.57°. Requested 40.11°, got −120.32°. Ratio exactly −3.000.

### Why it stayed hidden

`NAngles.generator` draws θ uniformly from [0, 2π). And `3θ mod 2π` is also uniform. So for j=0 and j=1 the *set* of conformations sampled is unchanged. Nothing looked wrong.

That protection does not apply to j=2 and j=3. There, the damage is structural, not angular.

### What it costs you

Every chi and epsilon sample — half the torsional search — measures the energy of a molecule with broken bonds. Those energies land in the 10⁴–10⁸ kJ/mol range. Then they feed the score, where per [C1](03-selection-criterion.md#c1) they help the candidate **win**.

### Fix — 30 minutes

Normalise all 3 indices first. Then rotate once.

```python
element = [
    (x + L if (x is not None and x < 0) else x)
    for x in self.structure.rotating_elements[resname][residue_element_index]
]
revised_element = [...]        # build once
self.rotate_element(revised_element, angle, reverse=rev)   # call once
```

Write the bond-length test first, or you will not know it worked.

---

<a id="a2"></a>
## A2 — CRITICAL — `reverse=True` pivots on the wrong atom

**Where:** [chain.py:283-304](../../maws/chain.py#L283-L304) and [complex.py:715-724](../../maws/complex.py#L715-L724)

**Used by:** every prepend move — `run.py:315-327`, `maws2023.py:381-393`

### What it should do

Rotate the *other* side of a torsion. You need this when you grow the chain at the 5' end: move the new residue, leave the rest alone.

### What it actually does

It encodes "the other side" in the index numbers instead of moving the pivot point.

```python
# Chain.rotate_element
if rev:
    if revised_element[2] is None:
        revised_element[2] = 0        # end := chain-local 0
    else:
        revised_element[2] = revised_element[1]
rev = False                            # never reaches Complex
```

Then `Complex.rotate_element` sees `element[2] <= element[0]` and swaps:

```python
vec_a = pos[element[1]] - pos[element[0]]      # axis: correct, this is the bond
if element[2] <= element[0]:
    element[1], element[2] = element[2], element[1]
self.rotate_global(element, vec_a, angle, reverse=False, glob=False)
```

`rotate_global` takes its pivot from `element[1]`. The swap just set that to **chain atom 0**.

So: the axis points the right way, but it passes through the chain's first atom (HO5') instead of through either atom of the bond.

Rotating around an axis that misses the bond moves one bonded atom and not the other. The bond stretches.

### Proof

Same dinucleotide, θ = 0.7 rad, `reverse=True` on residue 0:

| torsion | max bond error | bonds broken |
|---|---|---|
| j=0 | 0.000 Å | 0 / 71 — degenerates to a no-op |
| j=1 (beta) | **0.351 Å** | 1 / 71 |
| j=2 (chi) | 4.836 Å | 6 / 71 — this one is A1, not A2 |
| j=3 (epsilon) | **0.299 Å** | 1 / 71 |

0.3 Å is small next to A1's 8.5 Å. But it happens on **every** prepend sample, and 0.35 Å on a 1.4 Å bond is roughly 300 kJ/mol of fake strain energy. That is the same size as the differences the search is trying to measure.

Also: **j=0 with `reverse=True` does nothing at all.** The swap makes the rotated range `[0:1)` — one atom, which is also the pivot. One of the four prepend degrees of freedom is dead.

### Fix — 2 hours

Stop encoding "complement" in the indices. Pass `reverse` through to `rotate_global`. Have it rotate `[chain_start, element[0])` and pivot on `element[0]`, the bond atom.

That is what the docstring at `complex.py:754-756` already claims happens.

---

<a id="a3"></a>
## A3 — HIGH — only 3 of 6 backbone degrees of freedom are sampled

**Where:** `run.py:246-247, 315-327`; tables in `rna_structure.py:70-151`, `dna_structure.py`

The code says `N_BACKBONE_TORSIONS = 4` and calls them "4 backbone torsions per residue".

Here is what those 4 actually are, checked against real LEaP atom names:

| j | bond rotated | name | moves the new residue? |
|---|---|---|---|
| 3 (on residue −2) | C3'–O3' | epsilon | yes |
| 0 | P–O5' | alpha | yes |
| 1 | O5'–C5' | beta | yes |
| 2 | C1'–N9/N1 | chi | **no** — only spins the base |

Six torsions set the geometry between two nucleotides: epsilon, **zeta** (O3'–P), alpha, beta, **gamma** (C5'–C4'), delta (ribose pucker).

MAWS samples epsilon, alpha, beta. **Zeta and gamma are never touched.** And the fourth one (chi) contributes nothing to placement.

So the new nucleotide is placed with 3 degrees of freedom, not 5. Whole regions of backbone space cannot be reached, no matter how many samples you draw. 5000 samples covers a 3-torus densely. It cannot fill in a missing dimension.

Step 1 is worse. For a single nucleotide:

- j=3 rotates one hydroxyl hydrogen (the range starts at O3', which is also the pivot)
- j=0 rotates the whole residue around the HO5'–O5' axis — which duplicates the global orientation draw that just happened

Step 1 has 2 meaningful internal degrees of freedom, not 4.

This is a design limit, not a slip. But it is undocumented, and it caps what the method can find.

---

<a id="a4"></a>
## A4 — MEDIUM — DNA terminal residues rotate around a non-bond

**Where:** `dna_structure.py`, 4th `ROTATIONS` entry for `DGN`/`DAN`/`DTN`/`DCN` and `DG3`/`DA3`/`DT3`/`DC3`: `(-4, -2, None)`

Checked against real LEaP atom ordering:

| residue | atoms | index −4 | index −2 | should be |
|---|---|---|---|---|
| `DGN` | 32 | **H2'** | O3' | C3'(25) → O3'(30), i.e. `(-7, -2)` |
| `DA3` | 33 | **H2'** | O3' | C3'(26) → O3'(31), i.e. `(-7, -2)` |

The axis runs between a hydrogen and an oxygen that are not bonded to each other.

RNA gets this right (`(-8, -2)`). So do internal DNA (`DG`) and 5'-terminal DNA (`DG5`). This is a DNA-only transcription error.

**Why it is Medium, not Critical:**

- `DX3` never gets a j=3 rotation on the main path. Append applies j=3 to residue −2, which is never a 3' terminus.
- `DXN` does get one, in step 1. But the range starts at O3', which is also the pivot, so the net effect is spinning HO3' around a bogus axis. Wrong hydroxyl geometry, not a broken backbone.

Still fix it. It means DNA and RNA runs are not doing the same thing, and it would bite hard if anyone changed the torsion indexing in `run.py`.

**Fix — 5 minutes.** Change `(-4, -2, None)` to `(-7, -2, None)` for the 8 affected residues.

---

<a id="a5"></a>
## A5 — MEDIUM — the shared template gets mutated

**Where:** [chain.py:341-350](../../maws/chain.py#L341-L350)

```python
element = self.structure.rotating_elements[...][residue_element_index]
for i in range(len(element)):
    if element[i] and element[i] < 0:
        element[i] += self.structure.residue_length[...]      # writes in place
```

`element` is a live reference into `Structure.rotating_elements`. Not a copy. The first call permanently rewrites the module-level template.

It happens to be safe today: once positive, the branch never fires again, and residue length is keyed by name. But:

- The `Structure` object's state now depends on which rotations you have requested
- `Structure.torsions()` returns different values before and after a search runs
- Any future change to length resolution turns this into an order-dependent bug

**Fix — 2 minutes.** Copy the triple before normalising.

---

<a id="a6"></a>
## A6 — MEDIUM — `rigid_minimize` crashes, and its logic is backwards

**Where:** [complex.py:932-954](../../maws/complex.py#L932-L954)

```python
energy = None
...
free_E = self.get_energy()[0]
if free_E < energy or energy is None:      # (1)
    energy = free_E
    self.positions = positions[:]          # (2)
```

**Problem 1:** Python evaluates `free_E < energy` before the `or`. On the first pass `energy` is `None`, so it raises:

```
TypeError: '<' not supported between instances of 'float' and 'NoneType'
```

Verified on a real complex. The operands are backwards — it should be `if energy is None or free_E < energy`.

**Problem 2:** `positions` is the snapshot from *before* the trial rotation.

- Move improves the energy → the code saves the improved energy, then **throws the improved coordinates away**
- Move makes it worse → the worse coordinates are kept

Accept and reject are swapped.

Nothing calls `rigid_minimize`, so no results are affected. Delete it or fix it. As written it is a trap.

---

<a id="a7"></a>
## A7 — LOW — `rotate_element` silently rewrites its own arguments

**Where:** [complex.py:720-723](../../maws/complex.py#L720-L723)

```python
if revised_element[2] <= revised_element[0]:
    revised_element[1], revised_element[2] = revised_element[2], revised_element[1]
```

There is no such thing as a bad input here. Any triple whose end index falls at or below its start gets silently re-read as a different rotation.

This is the mechanism A2 rides on. It is also what turns A1's negative indices into large wrong rotations instead of an `IndexError` you would have noticed.

**Fix — 10 minutes.** Validate the triple and raise.
