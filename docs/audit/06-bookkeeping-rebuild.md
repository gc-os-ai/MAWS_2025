# Class F — Chain rebuild

Files: `maws/chain.py`, `maws/complex.py`, `maws/rna_structure.py`, `maws/dna_structure.py`

`Complex.rebuild` ([complex.py:467-689](../../maws/complex.py#L467-L689)) splices a new nucleotide onto the chain while keeping everyone else's coordinates. About 200 lines of index arithmetic with an inlined rotation matrix. It is the hardest code in the repo to read.

## Good news first

**The index bookkeeping is correct.** I traced `Chain.update_chains`, `append_sequence`, `prepend_sequence` and all three splice branches through a full append and a full prepend. `start`, `start_history`, `length` and `length_history` line up for both the aptamer and the ligand chain — including the awkward case where an empty aptamer chain and the ligand chain both have `start == 0`. The ligand's coordinates are carried from the right slice every time.

That machinery works. The problem is the chemistry in the connectivity table.

---

<a id="f1"></a>
## F1 — HIGH — the `CONNECT` table names the wrong atoms

**Where:** [rna_structure.py:155](../../maws/rna_structure.py#L155) and `dna_structure.py`

```python
CONNECT = [[[0, -1], [-2, 0], 1.6, 1.6] for _ in RESIDUE_NAMES]
```

Read as `[[append_new, append_old], [prepend_new, prepend_old], append_len, prepend_len]`, this says:

- appending: bond new residue's atom `0` to old residue's atom `-1`
- prepending: bond new residue's atom `-2` to next residue's atom `0`

The same triple is used for **all 16 residue types**. But the residue types do not share an atom layout.

Resolved against real LEaP atom ordering:

| residue | atoms | index −1 | index −2 | correct 3'-anchor |
|---|---|---|---|---|
| RNA `G5` | 32 | O3' ✓ | **HO2'** | O3' = −1 |
| RNA `G` | 34 | O3' ✓ | HO2' | O3' = −1 |
| RNA `GN` | 33 | **HO3'** | O3' ✓ | O3' = −2 |
| RNA `A3` | 34 | **HO3'** | O3' ✓ | O3' = −2 |
| DNA `DG5` | 31 | O3' ✓ | **H2''** | O3' = −1 |
| DNA `DGN` | 32 | **HO3'** | O3' ✓ | O3' = −2 |

`-1` is right for residues with a 5'-phosphate or a downstream neighbour. It is the terminal hydroxyl **hydrogen** for the `N` and `3` forms. `-2` is the mirror image.

One constant cannot be right for both.

### What goes wrong in `rebuild`

**Append** ([complex.py:662-666](../../maws/complex.py#L662-L666)):
`shift_back = chain_positions[-1]`. The old chain's last residue was an `XN` or `X3` form, so its last atom is **HO3'**, not O3'. The junction is anchored one atom off — about 0.96 Å.

**Prepend** ([complex.py:528-530](../../maws/complex.py#L528-L530), `:549-551`, `:581-586`):
`connect[prepend_history[-1]][1][0]` = `-2`, resolved against the new 5' residue (`G5`, `DG5`), gives **HO2'** for RNA and **H2''** for DNA. Those are ribose 2'-position atoms — several ångströms away from O3', on a different face of the sugar.

The prepend junction is anchored on an atom that has nothing to do with the bond being formed.

The hardcoded 1.6 Å bond length is fine. The Amber P–O3' bond is ~1.61 Å.

### Why nobody noticed

`run.py:307` calls `pert_min(size=0.5, iterations=50)` right after every `rebuild`. That is 50 rounds of kick-and-minimise, which drags a badly placed junction back toward something bonded.

So the damage shows up as slow, expensive relaxation rather than a crash. It also means `pert_min` is doing structural repair, not the conformational exploration its docstring describes — and on `main` it does that repair by moving the **target** too ([D3](04-system-setup.md#d3)).

### Fix — 3 hours

1. Make `CONNECT` per-residue-type instead of one shared triple
2. Derive it from atom names (`O3'`, `P`) resolved against the built topology, not from positional indices
3. Assert after every `rebuild` that the new inter-residue P–O3' distance is within 0.2 Å of 1.61 Å, **before** any minimisation runs

---

<a id="f2"></a>
## F2 — MEDIUM — the junction geometry cannot be reviewed as written

**Where:** [complex.py:524-599](../../maws/complex.py#L524-L599) (prepend), `:602-679` (append)

Beyond F1, this code cannot be checked with confidence and cannot be tested as it stands:

1. **The rotation matrix is inlined twice**, duplicating 20 lines already in `rotate_global` ([complex.py:788-806](../../maws/complex.py#L788-L806)). Three copies in one file is three places for a sign to drift.
2. **The sign convention is implicit.** `np.dot(vector, rot)` computes `Rᵀv` — the *inverse* of the matrix built just above it. Then `angle = -ang(...)` applies a second negation. They may cancel. But `helpers.angle` returns an **unsigned** angle in [0, π], so the rotation's sign is carried entirely by the `np.cross` axis direction, with no comment saying so. `helpers.directed_angle` returns a signed angle and is never used.
3. **`pre_positions[-1 + connect[...][1][0]]`** (`:550`) mixes a Python negative index with a connectivity offset in one expression. You cannot tell from the code whether that is deliberate or an off-by-one.
4. **The degenerate-case guard** `if all(axis == np.zeros(3))` (`:540`, `:620`) catches exactly-parallel vectors but not near-parallel ones. There, `axis /= np.linalg.norm(axis)` turns floating-point noise into an arbitrary rotation axis.

None of these is provably wrong on its own. But nothing tests them, and `pert_min` hides their output. The honest read: **this has never been verified and cannot be trusted until it is.**

### Fix — 30 minutes

Do not review this code. Bound it instead.

After every `rebuild`, before any minimisation, assert every covalent bond in the aptamer is within tolerance of its equilibrium length. About 10 lines using the prmtop's own bond list — see [`repro/01_rotation_geometry.py`](repro/01_rotation_geometry.py) for the pattern.

That turns this whole class from "unknown" into "known good" or "known broken".
