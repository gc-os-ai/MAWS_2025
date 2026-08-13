# Class E — PDB cleaning

File: `maws/pdb_cleaner.py`

## Does this affect you?

Only if you used `--clean-pdb` **and** `molecule_type == "protein"` ([pdb_cleaner.py:631-640](../../maws/pdb_cleaner.py#L631-L640)). It is off by default.

If you never used `--clean-pdb`, skip this file.

But the module exists because large proteins otherwise fail in LEaP. So the protein workflow is exactly where you were most likely to switch it on.

The cleaned file goes straight to `tleap loadpdb` via `prepare.make_lib(..., parameterized=True)`, which `saveoff`s the whole protein as one library entry. **Nothing between the cleaner and LEaP checks that the molecule survived.**

E1, E2 and E3 were verified directly against `data/1HAO.pdb` and synthetic inputs. The rest carry their verification status inline.

---

<a id="e1"></a>
## E1 — CRITICAL — one altLoc re-sorts the whole file

**Where:** [pdb_cleaner.py:551-557](../../maws/pdb_cleaner.py#L551-L557), `_resolve_altloc`

```python
df.sort_values(["Seq_Num", "ChainID", "AtomTyp", "__occ", "__rk"])
  .drop_duplicates(subset=["Seq_Num", "ChainID", "AtomTyp"], keep="first")
```

Three problems:

1. The sort hits the **whole DataFrame**, not just the altLoc rows
2. It is never undone before writing
3. `Seq_Num` is a **string**, so the order is alphabetical, not numeric

`TER` rows have an empty `AtomTyp`, so they get moved too.

**`tleap loadpdb` builds connectivity from file order.** Reordering the file rewrites the protein's covalent structure.

### Proof

A 4-residue chain (`ALA 8, GLY 9, SER 10, VAL 11`) with one altLoc pair on residue 10.

Input order: `8, 9, 10, 11`. Output:

```
ATOM  ...  C SER A  10        <- residue 10 first
ATOM  ...  CA SER A  10
ATOM  ...  N SER A  10
ATOM  ...  O SER A  10
TER    21      VAL A  11      <- TER before the residue it terminates
ATOM  ...  C VAL A  11
...
ATOM  ...  C ALA A   8        <- residues 8 and 9 last
...
ATOM  ...  O GLY A   9
```

Three corruptions in one file:

- residues emitted as `10, 11, 8, 9`
- `TER` placed **before** the residue it belongs to
- atoms alphabetised inside each residue (`C, CA, N, O` instead of `N, CA, C, O`)

Most X-ray structures at 2.0 Å or worse have at least one altLoc. So this fires on a typical protein target. `data/1HAO.pdb` has **zero** altLocs, which is why it does not reproduce there.

**Fix — 30 minutes.** Save the original row order, sort a copy for the dedup, reindex back.

---

<a id="e2"></a>
## E2 — CRITICAL — insertion-code residues are deleted, including their TER

**Where:** [pdb_cleaner.py:594-595](../../maws/pdb_cleaner.py#L594-L595)

```python
if pdb_df.InsCode.fillna(" ").nunique() > 1:
    pdb_df = pdb_df[pdb_df.InsCode.fillna(" ") == " "]
```

Every atom with an insertion code is dropped.

Insertion codes are how Kabat/Chothia-numbered antibodies (and many structures with inserted loops) express residues that sit *between* two sequential numbers. Deleting them removes residues from the **middle** of a chain. The residues either side become adjacent in the file, and LEaP bonds them directly.

Worse: if a chain's `TER` record carries an insertion code, that `TER` is deleted too. LEaP then bonds two separate chains into one molecule.

### Proof — on the repo's own data

`data/1HAO.pdb` (thrombin + 15mer):

| | original | cleaned |
|---|---|---|
| atoms | 2769 | **2281** (488 gone, 17.6%) |
| TER records | 3 | **2** (two chains fused) |
| chain D | 334 | 315 |
| chain L | 233 | **113** (52% of the chain gone) |
| chain H | 2202 | 1853 |
| insertion-coded atoms | 36 | 0 |

Chain L loses more than half its atoms. The `TER` separating it from chain H is gone. So the file handed to LEaP says chain L is covalently bonded to chain H.

No warning is printed.

**Fix — 2 hours.** Renumber insertion-coded residues into the sequence instead of deleting them. Never drop a `TER`.

---

<a id="e3"></a>
## E3 — CRITICAL — everything after the last TER is thrown away

**Where:** [pdb_cleaner.py:189-191](../../maws/pdb_cleaner.py#L189-L191) and `:582`

```python
last_ter = ter_idx[-1]
pdb_main = pdb_df.iloc[: last_ter + 1, :].copy()
ligand   = pdb_df.iloc[last_ter + 1 :, :].copy()
...
_method, pdb_df, _ligand = pdb_reader(input_path)     # _ligand is never used again
```

Two failure modes:

1. **Bound ligands, cofactors, metals and waters are dropped even with `drop_hetatm=False`.** Verified on `data/1HAO.pdb`: the 30-atom `0G6` inhibitor and 149 waters are gone from the cleaned file. For a thrombin structure the co-crystallised inhibitor is arguably the most interesting part of the target.
2. **A file whose last chain has no `TER` loses that whole chain.** Common from modelling and minimisation tools. No warning.

**Fix — 30 minutes.** Do not split at the last TER. Keep everything and filter by record type.

---

<a id="e4"></a>
## E4 — HIGH — the residue whitelist deletes MSE, Amber variants and caps

**Where:** [pdb_cleaner.py:61-108](../../maws/pdb_cleaner.py#L61-L108), applied at `:600-605`

Any residue not in the 20 standard 3-letter names or a small nucleic set is deleted whole. That includes:

- Amber protonation and disulfide variants: `HIE`, `HID`, `HIP`, `CYX`, `ASH`, `GLH`, `LYN`
- `MSE` (selenomethionine) — very common in crystal structures
- terminal caps `ACE`, `NME`
- Amber nucleic naming: `RA`, `RU`, `RG5`, `DA5`, `DG3`

*Verified by running:* a chain of `HIE/CYX/ALA/MSE/HID` reduces to `ALA` alone. A 5-residue Amber-named nucleic chain produces a cleaned file containing only `END`.

Deleting `MSE` punches holes in the middle of chains — same downstream damage as E2. Deleting `CYX` breaks disulfides.

**Fix — 1 hour.** Rename `MSE` to `MET` (standard practice). Whitelist the Amber variants.

---

<a id="e5"></a>
## E5 — HIGH — nothing checks the output, and failures fall back silently

**Where:** [pdb_cleaner.py:607-611](../../maws/pdb_cleaner.py#L607-L611) and `:642-655`

```python
except Exception as e:
    if logger is not None:
        logger.warning("PDB cleaning failed; using original file. Reason: %s", e)
    return pdb_path, original
```

Three ways this bites:

1. An empty or near-empty cleaned PDB (E4) goes to LEaP with no error. The run proceeds against a target that is not there.
2. If cleaning raises — say the output is written next to the input (`:608-609`) and the data directory is read-only — the run continues on the **uncleaned** file, at `warning` level, while you believe cleaning happened.
3. Nothing anywhere compares atom counts before and after.

**Fix — 20 minutes.** Assert the cleaned atom count is non-zero and within tolerance of the original. That one check would have caught E2, E3 and E4 the first time any of them fired.

---

<a id="e6"></a>
## E6 — HIGH — multi-model NMR files are stacked, not split

**Where:** [pdb_cleaner.py:167](../../maws/pdb_cleaner.py#L167)

```python
if line.startswith(("ATOM", "HETATM", "TER")):
```

`MODEL` and `ENDMDL` are ignored. And `pdb_reader` cuts at the *last* `TER`, which is the end of the last model (E3).

A 20-model NMR target comes out as 20 superimposed copies with duplicate atom serials and duplicate residue numbers. LEaP builds a 20× system with atoms on top of each other.

*Verified by running* on a synthetic 2-model file: both models present in the output, both `TER` lines, serials 1–8 twice.

**Fix — 20 minutes.** Take model 1 only, unless told otherwise.

---

<a id="e7"></a>
## E7 — MEDIUM — `drop_hetatm=False` does not keep HETATM

**Where:** `pdb_cleaner.py:588-589` vs `:600-605`. The docstring at `:574` promises the opposite.

Even with `drop_hetatm=False`, structural Zn/Mg/Ca ions, heme and NAD are removed by the non-standard-residue filter (E4) — and anything after the last `TER` by E3.

Removing a structural metal changes the charge state and the fold the force field sees.

*Verified by running:* a synthetic file with `HETATM ZN` and `HETATM FE/HEM` placed before the `TER`, cleaned with `drop_hetatm=False` — both gone.

---

<a id="e8"></a>
## E8 — MEDIUM — altLoc resolution ignores the residue name and works per atom

**Where:** [pdb_cleaner.py:555](../../maws/pdb_cleaner.py#L555): `drop_duplicates(subset=["Seq_Num", "ChainID", "AtomTyp"])`

Two failure modes:

1. **Two different residue types at one position.** The winner is picked independently per atom. Shared backbone atoms take one residue name; side-chain atoms unique to the loser keep the other. One `resSeq` ends up with two residue names. *Verified:* `SER` (altLoc A, occ 0.40) / `ALA` (altLoc B, occ 0.60) at residue 10 gives five `ALA A 10` atoms plus one `OG` labelled `SER A 10`.
2. **Non-uniform occupancies.** Occupancy is compared per atom, so a residue with different refined occupancies gets a mix of conformer A and conformer B atoms. That geometry cannot physically exist — typically the side chain clashes into itself.

**Fix — 30 minutes.** Resolve altLocs per **residue**. Pick one conformer for the whole residue.

---

<a id="e9"></a>
## E9 — LOW — five smaller issues

All read-only verification unless noted.

1. **`:464` — the altLoc letter is never cleared.** The output still declares altLoc `A`/`B` in column 17, despite claiming they were resolved. Harmless for LEaP, misleading for anything else. *Seen in the E1 output above.*
2. **`:518, :533` — chain IDs are upper-cased.** PDB chain IDs are case-sensitive. `keep_chains="A"` also keeps chain `a`. In assemblies that use both, you silently get two chains.
3. **`:522-529` — `keep="one"` picks the chain with the most *atoms*, not the longest chain.** A short chain surrounded by retained heteroatoms can beat the actual longest polypeptide. It also silently discards every other chain.
4. **`:591` — `remove_h` keys only on the element column.** Legacy PDBs with blank columns 77–78 keep all their hydrogens (the cleaner never infers the element). Deuterium is never removed. `check_hydrogen` at `:291-300` has a broader test but is never called.
5. **`:502-503` — the fallback `TER` reuses the last atom's serial** instead of last + 1. *Verified:* `TER 8` follows `ATOM 8`. Also `:368-377`: `_atom_field` ignores the `element` argument it is given and right-aligns atom names, breaking the column-13 convention.

---

## Not a bug

Skipping cleaning for `molecule_type in {"organic", "lipid"}` (`:634-640`) is the right call. E4 would delete an organic ligand outright.

The gap it leaves: there is then **no** path that strips crystallographic waters from an organic-ligand PDB before antechamber sees it. Antechamber will happily try to type them.
