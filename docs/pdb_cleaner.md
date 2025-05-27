### Reading `special_PDB_files.txt` — a quick reference

When you run **`pdb_cleaner.py`** it drops one text file in the same directory:
`special_PDB_files.txt`.

- Think of it as a checklist: every section is a potential _problem class_, and under each header the cleaner lists **only** the files that triggered that flag.

---

| Section header                 | Meaning (what was detected)                                                                                                                                 | Why it matters                                                                           | Typical follow-up                                                             |
| ------------------------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------- |
| `The files below have ligands` | The PDB contains **HETATM** records after the last `TER` of each chain (e.g. cofactors, ions, waters). The cleaner prints their three-letter residue codes. | Ligands are **not** parameterised in the default AMBER/OpenMM force-fields used by MAWS. | Either delete them or supply custom `.xml`/`.frcmod` parameters.              |
| `alternate locations`          | One or more atoms have alt-loc IDs (“A”, “B”…).                                                                                                             | Ambiguous coordinates; force-field sees overlapping atoms.                               | Keep only the first alt-loc or merge them.                                    |
| `non-standard residues`        | Amino acids outside the canonical 20 (e.g. MSE, HIP) or other polymers.                                                                                     | Missing parameters break system construction.                                            | Strip or patch them; alternatively load a force-field that supports them.     |
| `hydrogen atoms`               | Explicit H atoms remain in the file.                                                                                                                        | Mixed protonation and duplicate hydrogens after OpenMM’s `addHydrogens()`.               | Usually safe to let the script delete them; OpenMM will rebuild H later.      |
| `sequence gaps`                | Residue numbering jumps by > 1 (`…148, 150…`). The tuple shows the last residue _before_ and the first residue _after_ every gap.                           | Broken backbone; leads to large forces or NaN energies.                                  | Rebuild missing loops with **PDBFixer/Modeller** or truncate the region.      |
| `insertion code`               | Residues have insertion letters (148A, 148B…). All unique codes are listed.                                                                                 | Creates duplicate `(chain, resSeq)` keys; MAWS grouping logic fails.                     | Renumber residues to a clean monotonic series.                                |
| `multiple chains`              | More than one chain ID in the file.                                                                                                                         | MAWS assumes a single receptor chain unless you explicitly keep several.                 | Decide whether to keep all (`--keep all`) or only the longest (`--keep one`). |
| `negative sequence number`     | Any residue index < 0.                                                                                                                                      | Rare legacy numbering; some tools treat them as invalid.                                 | Renumber starting from 1.                                                     |

---

#### Example snippet explained

```text
-------------------------------------------------------------------------------
 The files below have ligands
examples/1hut.pdb    ligands: 0G7   HOH
```

- **`examples/1hut.pdb`** – path to the offending file
- **`0G7   HOH`** – 0G7 is a GDP analogue; HOH = water.
  → remove or parameterise before running MAWS.

```text
-------------------------------------------------------------------------------
 The files below have sequence gaps
examples/1hut.pdb    sequence gaps:('H:148', 'H:150')('H:217', 'H:219')
```

- Chain **H** jumps 148 → 150 and 217 → 219.
  → rebuild or trim these loop regions.

---

### Workflow after reading the report

1. **Open each flagged file in PyMOL/ChimeraX**
   Quickly inspect ligands, gaps, alternate-locs.
2. **Choose a remediation strategy** per row of the table above.
   (Most users simply rerun the cleaner with `--keep one --no-hydrogens` then fix gaps.)
3. **Re-run the cleaner or external tools** until the report prints _empty_ sections.
4. **Feed the cleaned PDB to MAWS** (`maws run --pdb your_clean.pdb …`).

Once every section in `special_PDB_files.txt` is blank, your PDBs are in a state that OpenMM-based pipelines—including MAWS—can load without surprises.
