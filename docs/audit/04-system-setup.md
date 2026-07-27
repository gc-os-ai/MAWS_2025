# Class D — Force field and system setup

Files: `maws/prepare.py`, `maws/complex.py`, `maws/tools.py`

---

<a id="d1"></a>
## D1 — CRITICAL (organic ligands only) — antechamber is never told the charge

**Where:** [prepare.py:127-146](../../maws/prepare.py#L127-L146)

```python
run([
    find_exe("antechamber"),
    "-i", ante_in, "-fi", ext,
    "-o", f"{name}.mol2", "-fo", "mol2",
    "-c", charges,            # "bcc" -> AM1-BCC
    "-rn", residue_name,
    "-at", ante_at,
], cwd=w)
```

No `-nc` flag. So antechamber assumes **net charge 0** for every ligand.

For a charged ligand, that gives wrong AM1-BCC partial charges across the **whole molecule** — not just the charged group. LEaP then gets the wrong total charge, which changes every electrostatic term and every GB solvation term in the run.

For aptamer design this is close to worst case. The aptamer backbone carries one negative charge per nucleotide, so target electrostatics dominate the interaction. A ligand that should be −1 modelled as neutral flips the sign of the leading term.

### The shipped default hits this

`maws2023.py:53` defaults to `--path ./data/pfoa.pdb`. That is perfluorooctanoic acid, which at any relevant pH is the perfluorooctanoate **anion** (−1).

Any run against the default ligand has wrong charges.

**Fix — 30 minutes.**

1. Add a `net_charge` parameter to `make_lib`
2. Thread it through `add_chain_from_pdb` to `MawsRunner` and the CLI
3. Pass it as `-nc`

Do not try to guess it. Require it for the organic path.

Also missing, same area: antechamber is never told the pH or protonation state, and nothing checks that the input PDB's hydrogens match the intended ionisation.

---

<a id="d2"></a>
## D2 — HIGH — AmberTools failures are not caught

**Where:** [tools.py:86](../../maws/tools.py#L86), [prepare.py:196](../../maws/prepare.py#L196), [complex.py:426-433](../../maws/complex.py#L426-L433)

```python
def run(cmd, cwd=None):
    subprocess.run(cmd, cwd=cwd, check=True)
```

`check=True` catches non-zero exits. But **tleap routinely exits 0 while printing fatal errors** — missing parameters, failed `check`, unrecognised residues.

Nothing parses stdout or `leap.log` for `FATAL`, `Errors =`, or `Could not find`.

Two gaps:

1. `make_lib` runs `check {residue}` and then throws the result away. It moves the `.lib` into place whatever it said. A `.lib` with missing parameters or zeroed charges moves cleanly.
2. `Complex.build` checks only that `.prmtop`/`.inpcrd` **exist**. Not that they describe the molecule you asked for.

During this audit, LEaP printed this and everything downstream ignored it:

```
teLeap: Warning! The unperturbed charge of the unit (-5.000000) is not zero.
Residues lacking connect0/connect1 - ... CASP 1, NGLU 1
```

The first is expected for a charged system. The second says two residues were built with no polymer connectivity. Neither reached the user.

**Fix — 1 hour.**

1. Capture stdout from `run()`
2. Scan for `Errors = [1-9]`, `FATAL`, `Could not find bond parameter` — raise on any hit
3. Assert the built atom count matches the sum of the chains' expected lengths

Step 3 alone would also have caught all of [Class E](05-input-preparation.md).

---

<a id="d3"></a>
## D3 — MEDIUM — `pert_min` moves the rigid target

**Where:** [complex.py:956-970](../../maws/complex.py#L956-L970) on `main`

```python
def pert_min(self, size=1e-1, iterations=50):
    for _repeat in range(iterations):
        for i in range(len(self.positions)):          # every atom, target included
            self.positions[i] += np.random.uniform(-size, size, 3) * unit.angstrom
        self.minimize()                                # relaxes every non-zero-mass particle
```

Called as `pert_min(size=0.5)` once per candidate per growth step (`run.py:307`, `maws2023.py:373`). It kicks and then relaxes the protein that MAWS treats as a rigid docking target.

The result flows into `best_positions → best_old_positions`, so the damage **accumulates across steps**. Each step inherits the previous step's damaged target and adds another 50 rounds.

### Already fixed on a branch

`fix/pert-min-perturbs-receptor`, commit `f89bd7b`. It adds an explicit `atoms` argument and zeroes the masses of everything outside it.

That commit's own measurement: **max target displacement 3.9 Å before, 1.8 × 10⁻¹⁵ Å after.**

Two things to carry over when you merge:

1. The fix must also land in `maws2023.py:373` — that CLI is what most runs go through. Check whether the branch updated it.
2. `minimize()` ([complex.py:876](../../maws/complex.py#L876)) uses `max_iterations=100`. That is not enough to relax a badly clashing junction. It returns silently whether or not it converged, and does not check for NaN.

---

<a id="d4"></a>
## D4 — MEDIUM — everything is named `LIG` and written to the current directory

**Where:** [complex.py:115-163](../../maws/complex.py#L115-L163) (`pdb_name: str = "LIG"`), `prepare.py:97, 199-201`

Every target — protein or organic — is parameterised as one residue named `LIG`, and `<input_dir>/LIG.lib` is written next to the input PDB.

Consequences:

- Two targets in the same directory **overwrite each other's `LIG.lib`**. The build cache hashes lib *contents* ([complex.py:306-353](../../maws/complex.py#L306-L353)), which correctly stops a stale topology being served. But it does not stop the file itself being clobbered mid-campaign.
- Concurrent runs in a shared data directory race on that file.
- `Complex.build` writes `out.in` and LEaP writes `leap.log` into the **current working directory** (`complex.py:409, 423`). `.maws_cache/` goes there too. The repo root already has `out.in` and `leap.log` from past runs. Start a run from a different directory and you silently get a cold cache.

**Fix — 1 hour.** Derive the residue name from the input (hash or stem). Put all artifacts under an explicit, configurable work directory.

---

<a id="d5"></a>
## D5 — LOW — dead and inconsistent `createSystem` arguments

**Where:** [complex.py:457-463](../../maws/complex.py#L457-L463)

```python
return self.prmtop.createSystem(
    nonbondedCutoff=5 * unit.angstrom,     # ignored: NoCutoff is set below
    nonbondedMethod=app.NoCutoff,
    constraints=None,
    implicitSolvent=app.OBC1,
    implicitSolventSaltConc=self.salt_conc * unit.molar,
)
```

Three notes:

1. `nonbondedCutoff` is silently ignored under `NoCutoff`. Harmless, but it reads as though a 5 Å cutoff is active — which would be badly wrong for a charged system in GB.
2. `constraints=None` plus the 2 fs Langevin timestep (`complex.py:444-446`) is unstable for unconstrained X–H bonds. Only matters for `Complex.step()`, which nothing on the main path calls. If you ever enable MD, use `constraints=app.HBonds` or drop to 1 fs.
3. `implicitSolvent=app.OBC1` uses whatever Born radii the prmtop carries (mbondi by default). OBC models are normally paired with **mbondi2**. Set it deliberately with `set default PBradii mbondi2` in the LEaP script.

The salt-screening term from #41 is a real improvement and is wired correctly — `salt_conc` flows from the CLI to `implicitSolventSaltConc`. Its docstring caveat that GB screening is monovalent-only and does not model Mg²⁺ is accurate and important. Most real aptamer folding is Mg²⁺-dependent.

---

## D6 — Note — `groups=1` is safe for now

`complex.py:857, 879, 907` pass `groups=1` to `getState`. That is a **bitmask selecting force group 0 only**, not "all groups".

Verified on a real system that every force is in group 0:

```
HarmonicBondForce 0, HarmonicAngleForce 0, PeriodicTorsionForce 0,
NonbondedForce 0, CustomGBForce 0, CMMotionRemover 0
```

So today the energies are complete.

It is fragile. Any force added to a different group later would be silently dropped from every energy MAWS computes, with no error.

**Fix — 2 minutes.** Use `getState(getEnergy=True)` with no `groups`, or `groups=-1`.
