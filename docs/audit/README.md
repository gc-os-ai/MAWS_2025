# MAWS code audit

## Fix these first

1. **Add a bond-length test.** After any rotation, check every covalent bond is within 0.05 Å of what it was. ~15 lines. Start from [`repro/01_rotation_geometry.py`](repro/01_rotation_geometry.py). **30 minutes.**
2. **Fix `Chain.rotate_in_residue`** ([chain.py:346-372](../../maws/chain.py#L346-L372)). Move the rotation call outside the normalise loop. **30 minutes.** See [A1](01-conformational-search.md#a1).
3. **Fix `reverse=True`** ([chain.py:283-304](../../maws/chain.py#L283-L304)). It pivots on the wrong atom. **2 hours.** See [A2](01-conformational-search.md#a2).
4. **Decide what the score should rank.** The current one rewards steric clashes. This is a design call, not a typo. **Half a day to think, then a day to build.** See [C1](03-selection-criterion.md#c1).
5. **Fix pose placement** ([run.py:243](../../maws/run.py#L243)). The sampled point is used as a displacement, not a target. **1 hour.** See [B1](02-pose-sampling.md#b1).

Do them in this order. Step 1 first, because you cannot tell if steps 2 and 3 worked without it.

---

## What this is

An adversarial review of MAWS. I read the code and ignored the docs, as you asked.

- **Date:** 2026-07-27
- **Reviewed:** `main` @ `4126531`, checked out as `refactor/scipy-entropy-score`
- **Also checked:** `fix/pert-min-perturbs-receptor` (fixes D3 only)
- **Verified how:** every finding marked *verified* was reproduced by running code in the `maws` conda env against real AmberTools output. No production run was started.

Line numbers point at the working tree I audited. `maws/run.py` had uncommitted edits, so its numbers may drift by a few lines. Every reference also names the function, so it stays findable. `chain.py`, `complex.py`, `space.py`, `prepare.py` and `pdb_cleaner.py` are identical across `main` and both branches — that is where all the Critical findings are.

---

## The three big ones

**1. The rotation code breaks chemical bonds.**

Every torsion is applied 3 times instead of 1. For 2 of the 4 torsions, the first 2 applications use array indices that are still negative, which wrap around to the wrong end of the coordinate array.

One chi rotation moves a bond by 8.5 Å and breaks 5 of 71 bonds. Every sample. → [A1](01-conformational-search.md#a1)

**2. The steric filter does not filter.**

The sampler picks a point. The code then uses that point as a *shift* instead of a *destination*. And the surface check tests the point, not the 30 atoms that actually get placed there.

Result on 1BRQ: 16% of "clash-free" poses have an atom within 2.6 Å of the protein. 7% within 1.0 Å. → [B1](02-pose-sampling.md#b1), [B2](02-pose-sampling.md#b2)

**3. The score picks the candidate that clashes most.**

At `beta = 0.01` the entropy score works out to `log(fraction of samples that did not clash)`. Add pure garbage to a candidate and it wins. → [C1](03-selection-criterion.md#c1)

---

## Can you trust past results?

**The sequence: no.** It is picked only by the entropy score. `free_E` is computed and logged but never used to choose anything.

**The final pose: partly.** It is the lowest-energy pose using the real force field, so there is some real signal. But it minimises over broken geometry, starts from an offset placement, and on `main` the protein itself drifts.

### Check your old logs — 20 minutes, no re-running

Bug C1 makes a specific prediction you can test. Clashes go up with residue size. Atom counts:

- RNA: G 33 > A 32 > C 30 > U 29
- DNA: DG 32 > DA 31 = DT 31 > DC 29

So MAWS should prefer **G, then A**, and avoid **U/T and C** — no matter what the target is.

Your `*_entropy.log` files have one line per candidate:

```
SEQUENCE: <seq> ENTROPY: <score> ENERGY: <min energy>
```

Three checks:

1. Count which base got added at each step, across all past runs. G/A heavy on unrelated targets = the bug, not the target.
2. Plot ENTROPY against ENERGY. They should correlate negatively. If they do not, the score was not tracking energy.
3. Look at the gap between the winner and runner-up in each step. Small gaps = noise.

---

## What is actually fine

Worth saying, since most of this report is bad news:

- The RNA torsion tables in `rna_structure.py` are correct. I checked them against real LEaP atom ordering.
- The GB salt-screening work (#41) fixes a real physical problem and is wired up correctly.
- The scipy entropy rewrite matches the old mpmath one to 1.1 × 10⁻¹⁶. Faithful refactor. It just faithfully reimplements a broken criterion.

---

## Why nobody caught these

`tests/test_chain.py:31` replaces the rotation function with a no-op:

```python
def rotate_element(self, element, angle, reverse=False):
    """Mock rotation - does nothing in unit tests."""
    pass
```

So every rotation test passes against nothing. No test anywhere checks a bond length. → [Class G](07-verification-gap.md)

---

## All findings

| Class | Area | Critical | High | Medium | Low |
|---|---|---|---|---|---|
| [A](01-conformational-search.md) | Rotation kernel | 2 | 1 | 3 | 1 |
| [B](02-pose-sampling.md) | Pose sampling | 2 | 2 | 1 | – |
| [C](03-selection-criterion.md) | Selection score | 2 | 1 | 2 | – |
| [D](04-system-setup.md) | Force field setup | 1 | 1 | 2 | 2 |
| [E](05-input-preparation.md) | PDB cleaning | 3 | 3 | 2 | 3 |
| [F](06-bookkeeping-rebuild.md) | Chain rebuild | – | 1 | 1 | – |
| [G](07-verification-gap.md) | Missing tests | – | 1 | – | – |

Severity:

- **Critical** — invalidates results
- **High** — biases results a lot
- **Medium** — wrong, but bounded or not on the main path
- **Low** — wrong, currently harmless

---

## Reproduce any number in this report

```
conda activate maws
python docs/audit/repro/01_rotation_geometry.py     # A1, A2, A4   (~10 s)
python docs/audit/repro/02_pose_placement.py        # B1, B2       (~1 min, needs data/LIG.lib)
python docs/audit/repro/03_entropy_criterion.py     # C1, C2, C3   (instant)
python docs/audit/repro/04_pdb_cleaner.py           # E1, E2, E3   (~5 s)
```
