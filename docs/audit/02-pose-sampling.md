# Class B — Pose sampling and placement

Files: `maws/space.py`, `maws/run.py`, `maws/maws2023.py`

`space.py` itself is good code. The cube-root radial law, the Gaussian-normalise trick for uniform directions, the Bondi radii table, the KD-tree surface test — all correct.

The bug is at the seam. The sampler returns a **point**. `run.py` uses it as a **shift**.

All numbers below: `data/1BRQ_cleaned.pdb` (1408 atoms) as target, RNA `G` aptamer, `reach=10`, `probe=1.4` (the shipped defaults), 300 accepted poses. The repro script seeds the RNG so numbers are stable. Unseeded they move 2–3 points between runs — which is finding [B5](#b5) in miniature.

---

<a id="b1"></a>
## B1 — CRITICAL — the pose is used as a shift, not a destination

**Where:** [run.py:243-244](../../maws/run.py#L243-L244), `maws2023.py:293-294`

```python
pose = sampler.generator()
cx.translate_global(aptamer.element, pose.position * unit.angstrom)
cx.rotate_global(aptamer.element, pose.axis * unit.angstrom, pose.angle)
```

`sampler.generator()` returns an **absolute point** in the target's coordinate frame. `Sphere.centre` is the ligand's centre of mass (`space.py:296-298`) and `Sample.position` is `centre + r·direction`.

But `Complex.translate_global` ([complex.py:838-842](../../maws/complex.py#L838-L842)) just adds the vector to every atom:

```python
for j in range(element[0], element[2]):
    pos[j] += vec_shift
```

So the aptamer lands at `where_it_already_was + pose.position`. Not at `pose.position`.

### The correct call

```python
current_com = <aptamer centroid>
cx.translate_global(aptamer.element, (pose.position - current_com) * unit.angstrom)
```

### Proof

LEaP builds the aptamer near the origin, so the error is the aptamer's own distance from the origin — not something enormous:

```
aptamer centroid as LEaP builds it : [ 2.36 -0.13  3.95]   (4.62 Å from origin)
target centroid                    : [20.03 46.50 -35.82]
envelope: centre [20.1 46.6 -35.8], radius 38.97 Å

|aptamer centroid after move − sampled point|   mean 4.77 Å   min 1.92   max 6.69
```

So: a **systematic ~5 Å offset, always in the same direction**. About a tenth of draws land outside the envelope the sampler was asked to sample.

Not fatal alone — the aptamer does end up in roughly the right region. But it is why [B2](#b2) is as bad as it is. And it would grow without limit if LEaP ever stopped building near the origin, which is not a guaranteed behaviour.

**Fix — 1 hour.** Subtract the current centroid before translating.

---

<a id="b2"></a>
## B2 — CRITICAL — the surface filter tests the wrong thing

**Where:** [space.py:348-357](../../maws/space.py#L348-L357) and [space.py:257-264](../../maws/space.py#L257-L264)

```python
sample = self.envelope.generator()
if self.excluder.is_clear(sample.position):
    return sample
```

`is_clear` asks: could a single 1.4 Å water-sized sphere sit at this point without overlapping any target atom?

That is a correct surface test. But the thing being placed is not a water sphere. It is a nucleotide:

- ~30 atoms
- spanning ~10 Å
- placed ~4.6 Å away from the tested point (B1)
- then rotated around an arbitrary axis (B3)

Nothing checks the atoms that actually get placed.

### Proof

Over 300 accepted ("surface-clear") poses, measuring the real minimum distance from any aptamer atom to any target atom:

```
min aptamer-atom → target-atom distance:  mean 9.09 Å,  min 0.25 Å

contact < 1.0 Å : 6.7%  of ACCEPTED poses
contact < 2.0 Å : 13.3%
contact < 2.6 Å : 16.3%
contact < 3.0 Å : 18.0%
```

**16% of the poses the "surface-aware" sampler certifies as clear are hard clashes.** 7% have atoms basically on top of each other.

Energies confirm it. Over 60 samples the potential energy ranged from 1.49 × 10⁴ to **1.11 × 10⁸ kJ/mol**.

This is worse than a wasted-sample rate, because those clash samples are not discarded. Per [C1](03-selection-criterion.md#c1) they are exactly what makes a candidate win.

**Fix — 1 hour.** Reject based on the placed aptamer's atoms, after translating and rotating. One `KDTree.query(aptamer_positions, k=1).min()` against the target tree. Microseconds. It is the check the module docstring already describes.

---

<a id="b3"></a>
## B3 — HIGH — the orientation draw pivots on the first atom

**Where:** [run.py:244](../../maws/run.py#L244), [complex.py:782-785](../../maws/complex.py#L782-L785)

```python
cx.rotate_global(aptamer.element, pose.axis * unit.angstrom, pose.angle)
```

With `glob=True` and `reverse=False`, `rotate_global` pivots on `pos[element[0]]` — the aptamer's **first atom** (HO5'). Not its centre of mass.

Rotating a rigid body around a point that is not its centre also **moves** it, by up to twice the centre-to-pivot distance.

For a mononucleotide that is another ±5–9 Å of uncontrolled displacement. Applied *after* the placement in B1, and *after* the surface check in B2. This is a large part of why the clash rate is 16% and not a few percent.

**Fix — 20 minutes.** Pivot on the mass-weighted centre. `helpers.mass_weighted_center` already exists and is already used by `space.compute_envelope_dims`.

---

<a id="b4"></a>
## B4 — MEDIUM — the sampling region is not aimed at a binding site

**Where:** [space.py:274-298](../../maws/space.py#L274-L298)

```python
return {"radius": float(dists.max()) + reach, "centre": com}
```

The region is a sphere on the target's centre of mass, radius `furthest_atom + reach`.

Fine for a small organic ligand. Not fine for a protein. On 1BRQ the radius is **38.97 Å**, so the sphere holds the whole protein plus a 10 Å shell. After surface rejection, the accepted region is "anywhere in solvent within 10 Å of the protein" — the entire surface, evenly, with no preference for pockets.

That is a coverage problem. Spread 5000 samples over the full surface of a 1400-atom protein and each candidate site gets a handful. Nowhere near enough to tell four nucleotides apart on how well they fit a pocket.

`SurfaceFollowingSampler` ([space.py:360-453](../../maws/space.py#L360-L453)) is a better shape and looks correctly written. But `make_sampler` defaults to `mode="sphere"` and neither `run.py` nor `maws2023.py` lets you pick the other one — `mode` is never passed. There is also no way to give a known site centre.

**Fix — 2 hours.** Expose `mode` and an optional site centre through `MawsRunner` and the CLI.

---

<a id="b5"></a>
## B5 — HIGH — nothing is seeded, so no run can be reproduced

`grep -rn "seed\|RandomState\|default_rng" maws/` returns nothing.

Every random draw uses the global `numpy.random` state: `Sphere.generator`, `_random_unit_axis`, `NAngles.generator`, `pert_min`'s kicks, `rigid_minimize`.

No seed argument on `MawsRunner`. No `--seed` flag. No seed in the log.

Two consequences:

1. **No past result can be reproduced.** You cannot re-derive a sequence from its log.
2. **You cannot measure run-to-run variance.** That is the one diagnostic that would have exposed most of this report. If two runs on the same target give different sequences with similar scores, the ranking is noise — but with no seed you cannot tell "different because random" from "different because something changed".

**Fix — 1 hour.** Add a `seed` argument. Thread it to a `numpy.random.Generator` instance, not the global state. Log it. Put it in `MawsResult`.
