# Class C — The selection score

Files: `maws/routines.py`, `maws/run.py`, `maws/maws2023.py`

**This class is different.** Nothing here is a typo. The code does exactly what it says. The scipy rewrite is a faithful copy of the old mpmath one.

The problem is that the thing being computed does not rank binders. And at the shipped `beta`, it ranks something actively harmful.

---

> **Naming.** On `main` the function is `routines.S`. The `refactor/scipy-entropy-score` branch renames it to `entropy_score`. This file uses the new name. The repro script handles both.

## What the score is

`routines.entropy_score` returns `−D_KL(P ‖ U)`, where:

- `P` = the Boltzmann distribution over the sampled energies
- `U` = uniform over the same `N` samples

Which works out to:

```
score = H(P) − log N
```

It is always ≤ 0. It hits 0 when all energies are identical. It hits its minimum `−log N` when one sample carries all the weight.

`run.py:265` and `run.py:348` pick the candidate with the **lowest** score.

---

<a id="c1"></a>
## C1 — CRITICAL — the score rewards steric clashing

### The mechanism

Say `m` of your `N` samples are physically sensible, and `N − m` are hard clashes.

A clash energy of ~10⁸ kJ/mol gives `exp(−βE) = exp(−10⁶)`, which is **exactly zero** in double precision. Those samples carry no weight at all. They contribute nothing to `H(P)`.

If the surviving `m` samples have roughly comparable energies — which they do at `beta = 0.01` — then `H(P) ≈ log m`. So:

```
score ≈ log m − log N = log(m/N) = log(fraction that did NOT clash)
```

**More clashes → lower score → MAWS picks it.**

### Proof

Take one fixed set of 5000 good, clash-free energies (σ = 50 kJ/mol around −1000). Replace a fraction with hard clashes at +10⁸ kJ/mol. The good poses never change. The candidate only gets worse:

| clash fraction | score | predicted `log(1−f)` |
|---|---|---|
| 0.000 — the good candidate | −0.1227 | — |
| 0.100 | **−0.2293** | −0.1054 |
| 0.500 | **−0.8172** | −0.6931 |
| 0.900 | **−2.4281** | −2.3026 |
| 0.990 | **−4.7533** | −4.6052 |
| 0.998 | **−6.3049** | −6.2146 |

A candidate that clashes 99.8% of the time beats one that never clashes, by 6.2 points.

The measured numbers track `log(1−f)` closely. So this confirms the mechanism, not just the direction.

### This is not hypothetical

- [B2](02-pose-sampling.md#b2) measures a real 16% clash rate among *accepted* poses on 1BRQ
- [A1](01-conformational-search.md#a1) manufactures extra 10⁴–10⁸ kJ/mol samples on every chi and epsilon move

### The fingerprint this leaves

Clash chance goes up with the size of the nucleotide. Atom counts:

| | biggest → smallest |
|---|---|
| RNA | G (33) > A (32) > C (30) > U (29) |
| DNA | DG (32) > DA (31) = DT (31) > DC (29) |

So C1 predicts MAWS prefers **G, then A**, and avoids **U/T and C** — mostly regardless of the target.

That is falsifiable. And it is the cheapest way to check how much a past run was driven by this bug.

<a id="what-to-check-in-your-old-logs"></a>
### Check your old logs — 20 minutes, no re-running

`maws2023.py:326` and `:416` write one line per candidate to `{JOB}_entropy.log`:

```
SEQUENCE: <seq> ENTROPY: <score> ENERGY: <min energy>
```

Three checks:

1. **Composition.** Count which base got added at each step, across all past runs. If G/A dominate on chemically unrelated targets, the ranking was driven by residue size, not by the target.
2. **Score vs energy.** Plot ENTROPY against ENERGY per step. A working score correlates negatively (better binder → lower energy → sharper peak). Uncorrelated or positive means the score was not tracking energy.
3. **Winning margins.** Look at the spread of ENTROPY across the 8 candidates in a step. If the winner barely beats the runner-up, the choice was noise. (You cannot measure the noise floor retrospectively — see [B5](02-pose-sampling.md#b5).)

---

<a id="c2"></a>
## C2 — CRITICAL — the score cannot see binding strength

`P(i) = exp(−βE_i) / Σ exp(−βE_j)` does not change if you add a constant to every energy.

So the score depends **only on the spread** of energies. Never on their actual value.

A candidate whose poses all sit at −5000 kJ/mol and one whose poses all sit at +5000 kJ/mol get **identical** scores.

The score cannot tell a strong binder from a weak one, or from a non-binder. It measures how sharply peaked the ensemble is. That is a statement about entropy, not affinity.

There is a reasonable intuition behind it — a well-fitted pose should have a narrower energy spread. But that argues for using it as a **tie-breaker** next to an affinity term. Not as the only criterion.

One upside worth noting: this same shift-invariance saves you from a different problem. Candidates in a step have different atom counts (G and U differ by 4 atoms), so their absolute energies are not comparable. A raw-energy score would need an unbound reference. This one dodges that by ignoring absolute energy entirely.

---

<a id="c3"></a>
## C3 — HIGH — `beta = 0.01` means 12,000 K

**Where:** `run.py:53`, `maws2023.py:73`, `routines.entropy_score(sample, beta=0.01)`

`beta` is in mol/kJ. At 300 K the physical value is `1/RT = 1/(0.0083145 × 300) = 0.401`.

The default of 0.01 is **40× smaller**. That is `T = 1/(0.0083145 × 0.01) ≈ 12,000 K`.

At that temperature the distribution over sensible poses is nearly flat — 100 kJ/mol of difference gives a weight ratio of only `e¹`. What is **not** flattened is the clash population, which underflows to zero at any `beta`.

So `beta = 0.01` is exactly the regime where the score turns into "count the clashes":

| beta | T | good poses only | 50% clashes | discriminates on |
|---|---|---|---|---|
| 0.401 (300 K) | 300 K | −8.2572 | −8.5172 | nothing (at the floor) |
| 0.100 | 1,203 K | −7.2102 | −7.9891 | clash rate |
| **0.010 (default)** | **12,027 K** | **−0.1227** | **−0.8157** | clash rate |
| 0.001 | 120,272 K | −0.0012 | −0.6944 | clash rate |

At a physical `beta`, every candidate sits near the `−log 5000 = −8.52` floor, because real energy gaps are much bigger than `RT`. Nothing discriminates.

At the default `beta`, it does discriminate — on clash rate.

**There is no `beta` that makes this criterion rank binders.** That is why C1 and C2 are design problems, not tuning problems.

---

<a id="c4"></a>
## C4 — MEDIUM — the reported `energy` is not a binding energy, and is never used

**Where:** `run.py:248-252, 330-334`, `MawsResult.energy`

`free_E` is the lowest **total potential energy of the whole complex** over the sampled poses.

Three things wrong with reporting it as an energy:

1. **Not a binding energy.** No unbound reference. It is dominated by the target's own internal energy — measured at ~1.49 × 10⁴ kJ/mol for 1BRQ plus one nucleotide, essentially all intramolecular.
2. **Not comparable between candidates.** They have different atom counts.
3. **Never used to choose anything.** `run.py:265` and `:348` branch only on `entropy`.

So the number in `MawsResult.energy`, and the `ENERGY:` field in every log, is not an affinity and did not influence any decision the algorithm made.

**Fix — 2 hours.** Either compute a real interaction energy (complex − aptamer alone − target alone, same geometry — the machinery is already there) or stop calling it an energy.

---

<a id="c5"></a>
## C5 — MEDIUM — greedy growth, no backtracking

**Where:** `run.py:278-351`

The sequence grows one nucleotide at a time. Each step tries 8 candidates (4 bases × append/prepend), commits the best, and never revisits it.

A wrong choice at step 2 constrains every later step. No beam. No restart. No way to escape.

This is faithful to the original MAWS method, so it is a known limit rather than a regression. But combined with C1, one clash-driven mistake early propagates through the whole design.

**Fix — 3 hours.** Keep the top-k candidates per step. `k = 3` triples the cost and makes the search a lot more robust.

---

## C6 — Note — the scipy refactor is correct

For the record, since this change is currently in flight on `refactor/scipy-entropy-score`:

```
new entropy_score([100,150,200,175], 0.01) = -0.07243283074488105
old mpmath S       (same input)           = -0.07243283074488094
```

The `logsumexp` version is numerically better than the old mpmath one and drops a dependency. It reproduces the old behaviour exactly.

That is the problem, not the fix.
