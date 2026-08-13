# Routines API Reference

`maws.routines` provides the scoring function that drives nucleotide selection in the MAWS aptamer-design loop.

## Overview

Each step of the MAWS search samples many conformations of a candidate aptamer against the ligand and records their potential energies. `maws.routines` turns that energy sample into a single number. The one public entry point is [`entropy_score`](#entropy_score); MAWS keeps the candidate whose score is **lowest**.

### How `routines` connects to the rest of MAWS

`entropy_score` is the single seam between this module and everything else. Both the CLI script (`maws.maws2023.main`) and the programmatic API class (`maws.run.MawsRunner.run`) call it once per candidate nucleotide, per search step, on the list of energies gathered in that step's sampling loop.

```mermaid
flowchart LR
    SAMP[sampling loop] -->|energies: list of float| ES[entropy_score]
    BETA[beta] --> ES
    ES -->|score| SEL[keep lowest-scoring candidate]
```

---

## `entropy_score`

```python
entropy_score(sample, beta=0.01) -> float
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| `sample` | array-like of `float` | *required* | Energy values (kJ/mol) from conformational sampling |
| `beta` | `float` | `0.01` | Inverse temperature (1/kT in mol/kJ). Higher `beta` = sharper distribution |

Returns a `float`. Raises `ValueError` if `sample` is empty — there is no distribution to score, so no return value would be meaningful.

The energies are converted to a Boltzmann distribution `P(i) = exp(-beta * E_i) / Z`, and the score is

```
-sum(P * log(P * N))
```

which is the negative Kullback–Leibler divergence of that distribution from the uniform distribution over the `N` samples.

### Interpreting the score

The score is at most `0`, reaching `0` exactly when every sampled energy is equal. It falls toward `-log(N)` as the distribution concentrates on a few conformations. So a **more negative score means a more peaked distribution**, which MAWS treats as evidence of stronger, more specific binding — hence selection by minimum.

Two consequences worth knowing:

- **It is shift-invariant.** Adding a constant to every energy leaves the score unchanged; only the *spread* of the energies matters, not their absolute level. The absolute binding energy is tracked separately, as `MawsResult.energy`.
- **It depends on `N`.** Because the reference is the uniform distribution over exactly the sampled points, scores are only comparable between candidates evaluated with the same sample count. MAWS satisfies this by drawing a fixed number of samples per search step (the chunk size).

### Choosing `beta`

`beta` sets how sharply energy differences are weighted. As `beta → 0` all conformations weigh equally and the score flattens toward `0` for every candidate, losing discrimination; as `beta` grows the score is dominated by the single lowest-energy conformation. The default of `0.01` mol/kJ is inherited from the original MAWS implementation and is exposed as `MawsRunner(beta=...)` and `--beta` on the CLI.

### Usage

```python
from maws.routines import entropy_score

energies = [-1500.0, -1490.0, -1450.0, -1200.0]
score = entropy_score(energies, beta=0.01)
```

Inside the search loop, the pattern is:

```python
best_score = None
for ntide in nucleotides:
    energies = [sample_one_pose() for _ in range(chunk_size)]
    score = entropy_score(energies, beta=self.beta)
    if best_score is None or score < best_score:
        best_score = score
        best_sequence = ntide
```

---

## Numerical notes

Weights are computed in log-space with `scipy.special.logsumexp`, and the score itself with `scipy.stats.entropy`. This matters because OpenMM energies for clashing poses routinely reach 1e5–1e6 kJ/mol: at `beta=0.01` the corresponding weights are around `exp(-1e4)`, which flushes to zero in double precision if exponentiated directly, while strongly negative energies overflow `Z` to infinity in the other direction. Shifting by the maximum log-weight removes both failure modes.

Conformations whose weights do underflow after the shift carry probability on the order of `exp(-1e4)` and contribute nothing measurable to the score, so no accuracy is lost by dropping them.

> **Note (2026): implementation change.** This module previously used `mpmath` at 60 decimal digits of precision to avoid the same underflow. Arbitrary precision was never the right tool — the log-sum-exp shift addresses the problem directly in double precision, and is ~250× faster. Scores from the two implementations agree to within 2e-15 absolute across energy ranges from 1e0 to 1e6 kJ/mol and `beta` from 0.001 to 0.1, so results are unchanged. `mpmath` is no longer a dependency.

## Internal helpers

`_boltzmann(sample, beta)` returns `(P, log_z)`: the normalised probability array and the natural log of the partition function. It is private and may change without notice. Note that it returns **log** `Z` rather than `Z`, because `Z` itself overflows double precision for strongly negative energies — the log is always finite.
