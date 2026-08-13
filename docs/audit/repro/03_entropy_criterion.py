"""Reproduce findings C1, C2 and C3: what the MAWS entropy score actually rewards.

No AmberTools or OpenMM needed - this is pure numerics on synthetic energy sets.

    conda activate maws
    python docs/audit/repro/03_entropy_criterion.py
"""

import numpy as np

# main calls it S; refactor/scipy-entropy-score renamed it to entropy_score.
# float() because main's mpmath version returns an mpf, which cannot be formatted.
import maws.routines as _routines

_SCORE = getattr(_routines, "entropy_score", None) or _routines.S


def entropy_score(sample, beta=0.01):
    return float(_SCORE(sample, beta))


N = 5000
R = 0.0083145  # kJ/mol/K


def main():
    rng = np.random.default_rng(0)
    good = rng.normal(-1000.0, 50.0, N)

    print("=" * 78)
    print("C1: adding pure steric clashes to a fixed set of good poses")
    print("=" * 78)
    print("Candidate A: 5000 clash-free poses, energies ~ N(-1000, 50) kJ/mol.")
    print("Candidate B: the SAME good poses, but a fraction f replaced by hard")
    print("             clashes at +1e8 kJ/mol - i.e. a strictly WORSE candidate.")
    print("MAWS selects the MINIMUM score (run.py:265, run.py:348).\n")
    print("  clash fraction | entropy_score |  log(1-f)  | verdict")
    print("  ---------------+---------------+------------+-----------------")
    base = entropy_score(good, 0.01)
    print(f"  0.000 (A)      |    {base:+.4f}    |     --     | the good candidate")
    for f in (0.10, 0.50, 0.90, 0.99, 0.998):
        e = good.copy()
        e[rng.choice(N, int(f * N), replace=False)] = 1e8
        s = entropy_score(e, 0.01)
        print(
            f"  {f:5.3f}          |    {s:+.4f}    |  {np.log(1 - f):+.4f}   |"
            f" {'BEATS A' if s < base else 'loses to A'}"
        )
    print(f"\n  theoretical minimum = -log(N) = {-np.log(N):.4f}")
    print("  -> score tracks log(fraction not clashing): more clashes wins.")

    print()
    print("=" * 78)
    print("C2: the score is shift-invariant, so it is blind to binding strength")
    print("=" * 78)
    for shift in (-5000.0, 0.0, +5000.0):
        print(
            f"  energies shifted by {shift:+8.1f} kJ/mol -> "
            f"score {entropy_score(good + shift, 0.01):+.6f}"
        )
    print("  -> a strong binder and a non-binder are indistinguishable.")

    print()
    print("=" * 78)
    print("C3: there is no beta at which this criterion ranks binders")
    print("=" * 78)
    half = good.copy()
    half[rng.choice(N, N // 2, replace=False)] = 1e8
    print("   beta   |      T (K) | good poses only |  50% clashes | discriminates?")
    print("  --------+------------+-----------------+--------------+---------------")
    for b in (0.401, 0.100, 0.010, 0.001):
        g, h = entropy_score(good, b), entropy_score(half, b)
        note = "clash rate" if abs(g - h) > 0.05 and g > -8 else "nothing (at floor)"
        star = " <- default" if b == 0.010 else ""
        print(
            f"   {b:<6} | {1 / (R * b):10.0f} |     {g:+.4f}     |   "
            f"{h:+.4f}    | {note}{star}"
        )
    print("\n  At beta=0.401 (300 K) every candidate is pinned near the -log(N) floor.")
    print("  At beta=0.01 it discriminates - on clash rate.")

    print()
    print("=" * 78)
    print("C6: the scipy refactor is a faithful reimplementation of the mpmath S()")
    print("=" * 78)
    e = np.array([100.0, 150.0, 200.0, 175.0])
    w = np.exp(-0.01 * e)
    p = w / w.sum()
    old = -(p * np.log(p * len(p))).sum()
    new = entropy_score(e, 0.01)
    print(f"  old S formula   : {old!r}")
    print(f"  new entropy_score: {new!r}")
    print(f"  difference       : {abs(old - new):.3e}")


if __name__ == "__main__":
    main()
