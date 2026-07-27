# Class G — Missing tests

The bugs in Classes A, B and F are not subtle once you look. The real question is why eleven test files caught none of them.

---

<a id="g1"></a>
## G1 — HIGH — the rotation code is mocked out of its own tests

**Where:** [tests/test_chain.py:19-34](../../tests/test_chain.py#L19-L34)

```python
class MockComplex:
    """Minimal mock of Complex for Chain unit tests.

    Chain only uses Complex for:
    - storing reference (self.complex = Complex)
    - iterating chains (Complex.chains) in update_chains()
    - calling rotate_element() which we don't test here
    """
    def __init__(self):
        self.chains: list[Chain] = []

    def rotate_element(self, element, angle, reverse=False):
        """Mock rotation - does nothing in unit tests."""
        pass
```

Every `Chain` rotation test runs against a no-op.

`Chain.rotate_in_residue` applies each torsion three times using broken indices ([A1](01-conformational-search.md#a1)). Its tests pass, because the function it calls does nothing.

`grep -rn "rotate_in_residue\|rotate_element\|rotate_global\|reverse=True" tests/` returns two hits. Both are the mock above and its docstring.

**No test anywhere calls the real rotation code. No test anywhere checks a bond length, a bond angle, or any geometric invariant.**

That one gap explains A1, A2, A4, A7, and why F2 cannot be checked.

---

## What to write, ranked by payoff

Top 5. Do them in this order.

| # | Test | Catches | Time |
|---|---|---|---|
| 1 | **Bond-length invariant.** After any rotation or rebuild, assert every covalent bond in the aptamer is within 0.05 Å of its pre-move length. Use the prmtop bond list. | A1, A2, A4, A7, F1, F2 | 30 min |
| 2 | **Applied-angle check.** Request θ on a known torsion, measure the resulting dihedral, assert it equals θ. | A1 (the 3× factor) | 20 min |
| 3 | **Placement check.** After `translate_global` + `rotate_global`, assert the aptamer's centroid is where the sampler said. | B1, B3 | 20 min |
| 4 | **Criterion property test.** Assert `entropy_score` does not improve when you add clash samples to a fixed good set. | C1 | 15 min |
| 5 | **Cleaner atom-count check.** Assert the cleaned PDB keeps a plausible fraction of the input atoms and the same number of `TER` records. | E2, E3, E4, E5 | 20 min |

Number 1 alone catches six findings across three classes.

Two more once those land:

- **LEaP output scan.** Assert `Errors = 0` and no `FATAL` in captured tleap stdout. Catches D2.
- **Seeded end-to-end determinism.** Two runs, same seed, identical output. Catches B5 and gives you regression cover for everything above.

Starting code for all of these is in [`repro/`](repro/). `repro/01_rotation_geometry.py` is close to being a real test already.

---

## About the docstrings

You said the docs cannot be trusted. That held up. But the failure mode is unusual and worth naming.

The docstrings are not stale or vague. They are **detailed, confident, and describe behaviour the code does not have.**

Five examples found during this audit:

1. **`complex.py:754-756`** documents `reverse=True` as rotating "the complement of the selected range relative to the pivot". The parameter never reaches that function — `chain.py:291` sets `rev = False` unconditionally first. ([A2](01-conformational-search.md#a2))
2. **`space.py:196-213`** gives a correct, well-cited account of the solvent-accessible surface, then applies it to the wrong object. ([B2](02-pose-sampling.md#b2))
3. **`routines.py:63-67`** says lower scores mean "a more peaked distribution, suggesting stronger binding affinity". First clause true. Second does not follow, and measurement contradicts it. ([C1](03-selection-criterion.md#c1))
4. **`pdb_cleaner.py:574`** documents `drop_hetatm=False` as preserving HETATM records. It does not. ([E7](05-input-preparation.md#e7))
5. **`complex.py:1025-1029`** (on the fix branch) accurately describes the `pert_min` rigidity problem. Included because it shows the difference — this one was written from behaviour.

The pattern: docstrings written from **intent** rather than from **behaviour**.

That is more dangerous than no docstring. It survives review. A reader checking `chain.py` against `complex.py`'s documented contract finds them consistent — and the contract is the thing that is wrong.

**The fix:** wherever a docstring makes a checkable claim, turn that claim into an assertion.
