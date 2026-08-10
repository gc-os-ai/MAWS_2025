# 07 — Migration Plan

**Start with step 1.** It is 2 hours and does not depend on anything else.

**Rule for every step:** one PR, `main` stays green, no step needs a later one to
have landed. No big rewrite.

**Total:** about 3 weeks of focused work.

---

## Do these first (about 5 days)

Additive and independent. Work them in any order, or in parallel.

| # | PR | Time | Risk | Changes results? |
| --- | --- | --- | --- | --- |
| 1 | Fix `rotate_in_residue` applying torsions 3 times | 2 h | low | **yes, numbers** |
| 2 | Add `AtomRange` / `Torsion` / `NucleotideSequence` value types | 1.5 d | low | no |
| 3 | Add `EnergyModel` protocol + `StubEnergy` | 1 d | low | no |
| 4 | Add `Builder` protocol + `FakeBuilder` | 1 d | low | no |

## Then these (about 2 weeks)

Each replaces working code. Ordered so the riskiest is last.

| # | PR | Time | Risk | Net lines |
| --- | --- | --- | --- | --- |
| 5 | Unify the two search loops | 3 d | medium | **−350** |
| 6 | Add `maws.cli` + `[project.scripts]` + `AptamerDesigner` | 3 d | medium | +260 |
| 7 | `Pose` as an `(N,3)` ndarray | 2 d | medium | **−200** |
| 8 | Split `ResidueLibrary` from `LeapResources` | 1.5 d | medium | −80 |
| 9 | `ChainView` replaces `Chain`; cycle removed | 4 d | **high** | **−400** |
| 10 | Remove the `Complex` shim | 0.5 d | low | −150 |

Step 5 pays for the whole effort. Step 9 needs a quiet week.

---

## Step 1: Fix `rotate_in_residue` first, and alone

**Time: 2 hours. Do this before anything else.**

If this is a real bug ([01-values.md](01-values.md#2-the-bug-this-makes-impossible)),
fixing it changes the numbers every run produces. Do it inside a larger refactor
and you cannot tell which change caused what.

### 1a. Confirm it (30 min)

The claim is a reading of the code, not a run. Check it first:

```python
def test_rotate_in_residue_applies_torsion_exactly_once():
    """One call should rotate by `angle`, not by 3 x angle."""
    cpx, chain = _minimal_built_complex()
    before = np.array([nostrom(p) for p in cpx.positions])

    chain.rotate_in_residue(0, 0, np.pi / 6)
    once = np.array([nostrom(p) for p in cpx.positions])

    cpx.positions = _restore(before)
    for _ in range(3):
        chain.rotate_in_residue(0, 0, np.pi / 6)
    thrice = np.array([nostrom(p) for p in cpx.positions])

    # If the loop bug is real, `once` already equals `thrice`.
    assert not np.allclose(once, thrice), "one call is behaving like three"


def test_rotate_in_residue_does_not_mutate_the_template():
    """The Structure template must be unchanged after a rotation."""
    cpx, chain = _minimal_built_complex()
    snapshot = copy.deepcopy(chain.structure.rotating_elements)
    chain.rotate_in_residue(0, 2, 0.1)          # a torsion with negative indices
    assert chain.structure.rotating_elements == snapshot
```

### 1b. Fix it (30 min)

Move the call out of the loop. Stop pointing at the template.

```python
element = list(self.structure.rotating_elements[resname][residue_element_index])  # copy

n = self.structure.residue_length[resname]
element = [i + n if (i is not None and i < 0) else i for i in element]   # `is not None`

offset = self.residues_start[revised_residue_index]
if element[2] is None:
    revised = [element[0] + offset, element[1] + offset, None]
elif element[2] == 0:
    revised = [element[0] + offset, element[1] + offset, 0]
else:
    revised = [e + offset for e in element]
    rev = False

self.rotate_element(revised, angle, reverse=rev)      # once, outside the loop
```

`is not None` matters. Today's `if element[i] and element[i] < 0` treats index `0`
as false, so a real atom index of 0 skips the check.

### 1c. Record the new baseline (1 h)

Run `--length 5` on `data/pfoa.pdb`. Put the resulting sequence and energy in the
PR body. That is the number every later step compares against.

**Done when:** both tests pass and the baseline is recorded.

---

## Step 2: Value types (1.5 days)

Add `maws/values.py`. Nothing imports it except its own tests.

Give `ResidueLibrary.from_tables(...)` the **existing** parallel-array format, so
[`rna_structure.py`](../../maws/rna_structure.py) and
[`dna_structure.py`](../../maws/dna_structure.py) need no edits:

```python
def test_rna_library_matches_legacy_structure():
    """The new library must agree with Structure on every residue."""
    legacy = load_rna_structure()
    lib = ResidueLibrary.rna()

    assert set(lib) == set(legacy.residue_names)
    for name in legacy.residue_names:
        assert lib[name].n_atoms == legacy.residue_length[name]
        assert len(lib[name].torsions) == len(legacy.torsions(name))
```

This test is the safety net for every later step. It pins the new types to the old
behaviour before anything depends on them.

---

## Step 3: `EnergyModel` protocol (1 day)

Add `maws/energy.py`. Make `Complex` satisfy the protocol by delegating, so
nothing breaks:

```python
class Complex:
    def evaluate(self, pose) -> float:            # new, satisfies EnergyModel
        return self.get_energy()[0]

    def get_energy(self):                          # unchanged, still works
        ...
```

Add `StubEnergy`. It has no production caller yet. It exists so step 5 can be
tested.

**Also delete `rigid_minimize` in this PR** (question 2). No callers, no tests,
and it raises `TypeError` on its first pass. Carrying it into a clean layer moves
a fault along with it.

---

## Step 4: `Builder` protocol (1 day)

Pull `Complex.build` and `_build_cache_key` into `LeapBuilder`. `Complex.build`
becomes 2 lines of delegation. Add `FakeBuilder`, which returns a canned
`BuiltSystem` from a small checked-in `.prmtop`/`.inpcrd` fixture.

Do the `LIG` naming fix here too: `make_lib`'s `residue_name` comes from a digest
of the source PDB.

**This invalidates the cache.** Existing `.maws_cache/` entries will miss and
rebuild once. Say so in the PR.

---

## Step 5: Unify the search (3 days)

Add `maws/search.py` with `grow_aptamer`. Rewrite both callers to consume it.

**This must produce identical results.** Two tests:

```python
@pytest.mark.integration
def test_unified_search_matches_legacy_runner(tmp_path):
    """grow_aptamer must reproduce MawsRunner.run exactly, same seed."""
    legacy = MawsRunner(num_nucleotides=3, aptamer_type="RNA",
                        molecule_type="organic", first_chunk_size=50,
                        second_chunk_size=50).run(pdb=FIXTURE, name="legacy")

    new = design(FIXTURE, length=3, aptamer="RNA", molecule="organic",
                 samples=50, first_samples=50, seed=_same_seed)

    assert new.sequence == legacy.sequence
    assert new.energy == pytest.approx(legacy.energy, rel=1e-9)


def test_cli_and_library_produce_identical_events():
    """The CLI reporter and design() must see the same event stream."""
```

Seeding has to come first, so fold the `rng` threading from
[05-search.md](05-search.md#4-two-changes-to-sampling) into this PR. Without it
the equivalence test cannot be written.

**Done when:** `run.py` and `maws2023.py` are each about 80 lines.

---

## Step 6: CLI and estimator (3 days)

Two deliverables, same PR, because both are the public API.

**6a. `AptamerDesigner` (1 day).** Add `maws/_designer.py` following the
scikit-learn estimator contract — see
[06-cli.md](06-cli.md#what-to-build-a-scikit-learn-estimator). Add `scikit-learn`
to **both** `environment.yml` and `pyproject.toml`.

Gate it on the conformance suite:

```python
from sklearn.utils.estimator_checks import parametrize_with_checks

@parametrize_with_checks([AptamerDesigner()])
def test_sklearn_compliance(estimator, check):
    check(estimator)
```

Skip inapplicable checks by name with a reason — `X` is a list of file paths, not
a numeric array, so array-shape checks will not apply. What the suite does enforce
for free: `clone()` correctness, `get_params` round-tripping, and "no work in
`__init__`".

**6b. CLI (2 days).** Add `maws/cli.py` and the `[project.scripts]` block. The CLI
builds an `AptamerDesigner` from parsed arguments and calls `fit`, so there is one
configuration path, not two. Keep `maws2023.py` as a deprecating shim:

```python
# maws/maws2023.py
"""Deprecated. Use `maws design`. Kept through v0.3."""

def main() -> None:
    warnings.warn(
        "`python -m maws.maws2023` is deprecated; use `maws design`. "
        "This module will be removed in v0.3.",
        DeprecationWarning, stacklevel=2,
    )
    from maws.cli import main as _main
    _main(["design", *_translate_legacy_argv(sys.argv[1:])])
```

Add `maws/py.typed` and the `package-data` entry in the same PR. The package is
fully annotated but ships no marker, so none of those hints reach users today. See
[08-patterns.md](08-patterns.md#5a-ship-pytyped).

`_translate_legacy_argv` maps `-nt` to `--length`, `-p` to `--target`, and so on.
Existing scripts and
[`tests/test_maws2023_cli.py`](../../tests/test_maws2023_cli.py) keep passing.

---

## Step 7: `Pose` (2 days)

Swap `list[Vec3]` for an `(N,3)` float64 array. Convert at the OpenMM edge only.

Two things must both pass:

```python
def test_pose_rotation_matches_legacy_within_tolerance():
    """The ndarray path must match the Vec3 path to float64 precision."""
    legacy = _legacy_rotate(positions, element, angle)
    new = Pose(xyz, system).rotate(torsion, angle)
    assert np.allclose(new.xyz, legacy, atol=1e-10)
```

```console
$ pytest tests/bench_pose.py --benchmark
```

Put the before/after numbers in the PR. [03-pose.md](03-pose.md) says the speedup
should be large but does not give a figure, because nothing has been measured yet.
This PR is where the real number gets written down.

---

## Step 8: Split `Structure` (1.5 days)

`ResidueLibrary` (pure) and `LeapResources` (files) become the real path.
`Structure` becomes a shim over both. Low risk, because step 2's equivalence test
already pins the behaviour.

---

## Step 9: `ChainView`, and the cycle goes (4 days)

The risky one, and the reason it is last. It touches the 180 lines of junction
geometry in `rebuild` ([complex.py:511-689](../../maws/complex.py#L511-L689)),
which is the least-tested code in the package.

Three sub-steps, not one:

**9a (1.5 d).** Pull the 3 `rebuild` branches into named free functions:
`_splice_head`, `_splice_tail`, `_splice_whole`. No other change. Add unit tests
against 20-atom poses built by hand. This is where those two coordinate-mapping
paths finally get direct coverage.

**9b (1 d).** Add the `Edit` value. `rebuild` takes it as a parameter instead of
reading `*_history`. The fields stay, written but unread, so nothing breaks.

**9c (1.5 d).** Add `ChainView`. Delete `Chain`, the history fields, and
`update_chains`. `Complex` becomes a thin wrapper over `BuiltSystem` and `Pose`.

Delete the 2 dead "historic" methods in 9c. Already confirmed by grep: only their
own definitions, no callers.

```console
$ grep -rn "rotate_historic\|rotate_in_historic" --include="*.py" .
```

---

## Step 10: Remove the shim (0.5 days)

Only once nothing inside the repo uses `Complex`. Waits on question 3: does
anything outside the repo import it?

---

## Tests

The suite has 13 files. Most cannot run without AmberTools. Invert that.

| Tier | Needs | Speed | Covers after migration |
| --- | --- | --- | --- |
| unit | nothing | ms | values, sequences, offsets, pose math, entropy, config merge |
| algorithm | `StubEnergy` + `FakeBuilder` | seconds | growth loop, selection, events, CLI wiring |
| integration | AmberTools + OpenMM | minutes | LEaP build, cache, real energies, end to end |

Keep the `integration` marker in [`pyproject.toml`](../../pyproject.toml). The
middle tier is the point. Nothing between "pure function" and "full toolchain"
exists today, so the algorithm has no coverage at all.

Per your `feedback_dual_dep_paths` note, a new dependency has to go in **both**
`environment.yml` and `pyproject.toml`. This plan adds none. `tomllib` is stdlib
on 3.11+, which is the reason to bump `requires-python` in step 6.

---

## Rollback

Steps 1 to 4 and 6 are additive. Revert the commit.

Steps 5, 7, and 9 replace working code, so each keeps the old path reachable for
one release:

- Step 5: `MawsRunner.run` keeps its signature. Only the body changes.
- Step 7: `Pose.from_openmm` and `.to_openmm` bracket the change. The old rotation
  stays in the test module as `_legacy_rotate`, the thing new results are checked
  against.
- Step 9: sub-steps 9a and 9b revert on their own. Only 9c deletes anything.

---

## Done looks like

```console
$ pip install -e .
$ maws design --target data/pfoa.pdb --length 20 --aptamer DNA --seed 42 -o runs/
Step 1/20: seeding
  G 3prime -> entropy=-0.072433 energy=-1421.88
  ...
Step 20/20 complete. Best: G A T C G A T C ...
Done after 20 steps (E=-1523.71, S=-0.084216)
Wrote runs/MAWS_aptamer_RESULT.pdb

$ pytest -m "not integration"          # the algorithm, no AmberTools
147 passed in 3.2s
```

**Next:** start step 1. Open [chain.py:341](../../maws/chain.py#L341).
</content>
