# 08 — Patterns and Library Conventions

**Read this if you want to know what the redesign is called** and which OSS Python
libraries it copies.

**The closest model is OpenMM itself** — the library MAWS already depends on. The
redesign is largely "do what OpenMM does, one level up."

---

## 1. OpenMM's own split, applied to MAWS

OpenMM separates description, parameters, state, and driver. MAWS currently
collapses all four into `Complex`.

| OpenMM | What it holds | MAWS equivalent (proposed) |
| --- | --- | --- |
| `Topology` | what atoms and bonds exist | `Assembly` — what chains exist |
| `System` | forces and parameters | `ForceField` + `BuiltSystem` |
| `State` | a snapshot of positions | `Pose` |
| `Context` | live mutable simulation state | `OpenMMEnergy` (the only owner) |
| `Simulation` | convenience wrapper over the above | `AptamerDesigner` |

This is the strongest argument for the layering. A MAWS user who already knows
OpenMM will recognise the shape. Right now they have to learn a `Complex` that is
`Topology`, `System`, `Context`, and `Simulation` at once.

Two OpenMM habits worth copying directly:

- **`State` is immutable and you ask for what you want.**
  `getState(getPositions=True, getEnergy=True)`. `Pose` follows this.
- **`Context` is the only mutable thing, and it is explicit.** You never wonder
  where the state lives. [04-energy.md](04-energy.md) says `OpenMMEnergy` is the
  only class holding a `Context` for the same reason.

---

## 2. scikit-learn conventions

The public API follows the scikit-learn estimator contract. Full design in
[06-cli.md](06-cli.md#what-to-build-a-scikit-learn-estimator). The eight rules
that matter, and what each fixes here:

| Rule | What it fixes in MAWS |
| --- | --- |
| `__init__` takes hyperparameters only, stores them verbatim | `Structure.__init__` reads files; `add_chain_from_pdb` shells out |
| No validation or computation in `__init__` | `MawsRunner.__init__` validates ad hoc ([run.py:65-74](../../maws/run.py#L65-L74)) |
| `fit` does the work and returns `self` | no single verb today; `run()` returns a result |
| Fitted attributes end with `_` | nothing marks what exists only after a run |
| Hyperparameters never end with `_` | — |
| `random_state`, resolved by `check_random_state` inside `fit` | no run is reproducible today |
| `get_params`/`set_params` free from `BaseEstimator` | no way to introspect or clone a configured run |
| Private module, public re-export | `maws/_designer.py` -> `from maws import AptamerDesigner` |

**"No work in `__init__`" is the rule this codebase breaks hardest.**
`Structure.__init__` reads the filesystem
([structure.py:72-87](../../maws/structure.py#L72-L87)) and
`Complex.add_chain_from_pdb` runs `antechamber`, `parmchk2`, and `tleap`
([complex.py:154](../../maws/complex.py#L154)). Constructors that shell out are
why nothing can be unit-tested.

### Which estimator shape

MAWS is unsupervised and transductive: every target needs a fresh search, and
there is no fitted model to apply to new data. That is the `TSNE` /
`SpectralClustering` shape — `fit` and `fit_predict`, no `predict`, declared via
`__sklearn_tags__`.

### Function and class, both

scikit-learn ships both forms of the same algorithm: `sklearn.cluster.k_means()`
next to `KMeans`, `ridge_regression()` next to `Ridge`. The class is a thin
wrapper over the function.

`grow_aptamer()` is the function; `AptamerDesigner` wraps it. That is why
[05-search.md](05-search.md) can keep the event-stream generator without it
becoming a second public API to maintain.

### Result objects

`scipy.optimize.minimize` returns an `OptimizeResult` with `.x`, `.fun`,
`.success`, `.message`. `sklearn.utils.Bunch` gives dict access plus attribute
access. Both beat a bare tuple. `Relaxed`, `Candidate`, and `MawsResult` follow
this — including `.success` and `.message`, which `MawsResult` lacks today.

### Other scientific-Python conventions

| Convention | Who does it | Where used |
| --- | --- | --- |
| Domain objects, not raw tuples | MDAnalysis `AtomGroup`, RDKit `Mol` | `AtomRange`, `Torsion`, `ChainView` |
| Views over a shared array | MDAnalysis, pandas, numpy | `ChainView`, `ResidueView` |
| NumPy-style docstrings | numpy, scipy, sklearn, pandas | already used throughout |

---

## 3. Why `AtomRange` and not `range`

An earlier draft of this doc argued for using the built-in `range`, since it is
already a frozen half-open integer span. That argument was about matching the
standard library. It is the wrong target.

Scientific Python libraries consistently define domain types rather than reuse
built-ins:

- MDAnalysis has `AtomGroup`, not `list[int]`
- RDKit has `Mol` and `Conformer`, not arrays
- Biopython has `SeqRecord`, not `str`
- pandas has `Index`, which is a labelled array, not `range`

The reason is that a domain type carries validation, a docstring, a repr, and room
to grow. `def place(span: AtomRange)` documents itself; `def place(span: range)`
does not.

**Keep `AtomRange`.** It validates inverted ranges, and `as_slice()` gives you the
numpy-view behaviour explicitly rather than by accident.

---

## 4. Architectural patterns

Two organizing patterns. Everything else is local.

| Pattern | What it means here |
| --- | --- |
| **Ports and Adapters** (Hexagonal) | Pure core, side effects behind protocols |
| **Functional core, imperative shell** | Layers 0-3 pure; Layer 4 and adapters do I/O |

| Part | Role |
| --- | --- |
| `values.py`, `topology.py` | domain core — no I/O, no OpenMM |
| `Builder`, `EnergyModel`, `Sampler`, `Scorer` | ports (protocols) |
| `LeapBuilder`, `OpenMMEnergy`, `SurfaceSampler` | adapters (real) |
| `FakeBuilder`, `StubEnergy`, `FixedSampler` | adapters (test) |
| `cli.py`, `_designer.py` | driving side |

The test adapters are not a bolt-on. They are the second implementation the
pattern requires you to have. **If only one implementation of a port ever exists,
the port was not worth adding.**

### Local patterns

| Component | Pattern | Library doing the same |
| --- | --- | --- |
| `AptamerDesigner` | Estimator | sklearn `TSNE`, `SpectralClustering` |
| `Torsion`, `ForceField` | Value Object | RDKit `Conformer`, frozen dataclasses |
| `ResidueLibrary` | Immutable Mapping | Biopython `CodonTable` |
| `Assembly` + `build()` | Immutable Builder | `pathlib.Path` joins |
| `EnergyModel`, `Scorer` | Strategy | sklearn `scoring=`, `logging.Formatter` |
| `CompositeEnergy` | Composite | sklearn `Pipeline`, `FeatureUnion` |
| `grow_aptamer()` | Generator | `sklearn.cluster.k_means()` functional form |
| `StepEvent` subclasses | Tagged union + `match` | `ast` nodes |
| `ChainView` | View | MDAnalysis `AtomGroup`, `dict.keys()` |

---

## 5. Library conventions the redesign is missing

These are things well-run OSS Python packages do that docs 01 through 07 do not
mention. Ranked by how much they matter here.

### 5a. Ship `py.typed`

The package is fully type-annotated but has no `py.typed` marker, so **none of
those hints reach users.** `mypy` and IDEs silently treat `maws` as untyped.

```toml
# pyproject.toml
[tool.setuptools.package-data]
maws = ["py.typed"]
```

Plus an empty `maws/py.typed`. 5 minutes. Do it regardless of the redesign.

### 5b. `__repr__` that round-trips

Convention in `attrs`, `pydantic`, and dataclasses: `repr` shows how to rebuild
the object. Frozen dataclasses give this free. `Pose` and `BuiltSystem` are not
dataclasses, so write theirs by hand:

```python
def __repr__(self) -> str:
    return f"<Pose {len(self._xyz)} atoms, {len(self._system.assembly.chains)} chains>"
```

Today `print(cpx)` gives `<maws.complex.Complex object at 0x7f...>`.

### 5c. A stated deprecation policy

[06-cli.md](06-cli.md) says `maws2023` is "kept through v0.3" but nothing says
what a version number means. Adopt the common scientific-Python policy: deprecate
for two minor releases, `DeprecationWarning` with the replacement named and the
removal version.

```python
warnings.warn(
    "Complex is deprecated since 0.2 and will be removed in 0.4. "
    "Use maws.build(Assembly(...)) instead.",
    DeprecationWarning, stacklevel=2,
)
```

### 5d. Entry points for plugins

If third parties should add scorers or samplers (a folding-ΔG model, a docking
score), the convention is `entry_points`, as used by pytest, sphinx, and flake8:

```toml
[project.entry-points."maws.scorers"]
entropy = "maws.routines:entropy_score"
folding = "maws.contrib.folding:folding_score"
```

Only worth doing if you expect outside contributions. Skip otherwise.

### 5e. `NumPy`-style docstrings — already correct

Worth stating because it is easy to lose. This repo already uses NumPy-style
docstrings throughout, matching numpy, scipy, sklearn, and pandas. Keep them. Do
not switch to Google style mid-redesign.

---

## 6. Patterns deliberately not used

| Pattern | Why not |
| --- | --- |
| **Visitor** for `StepEvent` | `match` on dataclasses covers it in 3.10+. Visitor adds a class per consumer. |
| **Observer** for events | Push callbacks invert control. A generator lets the caller `break`, which is what makes early stopping work in [06-cli.md](06-cli.md). |
| **Abstract Factory** for force fields | One classmethod and two dicts cover all 6 combinations. |
| **Singleton** for `ResidueLibrary` | `functools.cache` on the classmethod does it without global state. |

The risk is real: a redesign trying to look principled collects patterns it does
not need. **If a pattern here is not earning either a deleted concept or a second
implementation, drop it.**

---

## 7. What is not a pattern

Being honest about which parts are just decomposition:

- Splitting `Structure` into `ResidueLibrary` and `LeapResources` is separation of
  concerns.
- Merging the two search loops is removing duplication.
- Replacing `list[Vec3]` with an ndarray is a data representation choice.
- The exception hierarchy is convention. (It does follow the common idiom, though:
  `MissingLibrary(FileNotFoundError)` mirrors `json.JSONDecodeError(ValueError)`,
  so old `except FileNotFoundError` still catches it.)

Those four account for most of the line-count reduction. The patterns organize the
result. They are not what does the work.

---

## 8. Changes to apply to the other docs

| Doc | Change | Time |
| --- | --- | --- |
| 01, 02, 03, 05, 06 | `Sequence` -> `NucleotideSequence` (collides with `collections.abc`) | done |
| 04-energy | `Relaxed` becomes a `NamedTuple`, so `energy, pose = ...` still unpacks | 5 min |
| 05-search | `Sampler` implements `__next__`, so `islice`/`zip`/`for` work free | 10 min |
| 06-cli | `MawsResult` gains `.success` and `.message`, like `OptimizeResult` | 10 min |
| 06-cli | `AptamerDesigner(BaseEstimator)` replaces `design()`/`DesignRun` | done |
| 07-migration | Add `py.typed` to step 6 | 5 min |

**Next:** add `scikit-learn` to both `environment.yml` and `pyproject.toml`, then
write `AptamerDesigner.__init__` with nothing in the body but attribute
assignments. That one constructor is the whole convention in miniature.
</content>
