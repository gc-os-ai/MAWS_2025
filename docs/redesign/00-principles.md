# 00 — Six Rules

Six rules. Every other doc applies one of them. If a design question comes up that
the other docs do not answer, answer it from here.

---

## P1 — Make bad states impossible to write

If a caller has to ask "which kind of value is this?", the producer should have
used a different type.

### Now

`Structure.rotating_elements[res]` can be four different things:

- `[None]` — no torsions defined
- `[start, bond, end]` — a normal torsion
- `end` may be `None` — meaning "to the end of the chain"
- any index may be negative — meaning "counting from the end of the residue"

So four places all check for the same conditions:

```python
# structure.py:393
if triples == [None]:                      # sentinel for "no rotations defined"
    return []

# complex.py:939 — same check, different file, different wording
rots = chain.structure.rotating_elements[residue]
if rots == [None]:                         # skip residues with no torsions
    continue

# chain.py:293
if len(revised_element) == 3 and revised_element[2] is not None:

# chain.py:347
if element[i] and element[i] < 0:
```

That last one has a bug hidden inside it. `if element[i] and ...` treats index `0`
as false, so a real atom index of 0 skips the check.

### Proposed

`ResidueTemplate.torsions` is a `tuple`. Empty means no torsions. No sentinel.
`Torsion.moving` is an `AtomRange` with a real `stop`. No `None`. Indices are made
absolute when the template loads. No negatives.

All four checks above get deleted.

---

## P2 — One coordinate frame per type, and say which

### Now

The same `element` list holds residue-local, chain-local, or global indices. You
cannot tell which. Conversions are written out by hand at every boundary:

```python
# chain.py:294 — chain-local to global
revised_element = [idx + self.start for idx in revised_element]

# chain.py:354 — residue-local to chain-local
element[0] + self.residues_start[revised_residue_index]

# chain.py:391 — old-global to current-global
historic_element[0] + self.start_history - self.start
```

### Proposed

`Torsion` always holds global indices. The only conversion is
`Torsion.shifted(offset)`, called in one place. Three hand-written idioms become
one method.

---

## P3 — Values are frozen; state has one owner

### Now

`Chain.rotate_in_residue` grabs a live reference into the shared `Structure` and
edits it ([chain.py:341-348](../../maws/chain.py#L341-L348)):

```python
element = self.structure.rotating_elements[...][residue_element_index]  # shared ref
for i in range(len(element)):
    if element[i] and element[i] < 0:
        element[i] += self.structure.residue_length[...]                # edits in place
```

Every `Chain` using that `Structure` sees the edit. Call it twice and you get two
different behaviours: the second call finds the negatives already gone.

### Proposed

Everything in Layer 0 is `@dataclass(frozen=True)` with `tuple` fields, not `list`.
`Pose.rotate()` returns a new `Pose`; the old one is unchanged.

Mutable state exists in exactly two places, both documented: the `LeapBuilder`
cache and the `OpenMMEnergy` context.

---

## P4 — No cycles; children do not edit parents

### Now

`Chain` holds `self.complex`. `Complex` holds `self.chains`. So a chain edits its
siblings, then undoes the edit for itself
([chain.py:149-161](../../maws/chain.py#L149-L161)):

```python
for chain in self.complex.chains:
    chain.start_history = chain.start
    if chain.start >= start:
        chain.start += self.length - old_length      # edit siblings
        chain.element = [chain.start, chain.start + 1, chain.start + chain.length]

self.start -= self.length - old_length               # undo it for self
self.element = [self.start, self.start + 1, self.start + self.length]
```

That last subtraction is not doing any work. It cancels a side effect the loop
above just caused.

### Proposed

`ChainView` holds no reference back to anything. `BuiltSystem` computes offsets in
one pass from the chain lengths. Nothing propagates, so nothing needs cancelling.

**Offsets should be computed, not passed around.**

---

## P5 — Side effects behind protocols

### Now

Nothing in the search path runs without AmberTools on `PATH` and a working LEaP.
`Structure.__init__` reads files and can raise `FileNotFoundError`
([structure.py:81](../../maws/structure.py#L81)). You cannot even build a residue
table in a test without `.lib` files on disk.

### Proposed

Three protocols wrap every side effect:

| Protocol | Real | For tests |
| --- | --- | --- |
| `Builder` | `LeapBuilder` (runs `tleap`) | `FakeBuilder` (canned result) |
| `EnergyModel` | `OpenMMEnergy` | `StubEnergy` (analytic, deterministic) |
| `Sampler` | `SurfaceSampler` | `FixedSampler` (replays a list) |

The whole growth algorithm becomes testable with no external tools. See
[04-energy.md](04-energy.md).

---

## P6 — One algorithm, many reporters

### Now

The growth search exists twice:

- [run.py:224-350](../../maws/run.py#L224-L350) — 130 lines, returns a `MawsResult`
- [maws2023.py:271-445](../../maws/maws2023.py#L271-L445) — 175 lines, writes logs and PDBs

Same algorithm. The only difference is what happens to each candidate after
scoring.

Both also re-derive the same force field mapping
([run.py:164-172](../../maws/run.py#L164-L172),
[maws2023.py:213-223](../../maws/maws2023.py#L213-L223)) and both build the same
throwaway ligand-only `Complex`.

Any change — a folding-ΔG filter, a different growth strategy, early stopping —
has to be made twice or the CLI and the library drift apart.

### Proposed

One generator yields typed events. The CLI consumes them to write logs. The
library consumes them to keep the best result. See [05-search.md](05-search.md).

---

## Worked example: one function, all six rules

```python
# chain.py:283-304, shortened
def rotate_element(self, element, angle, reverse=False):
    revised_element = element[:]
    rev = reverse
    if rev:
        if revised_element[2] is None:
            revised_element[2] = 0                    # 0 means "complement"
        else:
            revised_element[2] = revised_element[1]   # reuse bond as end
    rev = False                                       # then throw the flag away
    if len(revised_element) == 3 and revised_element[2] is not None:
        revised_element = [idx + self.start for idx in revised_element]
        self.complex.rotate_element(revised_element, angle, reverse=rev)
    elif len(revised_element) == 3 and revised_element[2] is None:
        ...
```

| Rule | What it breaks |
| --- | --- |
| P1 | `end` is `int` or `None`; `0` also means "complement" |
| P2 | takes chain-local, emits global, no type says so |
| P3 | copies defensively because the caller may share the list |
| P4 | reaches through `self.complex` to do the work |
| P5 | needs a built `Complex` to run at all |

Under the proposal this function does not exist. `reverse` is resolved to a
different `AtomRange` when the template loads. By the time anything rotates there
is one path:

```python
def rotate(self, torsion: Torsion, angle: float) -> Pose:
    """Rotate torsion.moving about the pivot-to-bond axis."""
    sl = torsion.moving.as_slice()
    axis = self.xyz[torsion.bond] - self.xyz[torsion.pivot]
    origin = self.xyz[torsion.pivot]
    return Pose(_rodrigues(self.xyz, sl, origin, axis, angle), self.system)
```

No flags. No sentinels. No frame question. No back-reference. Testable against a
4-atom `Pose` built by hand.

**Next:** [01-values.md](01-values.md).
</content>
