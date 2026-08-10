---
name: docstring
description: Use when writing or editing ANY docstring, comment, or module header in the maws package - covers required structure, signature lines, math directives, cross-references, units and array shapes, reusable argument text, and the checklist every docstring must pass before it is considered done.
---

# MAWS Docstring Writing Guide

Adapted from the PyTorch docstring guide
(`pytorch/.claude/skills/docstring/SKILL.md`). The discipline is theirs; the
format is NumPy, because that is what numpy, scipy, scikit-learn and pandas use
and what this package uses throughout.

The section mapping, since the source guide is written in Sphinx/Google style:

| PyTorch | Here |
| --- | --- |
| `Args:` | `Parameters` + `----------` underline |
| `Keyword args:` | folded into `Parameters`, after the `*` |
| `Returns:` | `Returns` + `-------` underline |
| `Examples::` | `Examples` + `--------` underline |
| (none) | `See Also`, `Raises`, `Warns` — used here |

Everything else — the signature first line, math directives, cross-reference
roles, admonitions, the checklist — carries over unchanged.

## General Principles

- Use **raw strings** (`r"""..."""`) for any docstring containing a backslash,
  so LaTeX math and escapes survive.
- Follow **NumPy docstring format** (numpydoc), parsed by Sphinx.
- Be **concise but complete** — every essential fact, nothing else. A docstring
  is a reference entry, not an essay.
- Always include **examples**. They are run by `pytest --doctest-modules`, so a
  wrong example is a failing test.
- Use **cross-references** to related functions and classes rather than
  re-explaining them.
- **Line length 88 columns**, docstrings included. `ruff` enforces it.

## Who is reading

Someone opening this package for the first time. They know Python, know nothing
else, and may not be a chemist. They have never seen an earlier version of this
code. Three rules follow, and they are the ones most often broken:

1. **Never mention code the reader cannot see.** No "replaces X", no "the old
   code", no "previously", no names from a design that no longer exists. Design
   history belongs in `docs/redesign/`.
2. **Define each domain term the first time a module uses it**, in one
   sentence: aptamer, residue, torsion, force field, implicit solvent, LEaP,
   prmtop, 5'/3'.
3. **Banned words**: upstream, downstream, "above"/"below" for modules or
   layers, higher-level, lower-level, and a bare "the caller". Name the thing.

## Docstring Structure

### 1. Signature (first line)

Start callables whose signature is not self-evident with the call signature:

```
grow_aptamer(system, *, energy, sampler, n_nucleotides=15, beta=0.01) -> Iterator[StepEvent]
```

- Include the name, the `*` separator before keyword-only arguments, defaults,
  and the return annotation.
- **No trailing period.**
- Follow it with a blank line, then the summary.
- Skip it for one- or two-argument methods and for properties, where it only
  repeats what is already visible.

### 2. Brief description

One line, ending in a full stop, fitting on one line.

- **Function or method**: third-person verb. `Return the rotation matrix.`
- **Class**: noun phrase. `A run of consecutive atom indices.`
- **Property**: `type : description`, no verb. `int : Number of atoms.`

An extended description follows after a blank line, but only when the summary
leaves a real question open. Two or three sentences at most.

### 3. Mathematical formulas

Use a math directive for anything that is genuinely mathematics. Do not spell a
formula out in prose.

```
.. math::
    R = I + \sin(\theta) K + (1 - \cos(\theta)) K^2
```

Inline: :math:`\beta`. The docstring must be a raw string, or the backslashes
are read as Python escapes.

### 4. Cross-references

| Role | Links to |
| --- | --- |
| :class:`~maws.pose.Pose` | a class; `~` shows only the last part |
| :func:`maws.scoring.entropy_score` | a function |
| :meth:`~Pose.rotate` | a method |
| :attr:`n_atoms` | an attribute |
| :data:`VDW_RADII` | a module-level constant |
| ``literal`` | code, filenames, flags |

Plus a `See Also` section naming the two or three objects a reader wants next:

```
See Also
--------
Pose.rotate_all : Turn several bonds in one call.
maws.sampling.SurfaceSampler : Produces the placements this consumes.
```

### 5. Notes and warnings

`Notes` is for facts that change how the thing is used: units, a sign
convention, an ordering requirement, a performance trap, a numerical caveat.
Not for design commentary.

Use admonitions for anything a reader must not miss:

```
.. note::
    Building one of these starts an OpenMM simulation, which costs far more
    than a single evaluation. Build it once and reuse it.

.. warning::
    Positions are not checked against the topology. A pose from another
    system silently produces meaningless energies.
```

### 6. Parameters section

```
Parameters
----------
name : type, default=value
    What the value means, and what changes if you pick a different one.
```

Formatting rules:

- Type on the same line, after ` : `. Not in parentheses.
- `, default=value` when there is a default worth quoting; `, optional` when
  there is not. Never both.
- Combine parameters sharing a type and a meaning: `first, second : AtomRange`.
- Enumerated options in braces: `direction : {"3prime", "5prime"}`.
- Double backticks for inline code and literal values.
- Continuation lines indented four spaces under the name.
- **Never restate the type.** `probe : float` followed by "the probe radius" is
  a wasted line.

```
# no
reach : float, default=10.0
    The reach.

# yes
reach : float, default=10.0
    How far past the target's furthest atom the sampling region extends, in
    ångström. Larger values consider positions further from the target, at
    the cost of more proposals landing in empty space.
```

### 7. Keyword-only arguments

numpydoc has no separate section for these. Document them in `Parameters`, in
signature order, so the ones after the `*` come last. If the distinction
matters to the reader, the signature line in section 1 already shows it.

### 8. Returns section

```
Returns
-------
type
    What it is.
```

Name the value only when there is more than one, or when the docstring refers
to it by name later. Generators use `Yields` with the same shape.

`Raises` and `Warns` follow the same layout, one entry per exception or warning
type, describing the condition that triggers it. Every exception the body
raises must appear.

### 9. Examples section

Mandatory on every public function, method and class.

```
Examples
--------
>>> span = AtomRange(10, 25)
>>> len(span)
15
```

Formatting rules:

- `>>>` prompt for input, bare lines for output.
- **Show the output.** An example with no output demonstrates nothing.
- Runnable from a clean interpreter: include the imports it needs.
- Prefer a case whose answer can be checked by hand over a realistic but opaque
  one.
- Comments with `#` where a line needs one.
- Anything needing AmberTools or a large input gets `# doctest: +SKIP`.

### 10. External references

Cite the source of a constant, an algorithm, or a convention:

```
References
----------
.. [1] Bondi, A. (1964). "van der Waals Volumes and Radii".
       J. Phys. Chem. 68(3), 441-451.
```

Refer to it in the text as [1]_.

## Object Types

### Modules

```
"""
maws.modulename
===============

One-line statement of what this module is for.

Two or three paragraphs: the problem it solves, the main types, how they fit
together. Define every domain term used.

Examples
--------
>>> the shortest useful thing this module can do
"""
```

### Value types (frozen dataclasses)

Document the fields under `Parameters`, not `Attributes` — for a dataclass they
are the constructor arguments. Add `Attributes` only for something computed.

### Protocols

Say what a class must do to satisfy it, and name the implementations under
`See Also` — at least one real and one test stand-in.

### Properties

One line, `type : description`. No `Parameters`, no `Returns`.

### Module-level constants

A string directly beneath the assignment. Say where the value came from when it
is not obvious.

```python
DEFAULT_VDW = 1.70
"""Radius in ångström assumed for an element not listed in :data:`VDW_RADII`.

The value for carbon, the most common element in these molecules.
"""
```

### Tests

A test's docstring states the behaviour being pinned down, in one line, in the
present tense: `Rotating by a full turn returns every atom to where it was.`
Not `test rotate`.

## Common Patterns

### Array shapes and units

Every array parameter states its shape; every physical quantity states its
unit. Shapes go in double backticks, or in math notation when they carry
symbols:

```
xyz : numpy.ndarray
    Shape ``(N, 3)`` positions in ångström, one row per atom.
masses : numpy.ndarray
    Shape ``(N,)`` atomic masses in daltons.
```

Units used in this package: **ångström** for length, **radians** for angles,
**kJ/mol** for energy, **daltons** for mass, **mol/L** for concentration.

### Reusable argument text

Several arguments recur across the package. Word them the same way every time,
so a reader who has met one has met them all:

```
rng : numpy.random.Generator, optional
    Source of randomness. Pass one built with a fixed seed to make a run
    repeatable. Defaults to a fresh generator, so runs differ.

pose : maws.pose.Pose
    The atom positions to work from.

beta : float, default=0.01
    How sharply lower energies are favoured, in mol/kJ.

probe : float, default=1.4
    Radius in ångström of the ball rolled over the target to find its
    surface. 1.4 Å is the size of a water molecule, the usual choice.

max_iterations : int, default=100
    How many adjustment steps to allow before stopping.
```

### Immutability

Every method returning a new object instead of mutating says so, in these
words, as the last line of `Returns`:

```
Returns
-------
Pose
    A new pose. This one is unchanged.
```

## Comments

Comments explain **why**, never **what**. A comment restating the line below it
is noise.

```python
# no
# add 1 to index
index += 1

# yes
# tleap numbers residues from 1, the arrays here from 0.
index += 1
```

## Complete Example

```python
def entropy_score(energies: Sequence[float], *, beta: float = 0.01) -> float:
    r"""entropy_score(energies, *, beta=0.01) -> float

    Return how concentrated a candidate's energies are. Lower is better.

    A strand binds well when it settles into a small family of shapes rather
    than wandering between many. This measures that concentration: the
    energies are turned into probabilities and the spread of that distribution
    is returned.

    Parameters
    ----------
    energies : sequence of float
        Energies of the shapes tried for one candidate strand, in kJ/mol. At
        least one is required.
    beta : float, default=0.01
        How sharply lower energies are favoured, in mol/kJ. Raising it makes
        the score depend mostly on the few best shapes; at zero every shape
        weighs the same and the score is always zero.

    Returns
    -------
    float
        Zero when every shape is equally likely, and increasingly negative as
        the weight concentrates onto fewer shapes.

    Raises
    ------
    maws.errors.ConfigurationError
        If `energies` is empty, or `beta` is negative.

    See Also
    --------
    boltzmann_weights : The per-shape probabilities this is computed from.

    Notes
    -----
    Each shape *i* is weighted by its Boltzmann factor, normalised to a
    probability, and the result is the negative entropy of that distribution:

    .. math::
        p_i = \frac{e^{-\beta E_i}}{\sum_j e^{-\beta E_j}}
        \qquad
        S = -\sum_i p_i \ln(p_i N)

    The factor :math:`N` inside the logarithm puts the zero point at "every
    shape equally likely" whatever *N* is, so runs that sampled different
    numbers of shapes stay comparable.

    Sums are evaluated at 60 decimal digits. Whole-molecule energies run to
    thousands of kJ/mol, and :math:`e^{-\beta E}` for those overflows or
    underflows ordinary floating point.

    References
    ----------
    .. [1] Kalinowski, M. et al. (2016). "MAWS - Making Aptamers Without
           SELEX". iGEM Heidelberg.

    Examples
    --------
    Ten shapes, one far better than the rest, against ten that are much of a
    muchness. The first scores lower, meaning more promising.

    >>> one_clear_winner = entropy_score([0.0] + [1000.0] * 9)
    >>> nothing_to_choose = entropy_score([0.0] + [10.0] * 9)
    >>> one_clear_winner < nothing_to_choose
    True

    Equal energies score exactly zero.

    >>> abs(entropy_score([100.0, 100.0, 100.0])) < 1e-12
    True
    """
```

## Quick Checklist

When writing a MAWS docstring, ensure:

- [ ] Raw string (`r"""`) if it contains a backslash
- [ ] Signature first line, where the signature is not self-evident
- [ ] Summary line: one line, correct form for its kind, ends in a full stop
- [ ] No reference to code the reader cannot see
- [ ] Every domain term used is defined somewhere in the module
- [ ] No banned words (upstream, downstream, above, below, higher/lower-level)
- [ ] Every parameter documented, with type and `default=` where applicable
- [ ] Every parameter description says what the value *controls*, not its type
- [ ] Array shapes given; units given for every physical quantity
- [ ] `Returns` documented; "A new pose. This one is unchanged." where it applies
- [ ] Every exception the body raises appears under `Raises`
- [ ] `See Also` points at the obvious next thing to read
- [ ] Math in a `.. math::` directive
- [ ] Admonition for anything a reader must not miss
- [ ] At least one `Examples` block, showing its output
- [ ] The example runs and produces exactly that output
- [ ] Every line within 88 columns

## Common Sphinx Roles Reference

- :class:`~maws.pose.Pose` — class reference
- :func:`maws.scoring.entropy_score` — function reference
- :meth:`~Pose.rotate` — method reference
- :attr:`n_atoms` — attribute reference
- :data:`VDW_RADII` — module constant reference
- :math:`\beta` — inline math
- ``code`` — inline literal

## Additional Notes

- **Indentation**: four spaces for parameter descriptions and code.
- **Line length**: 88 columns, enforced by `ruff`.
- **Periods**: end sentences with them; not the signature line.
- **Backticks**: double for code and literals, single for a parameter named in
  prose.
- **Types**: prefer the name a reader can import — `numpy.ndarray`, `Pose`,
  `maws.values.AtomRange` — over an abstract description.
