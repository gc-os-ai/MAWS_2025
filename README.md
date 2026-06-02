Repository for refactoring the "MAWS 2023" repository by `dtu-denmark` for integration in the [pyaptamer](https://github.com/gc-os-ai/pyaptamer/) package

### Refactoring aims

* remove or replace binary dependencies
* remove or replace subprocess calls
* full python bindings

### Sampling region

MAWS samples the initial nucleotide pose around the ligand surface using
SAS-style rejection: candidate poses inside the protein bulk (within
`vdW + probe`, default `1.4 Å` water-equivalent) are skipped before the
energy evaluator sees them. The envelope is a single sphere auto-sized
from the ligand geometry (`radius = R_max + reach`) around its
mass-weighted centre of mass. Two CLI flags control behavior:

* `--reach FLOAT` — how far the envelope extends past the ligand's
  bounding radius, in Å (default `10.0`).
* `--probe FLOAT` — vdW probe radius for the SAS rejection, in Å
  (default `1.4`, water-equivalent).

The same options are available as keyword arguments on
`maws.run.MawsRunner` for programmatic use. An opt-in surface-following
sampling mode (`maws.space.make_sampler(..., mode="surface-following",
d_max=...)`) is also implemented for users who want accepted poses
concentrated near the molecular surface; see [docs/space.md](docs/space.md)
for the full API.

### Implicit-solvent salt screening

Energies are evaluated with the GB-OBC1 implicit solvent. The monovalent
salt concentration used for Debye–Hückel screening is configurable:

* `--salt-conc FLOAT` — monovalent salt concentration in mol/L (default
  `0.15`, ~physiological). Also available as the `salt_conc` keyword on
  `maws.run.MawsRunner` and `maws.complex.Complex`.

> **Behavior change:** earlier releases ran unscreened (effectively
> `0.0` mol/L). Because screening changes the GB energies that drive
> sequence selection, results will differ from prior versions unless you
> pass `--salt-conc 0`. The screening is monovalent only and does not
> model divalent ions such as Mg²⁺.

### original readme

https://github.com/gc-os-ai/MAWS_2025/blob/main/README_orig.md
