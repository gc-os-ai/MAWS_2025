Repository for refactoring the "MAWS 2023" repository by `dtu-denmark` for integration in the [pyaptamer](https://github.com/gc-os-ai/pyaptamer/) package

### Refactoring aims

* remove or replace binary dependencies
* remove or replace subprocess calls
* full python bindings

### Sampling region

MAWS samples the initial nucleotide pose around the ligand surface using
SAS-style rejection: candidate poses inside the protein bulk (within
`vdW + probe`, default `1.4 Å` water-equivalent) are skipped before the
energy evaluator sees them. The envelope is auto-sized from the ligand
geometry around its mass-weighted centre of mass. Three CLI flags control
behavior:

* `--shape {cube,sphere,shell}` — envelope geometry (default `shell`,
  which hugs the surface most tightly for both globular and elongated
  ligands).
* `--reach FLOAT` — how far the envelope extends beyond the ligand
  surface, in Å (default `10.0`).
* `--probe FLOAT` — vdW probe radius for the SAS rejection, in Å
  (default `1.4`, water-equivalent).

The same options are available as keyword arguments on
`maws.run.MawsRunner` for programmatic use.

### original readme

https://github.com/gc-os-ai/MAWS_2025/blob/main/README_orig.md
