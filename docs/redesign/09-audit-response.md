# 09 — What the audit found, and what this branch does about it

[PR #45](https://github.com/gc-os-ai/MAWS_2025/pull/45) is an adversarial
review of the scientific pipeline as it stood on `main`. This page says, for
each of its findings, what happened to it here.

Read it before comparing a run on this branch against an earlier one. Several
of these change the answer, and one of them changes which nucleotide the search
picks at every step.

## The short version

Nine of the findings the layered rewrite made impossible: they were consequences
of state being shared between things that had no business sharing it, and there
is no shared state left to consume. Twelve more were real defects in the
algorithm or the chemistry, and each is fixed below with the test that now
holds it in place. Every one of those tests is new — none of them could have
been written against the old code, because the quantity they check was not
something the old code had a name for.

## Fixed by the rewrite, with nothing extra to do

| Finding | Why it cannot happen now |
|---|---|
| A1 — rotations applied the transpose of the intended matrix | `pose.rodrigues` is written once, tested against hand-worked quarter turns, and every rotation in the package goes through it. |
| A5 — the pose and the topology could disagree about atom counts | A `Pose` carries the system it belongs to and refuses positions that do not match it. |
| A7 — a candidate could be scored against another candidate's structure | A `Candidate` carries its own `BuiltSystem` and `Pose`. Nothing is shared and nothing is mutated. |
| A6 — global state carried between runs | There is no global state. |
| D6 — the same file cache served two different molecules | The cache key is the contents, not the name. |
| D4 — settings read from three places that could disagree | One dataclass, `ForceField`, built once and passed down. |
| B5 — the sampler's random stream was entangled with the search's | Every source of randomness takes a `numpy.random.Generator`. |
| F1 — the aptamer could be built into the target's coordinates | Chains have spans and a chain's atoms are the only ones its moves touch. |

## Fixed here, with the test that pins each one

### The geometry

**B1 — a placement missed the position it was given, by up to 3.2 Å.**
Placing a chain slid it so its *first atom* landed on the proposed position and
then turned it about that atom, which moved every other atom away again. It
now turns about the chain's own centre and slides that centre onto the
position, so the offset is ~1e-15 Å whatever the angle.
`tests/unit/test_pose.py::TestPosePlace`

**§2.1 — turning a bond towards the 5' end could tear the molecule apart.**
The 5' reading of a bond swings "everything before it", which is only a
well-defined half of the molecule when cutting the bond would fall it into two
pieces. For a bond inside a ring it is not. Those bonds now read the same in
both directions.
`tests/integration/test_real_geometry.py::TestTurningOneBond`

**§2.2 — a bond that moves nothing was still handed to the search.**
At the 5' end of a strand the first bond has only the strand's own first atom
left to move, and that atom lies on the bond's own axis. `torsions()` leaves
such bonds out.
`tests/unit/test_geometry_invariants.py::TestTorsionsMoveSomething`

**A3 — growing at the 5' end explored less than growing at the 3' end.**
The 3' end got the new residue's own bonds plus the bond joining it to its
neighbour; the 5' end got only its own, one of which moved nothing. The search
picks an end by comparing the two scores, so the 5' end was losing for a reason
that had nothing to do with the molecule.
`tests/algorithm/test_search.py::TestDegreesOfFreedom`

**A4 — the DNA end residues turned about something that is not a bond.**
For the eight DNA residues that sit at the 3' end of a strand or stand alone,
the table pointed at a hydrogen and an oxygen 3.1 Å apart with no bond between
them. Turning about any line is still a rigid motion, so no bond changed length
and nothing caught it; what it did instead was wrench the bond angles.
`tests/integration/test_real_geometry.py::TestTheTurnableBondsAreBonds`

**§2.7 — splicing a regrown strand fitted it from the whole strand at once.**
A fresh build is fully extended while the shape being kept has been folded, so
no single rigid motion carries one onto the other and all the error landed on
the new residue. Fitting from the neighbouring residue alone puts the join at
0.53–1.48 Å, where it was 0.39–3.66.
`tests/unit/test_regrow.py::TestSplice`

### The sampling

**B2 — the clearance test looked at one point, not at the molecule.**
A placement says where the *middle* of the strand goes, and a nucleotide is
about ten ångström across, so a shape with half its atoms buried in the target
passed the test and was scored. `Sampler.accepts` is now asked once the atoms
are actually somewhere, and both sampling loops draw until they have as many
usable shapes as were asked for.
`tests/algorithm/test_search.py::TestShapesThatBuryTheStrandAreNeverScored`

### The selection criterion

**C1, C2 — the score preferred candidates that mostly clash.** This is the one
that changes which nucleotide is chosen.

`entropy_score` asks how *concentrated* the Boltzmann weight is across the
shapes tried. It divides each weight by their total, which throws away the scale
of the energies: adding the same amount to every energy leaves it untouched, so
a strand reaching −1000 kJ/mol and one that never gets below zero are
indistinguishable to it.

Concentration is also exactly what a clashing candidate produces. Sample a
thousand shapes, have 998 crash atoms into each other and two come out
ordinary, and the weight piles onto those two: score −6.21. A candidate whose
thousand shapes are all reasonable scores −0.13. Lower wins.

`free_energy_score` is the default instead:

<!-- markdownlint-disable-next-line MD033 -->
> F = −(1/β) · ln( (1/N) · Σᵢ exp(−β Eᵢ) )

Low energies push it down, and *many* low energies push it down further. It is
in kJ/mol, lies between the best energy and the average, and moves with the
energies rather than ignoring where they sit. `entropy_score` remains, for
comparing against published runs, carrying a warning that says what it does and
does not answer.
`tests/unit/test_scoring.py::TestTheTwoScoresDisagree`

The field carrying the number is `Candidate.score`, `MawsResult.score` and
`AptamerDesigner.scores_` — renamed from `entropy`, since it is no longer one
and its unit now depends on which scorer produced it.

### The set-up

**D3 — settling moved the target.** Twice over: the nudge-and-settle round
displaced every atom in the structure at random, and settling then relaxed the
target's own internal strain, which is worth hundreds of kJ/mol and lands
differently on each candidate. `perturb_and_minimize` now takes the atoms that
may move, and `OpenMMEnergy` takes atoms to freeze — zero mass, which is how
OpenMM is told an atom is held in place. Frozen atoms still contribute to the
energy, so the number shifts by a constant that cancels out of every
comparison.
`tests/integration/test_real_build.py::TestHoldingTheTargetStill`

**D1 — `antechamber` was never told the target's charge.** It shares a stated
total out across the atoms and was using its default of zero, so for a charged
target every atom came out with the wrong charge. `net_charge` is carried from
`maws design --net-charge` and from `design(net_charge=...)` through to
`antechamber -nc`. It still defaults to zero, which is right for a neutral
molecule and wrong for everything else — the docstrings say so.

### The input

**E1–E9 — the PDB cleaner corrupted its targets.** The whole module is
rewritten; see the commit message and
[maws/io/pdb_cleaner.py](../../maws/io/pdb_cleaner.py) for the detail. In
short: it sorted the file to resolve alternate positions and never undid the
sort, which rewrites the covalent structure LEaP reads from line order; it
deleted every residue carrying an insertion code, which are in the middle of
chains by definition, along with any `TER` line that carried one; and it threw
away everything after the last `TER`, losing bound ligands and cofactors.

On the repo's own `data/1HAO.pdb` the old code kept 1966 of 2769 atoms, wrote
one `TER` where the file has three, and lost the 30-atom inhibitor. The new
code keeps 2620 — the 2769 less exactly the 149 waters — keeps all three `TER`
records, keeps the inhibitor, and keeps all 309 insertion-coded atoms.
`tests/unit/test_pdb_cleaner.py`

### The verification gap

**§7 — nothing checked that a turn left the molecule in one piece.** This was
the audit's central point, and it was right: every one of the geometry defects
above is invisible to an energy, a score and a ranking, all of which come back
as ordinary numbers from a broken structure.

`tests/integration/test_real_geometry.py` now turns every bond of a real
AmberTools-built strand, in both directions, at three angles, and checks every
bond length against the bond list AmberTools produced. It also checks that
every bond the residue tables declare is a bond, which is what caught A4.

Guessing which atoms are bonded from how close they are does not work here, and
that is worth knowing: a structure straight out of `tleap` has not been settled
and holds unbonded atoms 1.1 Å apart, which a guess calls a bond and then
reports as broken the moment anything moves.

## Deliberately not changed

**The residue tables are still tables.** The audit suggests deriving turnable
bonds from the built topology rather than writing atom indices by hand. That
would be better, and it is a larger change than this branch is making. The new
check that every declared bond is a real bond closes the gap that made the
hand-written numbers dangerous.

**`beta` is still 0.01 mol/kJ.** It is not a physical temperature — 1/β is
100 kJ/mol, about forty times room temperature — but it is a tuning knob whose
value was chosen against runs that exist, and changing it is a separate
decision from fixing what it multiplies.
