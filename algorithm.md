# Step 1 — High-level overview (what MAWS2023.py does and why)

Think of MAWS as a **greedy sequence builder** for an aptamer around a ligand. It grows one nucleotide at a time (front or back), and for each candidate it **samples a bunch of poses**, scores them using a simple physics energy (OpenMM) + a **Boltzmann-style entropy** term, and keeps whichever nucleotide choice gives the “best” (lowest) entropy. Rinse, repeat, until you hit your target length.

Below is the mental model first, then where each variable fits.

---

## The goal (one sentence)

Given a ligand structure (PDB), **build an RNA/DNA sequence** that is likely to bind it by iteratively adding one nucleotide at a time and choosing the option that looks most favorable under quick-and-dirty sampling.

---

## Inputs & Outputs

**Inputs**

- Ligand 3D structure: `-p /path/to/ligand.pdb`
- Polymer type: `RNA` or `DNA` (chooses residue templates + force field)
- Ligand type: `protein`, `organic`, or `lipid` (chooses force field)
- Target aptamer length: `-nt` (e.g., 15)
- Sampling budget knobs: `-c1` (first step) and `-c2` (subsequent steps)
- Entropy temperature knob: `-b` (beta)

**Outputs (files)**

- Per-candidate snapshots: `MAWS_aptamer_1_G.pdb` … etc.
- Step-best snapshots: `MAWS_aptamer_best_k_X.pdb`
- Final structure: `MAWS_aptamer_RESULT.pdb`
- Logs: `MAWS_aptamer_output.log`, `MAWS_aptamer_entropy.log`, `MAWS_aptamer_step_cache.pdb`

---

## Big picture pipeline

1. **Setup chemistry**

   - Load aptamer residue templates (`RNA.xml` / `DNA.xml`) and select Amber force fields for (i) aptamer and (ii) ligand.
   - Build a **template complex** with an empty aptamer chain plus the ligand (from PDB). This template is deep-copied per trial.

2. **Define a sampling space**

   - Place a cube (20 Å wide) centered at the ligand’s center of mass.
   - Also define a small “rotation vector” space for internal torsions of nucleotides (4 angles = `N_ELEMENTS`).

3. **Seed (first nucleotide)**

   - Try each letter (RNA: G, A, U, C; DNA: G, A, T, C). For each:

     - Build the complex with that 1-mer aptamer.
     - **Sample** many random rigid-body placements inside the cube, plus internal rotations.
     - For each pose: get energy via OpenMM; keep the **lowest energy pose** as that letter’s representative; collect the full energy list for **entropy**.
     - Convert that list of energies into a **Boltzmann probability distribution** using β, compute the **entropy S** (lower is “sharper”/more confident).

   - Pick the nucleotide with the **lowest S** as the first residue.

4. **Grow the chain**

   - For each remaining position (until `-nt` nucleotides):

     - Try adding each letter **either at the front or the back** (prepend/append).
     - Rebuild the complex, lightly “shake”/minimize it, sample internal torsions, record energies, compute entropy S.
     - Keep the choice with the **lowest S**; carry its coordinates forward to the next step.

5. **Write outputs**

   - Dump per-step best structure, and finally `*_RESULT.pdb` with the best full sequence and coordinates.

---

## Why the score looks like “entropy of energies”

- For each candidate (e.g., “add A at the back”), we generate a **set of energies** from random conformations/rotations.
- We transform energies $E_i$ into **Boltzmann weights**:
  $p_i = \frac{e^{-\beta E_i}}{\sum_j e^{-\beta E_j}}$
- Then compute **Shannon entropy**: $S = -\sum_i p_i \log p_i$.
- **Lower entropy** means the distribution is more **peaked**: a lot of random trials end up favoring a **few very good poses**. That’s desirable: it suggests the nucleotide choice tends to land in good interactions rather than only rarely hitting a good pose.
- β (`-b`) acts like **inverse temperature**: higher β = more sensitive to energy differences (more peaky), lower β = flatter, more tolerant exploration.

> tl;dr: we reward candidates whose sampled pose-energies consistently cluster near low values (not just a single lucky low-energy outlier).

---

## Key variables and what they control

- `N_NTIDES`: how many nucleotides to build in total.
- `nt_list`: the alphabet to try per position (RNA: `"GAUC"`, DNA: `"GATC"`).
- `FIRST_CHUNK_SIZE` (`-c1`): number of random samples for the **first** letter. Bigger = more robust but slower.
- `SECOND_CHUNK_SIZE` (`-c2`): samples per candidate in later steps. This scales linearly with 2×alphabet size × (N_NTIDES-1).
- `BETA` (`-b`): inverse temperature for Boltzmann weighting. Start with \~0.01; tweak if you find entropy isn’t discriminating enough.
- `N_ELEMENTS = 4`: number of **internal torsional angles** sampled per iteration for nucleotides (e.g., sugar-phosphate torsions). This matches what the `Chain.rotate_in_residue` expects.
- `cube = Space.Cube(20.0, COM)`: rigid-body search volume for the current aptamer relative to the ligand (translate + global rotation).
- `Complex`, `Chain`: your “system” object and the mutable aptamer chain inside it.
- `get_energy()`: delegates to OpenMM to get the **potential energy** (kJ/mol) for the current coordinates.

---

## Pseudocode (so you can “see” the loop)

```text
load XML + choose force fields
template = Complex(aptamerFF, ligandFF)
template.add_chain("", aptamerXML)          # empty aptamer chain
template.add_chain_from_PDB(ligand.pdb)     # ligand chain

# separate complex to compute ligand center for sampling
c = Complex(...)
c.add_chain_from_PDB(ligand.pdb)
c.build()
cube_center = centerOfMass(c.positions)

best_sequence, best_positions = None, None

# Step 1: pick first nucleotide
for nt in nt_list:
    cx = deepcopy(template)
    cx.chains[0].create_sequence(nt)
    cx.build()
    energies = []
    best_pose_for_nt = None

    for t in range(FIRST_CHUNK_SIZE):
        random rigid-body move in cube
        random internal torsions (N_ELEMENTS)
        E = cx.get_energy()
        track lowest E pose
        append E to energies
        reset positions

    S = entropy(energies, beta=BETA)
    keep nt with min S (record positions/topology)

# Steps 2..N_NTIDES: prepend or append each letter; pick min-entropy choice
for k in 2..N_NTIDES:
    for nt in nt_list:
        for append_or_prepend in [True, False]:
            cx = deepcopy(template)
            cx.chains[0].create_sequence(best_sequence)
            cx.build()
            cx.positions = best_positions
            add nt (append or prepend)
            cx.rebuild()
            cx.pert_min()

            energies = []
            for t in range(SECOND_CHUNK_SIZE):
                random internal torsions
                E = cx.get_energy()
                track lowest E pose
                append E
                reset positions

            S = entropy(energies, beta=BETA)
            keep mutation with min S

# write RESULT.pdb with final sequence + best_positions
```

---

## What to tweak first when you run it

- **Speed vs. quality**: lower `-c1/-c2` to get a quick run; raise for better stability.
- **Temperature sensitivity**: if candidates look indistinguishable, try increasing `-b` slightly; if it’s too sensitive (always picks the same letter), reduce `-b`.
- **Cube size**: 20 Å is conservative around the ligand COM; if your ligand is huge or buried, you might go smaller/larger.
