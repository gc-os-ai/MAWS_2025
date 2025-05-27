# MAWS (Modified)

> Modified version of [iGEM MAWS](https://github.com/iGEM-NU-Kazakhstan/MAWS-Heidelberg-x-NU_Kazakhstan) for easier execution and better documentation.

⚠️ This is a work in progress. Some components may be incomplete or require adaptation.

## What is MAWS?

MAWS is a molecular modeling toolkit for **de novo design of aptamers** (DNA/RNA) without SELEX. It uses **entropy-based scoring** to iteratively select nucleotides that are most likely to bind a target protein with high specificity and affinity.

## Core Capabilities

- Generates DNA or RNA aptamer sequences for specific protein targets
- Performs conformational sampling using global and internal transformations
- Calculates binding energies using OpenMM
- Applies entropy minimization for nucleotide selection
- Outputs 3D aptamer-protein complexes in `.pdb` format

## Conformational Sampling Strategy

MAWS explores aptamer conformations using two sampling spaces:

| Sampling Space           | Class           | Purpose                          | Parameters                             |
| ------------------------ | --------------- | -------------------------------- | -------------------------------------- |
| Translational/Rotational | `Space.cube`    | Global aptamer positioning       | 20 Å cube around protein center        |
| Internal Rotations       | `Space.Nangles` | Bond dihedral angle perturbation | N_ELEMENTS = number of rotatable bonds |

### Sampling Process

1. `cube.generator()` randomly translates and rotates the aptamer globally
2. `rotations.generator()` samples internal bond angle rotations
3. `aptamer.translate_global()` and `aptamer.rotate_global()` apply global transformations
4. `aptamer.rotate_in_residue()` applies internal conformational changes

## Entropy-Based Selection

The aptamer is built nucleotide-by-nucleotide by selecting the base with the **lowest entropy**, derived from energy distributions across sampled conformations.

### Algorithm Logic

- For each candidate base (G, A, T/U, C), sample conformations
- Compute energy distribution across conformers
- Use entropy formula: `entropy = S(energies, beta=BETA)`
- Select base with **minimum entropy**
- Append base to growing aptamer and repeat

This results in a sequence with high likelihood of favorable and stable binding.

## CLI Parameters and Configuration

You can run MAWS with the following command-line arguments:

| Parameter        | Flag(s)                    | Purpose                                 |
| ---------------- | -------------------------- | --------------------------------------- |
| Job Name         | `-n`, `--name`             | Prefix for output files                 |
| Beta             | `-b`, `--beta`             | Inverse temperature for entropy         |
| First Chunk Size | `-c1`, `--firstchunksize`  | Number of samples for first nucleotide  |
| Chunk Size       | `-c2`, `--secondchunksize` | Samples per nucleotide extension step   |
| Nucleotides      | `-t`, `--ntides`           | Total aptamer length                    |
| PDB Path         | `-p`, `--path`             | Input protein PDB file path (required)  |
| Aptamer Type     | `-a`, `--aptamer_type`     | Select aptamer type: `"DNA"` or `"RNA"` |

### Nucleotide Sets

- **DNA**: `["DGN", "DAN", "DTN", "DCN"]` (from `DNA.xml`)
- **RNA**: `["G", "A", "U", "C"]` (from `RNA.xml`)

## Output Files

| File Pattern             | Content                       | Purpose                       |
| ------------------------ | ----------------------------- | ----------------------------- |
| `{JOB_NAME}_output.log`  | Algorithm execution trace     | Execution tracking            |
| `{JOB_NAME}_entropy.log` | Sequence entropy calculations | Selection diagnostics         |
| `{JOB_NAME}_RESULT.pdb`  | Final aptamer-protein complex | Optimized 3D structure output |

The `.pdb` file contains the **final aptamer bound to the protein** with the most thermodynamically favorable sequence configuration.

---

## Usage Example

1. Clean the target PDB file using the script:

```bash
python MAWS/utils/pdb_cleaner.py
```

2. Then run the aptamer design pipeline:

```bash
python -m MAWS.src.maws_cli \
  --name quick_test \
  --pdb examples/1hut_cleaned_keep_one_chain.pdb \
  --aptamer-type DNA \
  --firstchunksize 200 \
  --secondchunksize 200 \
  --ntides 5
```

---

## To Do

- Integrate Biopython's `PDBParser` for internal cleanup
- Add automatic parameter detection and logging
- Build GUI or web interface for easier usability
- Validate and benchmark outputs against known aptamer structures

---

## Credits

- Original authors: [iGEM NU Kazakhstan Team](https://github.com/iGEM-NU-Kazakhstan)
- Modifications: Ongoing efforts for accessibility
