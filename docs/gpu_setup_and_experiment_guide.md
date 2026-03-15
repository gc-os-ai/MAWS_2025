# GPU Cluster Setup and Experiment Guide (MAWS 2025)

This guide covers how to set up and run MAWS experiments on a GPU cluster, using explicit platform tuning to manage computational memory and precision efficiently for large protein targets.

## Cluster Setup

1. **Request GPU Resources:**
   Use your cluster's job scheduler (like SLURM) to allocate a GPU node.
   ```bash
   # Example SLURM command for an interactive session:
   srun --partition=gpu --gres=gpu:1 --time=24:00:00 --pty bash
   ```

2. **Load Required Modules:**
   Different clusters have different module names, but typically you need CUDA and a Conda environment manager.
   ```bash
   module load cuda/11.8
   module load miniconda3
   ```

3. **Install Dependencies:**
   Clone the repository and spin up the environment from the provided configuration.
   ```bash
   git clone https://github.com/gc-os-ai/MAWS_2025.git
   cd MAWS_2025
   conda env create -f environment.yml
   conda activate maws
   ```

## Running the Experiment

To leverage the GPU capabilities natively, MAWS now supports the `--platform` and `--platform-props` arguments. When dealing with protein PDB files, using `mixed` precision and explicit CUDA device indices prevents Out-of-Memory (OOM) errors and optimizes execution limits.

```bash
# Execute MAWS targeting the CUDA platform
python maws/maws2023.py \
    --name target_protein_experiment \
    --path data/target_protein.pdb \
    --moleculetype protein \
    --clean-pdb \
    --platform CUDA \
    --platform-props '{"Precision": "mixed", "DeviceIndex": "0"}' \
    --ntides 15
```

## Expected Outputs

- `target_protein_experiment_output.log`: Detailed logs containing parameters, configuration, and progression steps.
- `target_protein_experiment_entropy.log`: Final logged entropy and energy values per progression cycle.
- `target_protein_experiment_RESULT.pdb`: Handled and constructed PDB representation of the protein-aptamer complex.
