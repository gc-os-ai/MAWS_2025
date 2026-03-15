#!/usr/bin/env bash

# run_gpu_protein_experiment.sh
# 
# A helper script to quickly kick off a protein targeted MAWS experiment on a GPU cluster.
# Before running this script, ensure you have allocated a GPU node and activated the `maws` conda environment.

# Exit immediately if a command exits with a non-zero status
set -e

# Define Experiment Parameters (Can be overridden via environment variables)
TARGET_PDB=${TARGET_PDB:-"data/target_protein.pdb"}
EXPERIMENT_NAME=${EXPERIMENT_NAME:-"gpu_protein_experiment"}
MOLECULE_TYPE=${MOLECULE_TYPE:-"protein"}
NTIDES=${NTIDES:-15}
PLATFORM=${PLATFORM:-"CUDA"}
PLATFORM_PROPS=${PLATFORM_PROPS:-'{"Precision": "mixed", "DeviceIndex": "0"}'}

echo "=========================================="
echo "Starting MAWS GPU Experiment..."
echo "Target PDB: $TARGET_PDB"
echo "Experiment Name: $EXPERIMENT_NAME"
echo "Platform: $PLATFORM"
echo "Platform Properties: $PLATFORM_PROPS"
echo "=========================================="

python maws/maws2023.py \
    --name "$EXPERIMENT_NAME" \
    --path "$TARGET_PDB" \
    --moleculetype "$MOLECULE_TYPE" \
    --clean-pdb \
    --platform "$PLATFORM" \
    --platform-props "$PLATFORM_PROPS" \
    --ntides "$NTIDES"

echo "=========================================="
echo "Experiment completed."
echo "Results saved as ${EXPERIMENT_NAME}_RESULT.pdb"
echo "Logs saved as ${EXPERIMENT_NAME}_output.log"
echo "=========================================="
