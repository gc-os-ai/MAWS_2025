"""
MAWS - Making Aptamers With Software
=====================================

A Python package for computational aptamer design.

Public API
----------
run_maws : Run the MAWS algorithm.
MAWSConfig : Configuration for a MAWS run.
MAWSResult : Result of a MAWS run.

Examples
--------
>>> from maws import run_maws, MAWSConfig
>>> config = MAWSConfig(
...     pdb_path="data/ligand.pdb",
...     num_nucleotides=15,
...     aptamer_type="RNA",
... )
>>> result = run_maws(config)  # doctest: +SKIP
"""

from maws.run import MAWSConfig, MAWSResult, run_maws

__all__ = ["run_maws", "MAWSConfig", "MAWSResult"]
