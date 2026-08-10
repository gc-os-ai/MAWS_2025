"""
maws.libraries
==============

Residue tables for the nucleic acids MAWS can design, plus the factories that
turn them into validated :class:`~maws.values.ResidueLibrary` objects.

Each submodule is pure data. Nothing here reads a file, imports OpenMM, or
shells out, so a library can be built and inspected in a test with no
AmberTools installed.

Examples
--------
>>> from maws.libraries import rna
>>> lib = rna()
>>> lib["G5"].n_atoms
32
>>> lib.alias("G").start
'G5'
"""

from maws.libraries.dna import dna
from maws.libraries.rna import rna

__all__ = ["dna", "rna"]
