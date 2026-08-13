"""
maws.libraries
==============

Residue tables for the nucleic acids MAWS can design, plus the factories that
turn them into validated :class:`~maws.values.ResidueLibrary` objects.

MAWS designs an *aptamer*: a short strand of RNA or DNA that folds up against a
target molecule and sticks to it. A strand is built from *residues* — single
nucleotides — and a library says what each of those residues is made of: how
many atoms it has, which of its bonds can be turned, and how it joins to its
neighbours. Two libraries are provided, :func:`rna` and :func:`dna`, and a
design run uses one or the other, never both.

Each submodule is pure data. Nothing here reads a file, imports OpenMM, or
starts another program, so a library can be built and inspected in a test with
nothing installed. Both factories cache their result, so calling one repeatedly
costs nothing.

See Also
--------
maws.values.ResidueLibrary : The type both factories return.
maws.forcefield.ForceField : Chooses which of the two a run uses.

Examples
--------
>>> from maws.libraries import dna, rna
>>> lib = rna()
>>> lib["G5"].n_atoms
32
>>> lib.alias("G").start
'G5'

DNA writes ``T`` where RNA writes ``U``, and prefixes its residue names:

>>> dna().alias("T").middle
'DT'
"""

from maws.libraries.dna import dna
from maws.libraries.rna import rna

__all__ = ["dna", "rna"]
