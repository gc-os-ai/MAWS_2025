"""
maws.io
=======

Everything in MAWS that touches the outside world.

Reading a file, writing one, and starting another program all happen here and
nowhere else. Keeping them in one place is what lets the rest of the package be
tested without AmberTools installed and without leaving anything behind on
disk.

Three modules:

:mod:`maws.io.tools`
    Runs the AmberTools programs and turns their failures into MAWS errors.
:mod:`maws.io.prepare`
    Works out the parameters describing a target molecule, so that a structure
    containing it can be built.
:mod:`maws.io.pdb_cleaner`
    Tidies a downloaded PDB file into something the structure builder will
    accept: picking chains, dropping unwanted records, fixing hydrogens.

See Also
--------
maws.build : Uses all three to turn a description into a structure.
"""

from __future__ import annotations

__all__: list[str] = []
