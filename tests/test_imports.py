import os
import pytest

# simply import all top-level modules to ensure they load without errors
from maws import (
    chain,
    complex,
    dna_structure,
    helpers,
    kernels,
    maws2023,
    pdb_cleaner,
    prepare,
    rna_structure,
    routines,
    space,
    structure,
    tools,
)


def test_load_dna_structure():
    struct = dna_structure.load_dna_structure()
    # verify attributes from module
    assert hasattr(struct, "residue_names")
    assert struct.residue_names[0] == dna_structure.RESIDUE_NAMES[0]


def test_load_rna_structure():
    struct = rna_structure.load_rna_structure()
    assert hasattr(struct, "residue_names")


def test_maws2023_resolve_path(tmp_path, capsys):
    # create a dummy pdb file
    fake = tmp_path / "dummy.pdb"
    fake.write_text("HEADER\n")
    out, orig = maws2023._resolve_pdb_path(str(fake), "protein", clean_pdb=False, keep_chains="all", remove_h=False, drop_hetatm=False)
    assert out == str(fake)
    assert orig == str(fake)

