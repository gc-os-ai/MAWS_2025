"""
load_from.py – utilities to build :class:`structure.Structure` objects
from XML configuration files or from a PDB file alone.

The XML dialect is the same as the legacy MAWS repository used:
    <root>
        <residuePath>./amber_params</residuePath>
        <alias> … </alias>
        <residue name="ALA" length="10"> … </residue>
        …
    </root>
"""

from __future__ import annotations

import xml.etree.ElementTree as ET
from pathlib import Path
from typing import List, Sequence, Tuple

from openmm import app

import MAWS.src.Structure as structure

# --------------------------------------------------------------------------- #
# Type aliases
# --------------------------------------------------------------------------- #
RotationSpec   = Tuple[str, int, int, int | None]
BackboneSpec   = Tuple[str, int, int, int, int, int]
ConnectSpec    = Tuple[Tuple[int, int], Tuple[int, int], float, float]
AliasSpec      = List[str]                        # [name, alone, start, middle, end]
XMLParseReturn = Tuple[
    List[str],           # residue_names
    List[int],           # residue_lengths
    List[RotationSpec],  # rotations
    List[BackboneSpec],  # backbone
    List[ConnectSpec],   # connect
    str | None,          # residue_path
    List[AliasSpec],     # aliases
]


# --------------------------------------------------------------------------- #
# XML helpers
# --------------------------------------------------------------------------- #
def _parse_xml_structure(path: str | Path) -> XMLParseReturn:
    """Parse MAWS-style XML into raw arrays suitable for `structure.Structure`."""
    file_path = Path(path).expanduser().resolve()
    root = ET.parse(file_path).getroot()

    # ---- high-level containers ----------------------------------------- #
    residue_names: List[str] = []
    residue_lengths: List[int] = []
    rotations: List[RotationSpec] = []
    backbones: List[BackboneSpec] = []
    connects: List[ConnectSpec] = []
    aliases: List[AliasSpec] = []

    residue_path = root.findtext("residuePath", default=None)

    # ---- aliases -------------------------------------------------------- #
    for alias in root.findall("alias/element"):
        if len(alias) != 4:  # expect <alone> <start> <middle> <end>
            continue
        aliases.append(
            [
                alias.get("name", ""),
                alias.findtext("alone", ""),
                alias.findtext("start", ""),
                alias.findtext("middle", ""),
                alias.findtext("end", ""),
            ]
        )

    # ---- residues ------------------------------------------------------- #
    for res in root.findall("residue"):
        name = res.get("name")
        length = int(res.get("length", "0"))
        if name is None:
            raise ValueError("Residue without a 'name' attribute found in XML.")

        residue_names.append(name)
        residue_lengths.append(length)

        # rotations
        for rot in res.findall("rotation"):
            rotations.append(
                (
                    name,
                    int(rot.findtext("start")),
                    int(rot.findtext("bond")),
                    None if rot.findtext("end") == "end" else int(rot.findtext("end")),
                )
            )

        # backbone
        bb = res.find("backbone")
        if bb is None:
            raise ValueError(f"<backbone> missing for residue {name}")
        backbones.append(
            (
                name,
                int(bb.findtext("start")),
                int(bb.findtext("middle_pre")),
                int(bb.findtext("bond")),
                int(bb.findtext("middle_post")),
                int(bb.findtext("end")),
            )
        )

        # connect
        appnd = res.find("append")
        prep = res.find("prepend")
        if appnd is None or prep is None:
            raise ValueError(f"<append>/<prepend> missing for residue {name}")

        connect_spec: ConnectSpec = (
            (
                int(appnd.findtext("newFirstAtom")),
                int(appnd.findtext("oldLastAtom")),
            ),
            (
                int(prep.findtext("newLastAtom")),
                int(prep.findtext("oldFirstAtom")),
            ),
            float(appnd.findtext("bondLength")),
            float(prep.findtext("bondLength")),
        )
        connects.append(connect_spec)

    return (
        residue_names,
        residue_lengths,
        rotations,
        backbones,
        connects,
        residue_path,
        aliases,
    )


def xml_structure(path: str | Path) -> structure.Structure:
    """Return a ready :class:`structure.Structure` built from *path*."""
    parsed = _parse_xml_structure(path)
    return structure.Structure(
        parsed[0],                       # names
        parsed[1],                       # lengths
        rotating_elements=parsed[2],
        backbone_elements=parsed[3],
        connect=parsed[4],
        residue_path=parsed[5],
        alias=parsed[6],
    )


# --------------------------------------------------------------------------- #
# PDB helpers
# --------------------------------------------------------------------------- #
def _parse_pdb_structure(path: str | Path) -> Tuple[List[str], List[int]]:
    """Return residue names and lengths from a PDB file."""
    pdb = app.PDBFile(str(path))
    names: List[str] = []
    lengths: List[int] = []

    for res in pdb.topology.residues():
        if res.name not in names:
            names.append(res.name)
            lengths.append(sum(1 for _ in res.atoms()))
    return names, lengths


def pdb_structure(path: str | Path) -> structure.Structure:
    """Shortcut that calls :func:`_parse_pdb_structure` and builds a Structure."""
    names, lengths = _parse_pdb_structure(path)
    return structure.Structure(names, lengths)


# --------------------------------------------------------------------------- #
# Backward-compatibility camelCase wrappers
# --------------------------------------------------------------------------- #
XMLParseStructure = _parse_xml_structure
XMLStructure = xml_structure
PDBParseStructure = _parse_pdb_structure
PDBStructure = pdb_structure

__all__ = [
    # modern names
    "xml_structure",
    "pdb_structure",
    # raw parse helpers
    "_parse_xml_structure",
    "_parse_pdb_structure",
    # legacy camelCase
    "XMLParseStructure",
    "XMLStructure",
    "PDBParseStructure",
    "PDBStructure",
]
