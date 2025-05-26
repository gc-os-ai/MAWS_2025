"""
write_structure_xml.py – serialize a `structure.Structure` instance to XML.

Usage
-----
>>> import structure, load_from
>>> dna = load_from.xml_structure("DNA.xml")   # read existing
>>> # ... modify `dna` ...
>>> write_structure_xml(dna, "DNA_new.xml")
"""

from __future__ import annotations

from pathlib import Path
from typing import Sequence

from lxml import etree

import MAWS.src.Structure as structure


# --------------------------------------------------------------------------- #
# XML helpers
# --------------------------------------------------------------------------- #
def _add_child(parent: etree._Element, tag: str, text: str) -> None:
    child = etree.SubElement(parent, tag)
    child.text = text


# --------------------------------------------------------------------------- #
# Public writer
# --------------------------------------------------------------------------- #
def write_structure_xml(struct: structure.Structure, out_path: str | Path) -> None:
    """
    Dump *struct* to ``out_path`` in the same dialect understood by
    `load_from.xml_structure`.
    """
    root = etree.Element("structure")

    # Residue path (optional)
    if struct.residue_path is not None:
        _add_child(root, "residuePath", str(struct.residue_path))

    # ------------------------------------------------------------------ #
    # Alias block
    # ------------------------------------------------------------------ #
    alias_block = etree.SubElement(root, "alias")
    for name, (alone, start, middle, end) in struct.alias.items():
        elem = etree.SubElement(alias_block, "element", name=name)
        _add_child(elem, "alone",   alone)
        _add_child(elem, "start",   start)
        _add_child(elem, "middle",  middle)
        _add_child(elem, "end",     end)

    # ------------------------------------------------------------------ #
    # Residues
    # ------------------------------------------------------------------ #
    for res in struct.residue_names:
        res_elem = etree.SubElement(
            root,
            "residue",
            name=res,
            length=str(struct.residue_lengths[res]),
        )

        # append / prepend connect
        app, prep, bl_app, bl_prep = struct.connect[res]
        app_elem = etree.SubElement(res_elem, "append")
        _add_child(app_elem, "newFirstAtom", str(app[0]))
        _add_child(app_elem, "oldLastAtom",  str(app[1]))
        _add_child(app_elem, "bondLength",   str(bl_app))

        prep_elem = etree.SubElement(res_elem, "prepend")
        _add_child(prep_elem, "newLastAtom",  str(prep[0]))
        _add_child(prep_elem, "oldFirstAtom", str(prep[1]))
        _add_child(prep_elem, "bondLength",   str(bl_prep))

        # backbone
        bb = struct.backbone_elements.get(res)
        if bb:
            # by construction bb has exactly one entry
            ((s, mp, bond), (mpost, end)) = bb
            bb_elem = etree.SubElement(res_elem, "backbone")
            for tag, val in zip(
                ("start", "middle_pre", "bond", "middle_post", "end"),
                (s, mp, bond, mpost, end),
            ):
                _add_child(bb_elem, tag, str(val))

        # rotations
        for idx, (start, bond, end) in enumerate(struct.rotating_elements[res]):
            rot = etree.SubElement(res_elem, "rotation", index=str(idx))
            _add_child(rot, "start", str(start))
            _add_child(rot, "bond",  str(bond))
            _add_child(rot, "end",   "end" if end is None else str(end))

    # pretty-print & write
    xml_bytes = etree.tostring(root, pretty_print=True, encoding="utf-8")
    Path(out_path).write_bytes(xml_bytes)


# --------------------------------------------------------------------------- #
# Legacy alias (camelCase)
# --------------------------------------------------------------------------- #
WriteStructureXML = write_structure_xml
