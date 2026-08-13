import maws.structure as structure

# ---------------------------------------------------------------------------
# Core tables (mirrors the order & values in DNA.xml)
# ---------------------------------------------------------------------------

# 1) Residue names (must align with all per-residue lists below)
RESIDUE_NAMES: list[str] = [
    "DGN",
    "DAN",
    "DTN",
    "DCN",
    "DG",
    "DA",
    "DT",
    "DC",
    "DG5",
    "DA5",
    "DT5",
    "DC5",
    "DG3",
    "DA3",
    "DT3",
    "DC3",
]

# 2) Atom counts per residue (length=… in XML), aligned with RESIDUE_NAMES
RESIDUE_LENGTH: list[int] = [
    32,
    31,
    31,
    29,
    33,
    32,
    32,
    30,
    31,
    30,
    30,
    28,
    34,
    33,
    33,
    31,
]

# 3) Backbone tuples: (residue, start, middle_pre, bond, middle_post, end)
BACKBONE: list[tuple[str, int, int, int, int, int]] = [
    ("DGN", 0, 8, 10, 25, 31),
    ("DAN", 0, 8, 10, 24, 30),
    ("DTN", 0, 8, 10, 24, 30),
    ("DCN", 0, 8, 10, 22, 28),
    ("DG", 0, 10, 12, 27, 32),
    ("DA", 0, 10, 12, 26, 31),
    ("DT", 0, 10, 12, 26, 31),
    ("DC", 0, 10, 12, 24, 29),
    ("DG5", 0, 8, 10, 25, 30),
    ("DA5", 0, 8, 10, 24, 29),
    ("DT5", 0, 8, 10, 24, 29),
    ("DC5", 0, 8, 10, 22, 27),
    ("DG3", 0, 10, 12, 27, 33),
    ("DA3", 0, 10, 12, 26, 32),
    ("DT3", 0, 10, 12, 26, 32),
    ("DC3", 0, 10, 12, 24, 30),
]

# 4) Rotations: (residue, first_atom, second_atom) naming the bond to turn.
# Negative indices count back from the end of the residue and are
# resolved at runtime against RESIDUE_LENGTH.
ROTATIONS: list[tuple[str, int, int]] = [
    # DGN
    ("DGN", 0, 1),
    ("DGN", 1, 2),
    ("DGN", 8, 10),
    ("DGN", -7, -2),
    # DAN
    ("DAN", 0, 1),
    ("DAN", 1, 2),
    ("DAN", 8, 10),
    ("DAN", -7, -2),
    # DTN
    ("DTN", 0, 1),
    ("DTN", 1, 2),
    ("DTN", 8, 10),
    ("DTN", -7, -2),
    # DCN
    ("DCN", 0, 1),
    ("DCN", 1, 2),
    ("DCN", 8, 10),
    ("DCN", -7, -2),
    # DG
    ("DG", 0, 3),
    ("DG", 3, 4),
    ("DG", 10, 12),
    ("DG", -6, -1),
    # DA
    ("DA", 0, 3),
    ("DA", 3, 4),
    ("DA", 10, 12),
    ("DA", -6, -1),
    # DT
    ("DT", 0, 3),
    ("DT", 3, 4),
    ("DT", 10, 12),
    ("DT", -6, -1),
    # DC
    ("DC", 0, 3),
    ("DC", 3, 4),
    ("DC", 10, 12),
    ("DC", -6, -1),
    # DG5
    ("DG5", 0, 1),
    ("DG5", 1, 2),
    ("DG5", 8, 10),
    ("DG5", -6, -1),
    # DA5
    ("DA5", 0, 1),
    ("DA5", 1, 2),
    ("DA5", 8, 10),
    ("DA5", -6, -1),
    # DT5
    ("DT5", 0, 1),
    ("DT5", 1, 2),
    ("DT5", 8, 10),
    ("DT5", -6, -1),
    # DC5
    ("DC5", 0, 1),
    ("DC5", 1, 2),
    ("DC5", 8, 10),
    ("DC5", -6, -1),
    # DG3
    ("DG3", 0, 3),
    ("DG3", 3, 4),
    ("DG3", 10, 12),
    ("DG3", -7, -2),
    # DA3
    ("DA3", 0, 3),
    ("DA3", 3, 4),
    ("DA3", 10, 12),
    ("DA3", -7, -2),
    # DT3
    ("DT3", 0, 3),
    ("DT3", 3, 4),
    ("DT3", 10, 12),
    ("DT3", -7, -2),
    # DC3
    ("DC3", 0, 3),
    ("DC3", 3, 4),
    ("DC3", 10, 12),
    ("DC3", -7, -2),
]

# 5) Connectivity: per residue:
# [[append_first, append_last], [prepend_last, prepend_first], append_len, prepend_len]
# In this XML in previous verisons"
# all residues share the same connectivity (0,-1) and (-2,0) with bond lengths 1.6.
CONNECT: list[list] = [[[0, -1], [-2, 0], 1.6, 1.6] for _ in RESIDUE_NAMES]

# 6) Alias table: [name, alone, start, middle, end]
ALIASES: list[list[str]] = [
    ["DCN", "DCN", "DC5", "DC", "DC3"],
    ["A", "DAN", "DA5", "DA", "DA3"],
    ["C", "DCN", "DC5", "DC", "DC3"],
    ["DTN", "DTN", "DT5", "DT", "DT3"],
    ["G", "DGN", "DG5", "DG", "DG3"],
    ["DG3", "DG3", "DG", "DG", "DG3"],
    ["DG", "DG", "DG", "DG", "DG"],
    ["DAN", "DAN", "DA5", "DA", "DA3"],
    ["DA3", "DA3", "DA", "DA", "DA3"],
    ["DGN", "DGN", "DG5", "DG", "DG3"],
    ["DC", "DC", "DC", "DC", "DC"],
    ["DA", "DA", "DA", "DA", "DA"],
    ["DA5", "DA5", "DA5", "DA", "DA"],
    ["T", "DTN", "DT5", "DT", "DT3"],
    ["DG5", "DG5", "DG5", "DG", "DG"],
    ["DT3", "DT3", "DT", "DT", "DT3"],
    ["DT", "DT", "DT", "DT", "DT"],
    ["DC5", "DC5", "DC5", "DC", "DC"],
    ["DC3", "DC3", "DC", "DC", "DC3"],
    ["DT5", "DT5", "DT5", "DT", "DT"],
]

__all__ = [
    "RESIDUE_NAMES",
    "RESIDUE_LENGTH",
    "BACKBONE",
    "ROTATIONS",
    "CONNECT",
    "ALIASES",
]

# ---------------------------------------------------------------------------
# Factory
# ---------------------------------------------------------------------------


def load_dna_structure(residue_path: str | None = None) -> structure.Structure:
    """
    Build the DNA Structure object (inline replacement for XMLStructure('DNA.xml')).

    Parameters
    ----------
    residue_path : str | None
        If provided, LEaP init_string will include:
            loadoff {residue_path}/{name}.lib
            loadamberparams {residue_path}/{name}.frcmod
        If None, no LEaP init_string is generated (mirrors XML with no <residuePath>).

    Returns
    -------
    Structure.Structure
    """
    return structure.Structure(
        RESIDUE_NAMES,
        RESIDUE_LENGTH,
        rotating_elements=ROTATIONS,
        backbone_elements=BACKBONE,
        connect=CONNECT,
        residue_path=residue_path,  # was None in the XML
        alias=ALIASES,
    )
