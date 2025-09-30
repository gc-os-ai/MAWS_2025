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

# 4) Rotations: (residue, start, bond, end_or_None)
# Note: 'end' in XML becomes None here; negative indices are kept (normalized at runtime)
ROTATIONS: list[tuple[str, int, int, int | None]] = [
    # DGN
    ("DGN", 0, 1, None),
    ("DGN", 1, 2, None),
    ("DGN", 8, 10, -7),
    ("DGN", -4, -2, None),
    # DAN
    ("DAN", 0, 1, None),
    ("DAN", 1, 2, None),
    ("DAN", 8, 10, -7),
    ("DAN", -4, -2, None),
    # DTN
    ("DTN", 0, 1, None),
    ("DTN", 1, 2, None),
    ("DTN", 8, 10, -7),
    ("DTN", -4, -2, None),
    # DCN
    ("DCN", 0, 1, None),
    ("DCN", 1, 2, None),
    ("DCN", 8, 10, -7),
    ("DCN", -4, -2, None),
    # DG
    ("DG", 0, 3, None),
    ("DG", 3, 4, None),
    ("DG", 10, 12, -6),
    ("DG", -6, -1, None),
    # DA
    ("DA", 0, 3, None),
    ("DA", 3, 4, None),
    ("DA", 10, 12, -6),
    ("DA", -6, -1, None),
    # DT
    ("DT", 0, 3, None),
    ("DT", 3, 4, None),
    ("DT", 10, 12, -6),
    ("DT", -6, -1, None),
    # DC
    ("DC", 0, 3, None),
    ("DC", 3, 4, None),
    ("DC", 10, 12, -6),
    ("DC", -6, -1, None),
    # DG5
    ("DG5", 0, 1, None),
    ("DG5", 1, 2, None),
    ("DG5", 8, 10, -6),
    ("DG5", -6, -1, None),
    # DA5
    ("DA5", 0, 1, None),
    ("DA5", 1, 2, None),
    ("DA5", 8, 10, -6),
    ("DA5", -6, -1, None),
    # DT5
    ("DT5", 0, 1, None),
    ("DT5", 1, 2, None),
    ("DT5", 8, 10, -6),
    ("DT5", -6, -1, None),
    # DC5
    ("DC5", 0, 1, None),
    ("DC5", 1, 2, None),
    ("DC5", 8, 10, -6),
    ("DC5", -6, -1, None),
    # DG3
    ("DG3", 0, 3, None),
    ("DG3", 3, 4, None),
    ("DG3", 10, 12, -7),
    ("DG3", -4, -2, None),
    # DA3
    ("DA3", 0, 3, None),
    ("DA3", 3, 4, None),
    ("DA3", 10, 12, -7),
    ("DA3", -4, -2, None),
    # DT3
    ("DT3", 0, 3, None),
    ("DT3", 3, 4, None),
    ("DT3", 10, 12, -7),
    ("DT3", -4, -2, None),
    # DC3
    ("DC3", 0, 3, None),
    ("DC3", 3, 4, None),
    ("DC3", 10, 12, -7),
    ("DC3", -4, -2, None),
]

# 5) Connectivity: per residue: [[append_first, append_last], [prepend_last, prepend_first], append_len, prepend_len]
# In this XML, all residues share the same connectivity (0,-1) and (-2,0) with bond lengths 1.6.
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
