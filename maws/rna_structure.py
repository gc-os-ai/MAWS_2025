# RNA.py
# Inline RNA residue templates for MAWS (formerly RNA.xml)
# Produces a Structure.Structure identical to Loadfrom.XMLStructure("RNA.xml"), but without I/O.


import maws.structure as structure

# ---------------------------------------------------------------------------
# Core tables (mirrors the order & values in RNA.xml)
# ---------------------------------------------------------------------------

# 1) Residue names (must align with all per-residue lists below)
RESIDUE_NAMES: list[str] = [
    "GN",
    "AN",
    "UN",
    "CN",
    "G",
    "A",
    "U",
    "C",
    "G5",
    "A5",
    "U5",
    "C5",
    "G3",
    "A3",
    "U3",
    "C3",
]

# 2) Atom counts per residue (length=… in XML), aligned with RESIDUE_NAMES
RESIDUE_LENGTH: list[int] = [
    33,
    32,
    29,
    30,
    34,
    33,
    30,
    31,
    32,
    31,
    28,
    29,
    35,
    34,
    31,
    32,
]

# 3) Backbone tuples: (residue, start, middle_pre, bond, middle_post, end)
BACKBONE: list[tuple[str, int, int, int, int, int]] = [
    ("GN", 0, 8, 10, 25, 32),
    ("AN", 0, 8, 10, 24, 31),
    ("UN", 0, 8, 10, 21, 28),
    ("CN", 0, 8, 10, 22, 29),
    ("G", 0, 10, 12, 27, 33),
    ("A", 0, 10, 12, 26, 32),
    ("U", 0, 10, 12, 23, 29),
    ("C", 0, 10, 12, 24, 30),
    ("G5", 0, 8, 10, 25, 31),
    ("A5", 0, 8, 10, 24, 30),
    ("U5", 0, 8, 10, 21, 27),
    ("C5", 0, 8, 10, 22, 28),
    ("G3", 0, 10, 12, 27, 34),
    ("A3", 0, 10, 12, 26, 33),
    ("U3", 0, 10, 12, 23, 30),
    ("C3", 0, 10, 12, 24, 31),
]

# 4) Rotations: (residue, start, bond, end_or_None)
# Note: 'end' in XML becomes None here; negative indices are kept (normalized at runtime)
ROTATIONS: list[tuple[str, int, int, int | None]] = [
    # GN
    ("GN", 0, 1, None),
    ("GN", 1, 2, None),
    ("GN", 8, 10, -8),
    ("GN", -8, -2, None),
    # AN
    ("AN", 0, 1, None),
    ("AN", 1, 2, None),
    ("AN", 8, 10, -8),
    ("AN", -8, -2, None),
    # UN
    ("UN", 0, 1, None),
    ("UN", 1, 2, None),
    ("UN", 8, 10, -8),
    ("UN", -8, -2, None),
    # CN
    ("CN", 0, 1, None),
    ("CN", 1, 2, None),
    ("CN", 8, 10, -8),
    ("CN", -8, -2, None),
    # G
    ("G", 0, 3, None),
    ("G", 3, 4, None),
    ("G", 10, 12, -7),
    ("G", -7, -1, None),
    # A
    ("A", 0, 3, None),
    ("A", 3, 4, None),
    ("A", 10, 12, -7),
    ("A", -7, -1, None),
    # U
    ("U", 0, 3, None),
    ("U", 3, 4, None),
    ("U", 10, 12, -7),
    ("U", -7, -1, None),
    # C
    ("C", 0, 3, None),
    ("C", 3, 4, None),
    ("C", 10, 12, -7),
    ("C", -7, -1, None),
    # G5
    ("G5", 0, 1, None),
    ("G5", 1, 2, None),
    ("G5", 8, 10, -7),
    ("G5", -7, -1, None),
    # A5
    ("A5", 0, 1, None),
    ("A5", 1, 2, None),
    ("A5", 8, 10, -7),
    ("A5", -7, -1, None),
    # U5
    ("U5", 0, 1, None),
    ("U5", 1, 2, None),
    ("U5", 8, 10, -7),
    ("U5", -7, -1, None),
    # C5
    ("C5", 0, 1, None),
    ("C5", 1, 2, None),
    ("C5", 8, 10, -7),
    ("C5", -7, -1, None),
    # G3
    ("G3", 0, 3, None),
    ("G3", 3, 4, None),
    ("G3", 10, 12, -8),
    ("G3", -8, -2, None),
    # A3
    ("A3", 0, 3, None),
    ("A3", 3, 4, None),
    ("A3", 10, 12, -8),
    ("A3", -8, -2, None),
    # U3
    ("U3", 0, 3, None),
    ("U3", 3, 4, None),
    ("U3", 10, 12, -8),
    ("U3", -8, -2, None),
    # C3
    ("C3", 0, 3, None),
    ("C3", 3, 4, None),
    ("C3", 10, 12, -8),
    ("C3", -8, -2, None),
]

# 5) Connectivity: per residue: [[append_first, append_last], [prepend_last, prepend_first], append_len, prepend_len]
CONNECT: list[list] = [[[0, -1], [-2, 0], 1.6, 1.6] for _ in RESIDUE_NAMES]

# 6) Alias table: [name, alone, start, middle, end]
ALIASES: list[list[str]] = [
    ["CN", "CN", "C5", "C", "C3"],
    ["A", "AN", "A5", "A", "A3"],
    ["C", "CN", "C5", "C", "C3"],
    ["UN", "UN", "U5", "U", "U3"],
    ["G", "GN", "G5", "G", "G3"],
    ["G3", "G3", "G", "G", "G3"],
    ["AN", "AN", "A5", "A", "A3"],
    ["A3", "A3", "A", "A", "A3"],
    ["GN", "GN", "G5", "G", "G3"],
    ["A5", "A5", "A5", "A", "A"],
    ["U", "UN", "U5", "U", "U3"],
    ["G5", "G5", "G5", "G", "G"],
    ["U3", "U3", "U", "U", "U3"],
    ["C5", "C5", "C5", "C", "C"],
    ["C3", "C3", "C", "C", "C3"],
    ["U5", "U5", "U5", "U", "U"],
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


def load_rna_structure(residue_path: str | None = None) -> structure.Structure:
    """
    Build the RNA Structure object (inline replacement for XMLStructure('RNA.xml')).

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
        residue_path=residue_path,  # the XML had no <residuePath>
        alias=ALIASES,
    )
