"""
PDB cleaner (programmatic edition)

Original author (2017): Dr. Huan Wang <huan.wang@mail.huji.ac.il>
Copyright: The Hebrew University of Jerusalem,
           and The Open University of Israel.

Light modifications (2025): Siddharth
- Import-safe: removed top-level input prompts; everything is callable
  as functions.
- Fixed slice bug for charge field in PDB parsing (cols 79–80).
- Robust parsing: tolerates short lines (pads to 80 cols) and missing
  TER lines. Treats NMR/X-ray uniformly.
- AltLoc resolution: keeps per-residue/chain/atom the entry with the
  highest occupancy; tie-breaks by blank altLoc, then alphabetic.
- DNA/RNA preserved by whitelist; only truly non-standard residues are
  dropped.
- Consistent output naming: "{stem}_cleaned.pdb" regardless of chain
  policy.
- Programmatic API:
    * clean_one_file(input_path, keep="all"/"one"/"A"/"A,B",
      remove_h=False, drop_hetatm=False)
    * batch_clean(dir_path, **kwargs)
- Optional removal of HETATM (ligands/solvents/metals) and hydrogens.
- Strict writer: exact PDB column widths, proper chain/res split,
  %8.3f XYZ, trailing TER and END.
"""

from __future__ import annotations

import os
import re
from collections.abc import Iterable, Sequence
from datetime import datetime

import numpy as np
import pandas as pd
from numpy import char as _npchar

# --------------------------- Constants & Helpers ---------------------------

PDB_ITEMS: list[str] = [
    "Records",  # "ATOM", "HETATM", "TER"
    "AtomSeq",  # serial (string)
    "AtomTyp",  # atom name (aligned)
    "Alt_Loc",  # altloc flag: ' ', 'A', 'B', ...
    "ResName",  # residue name
    "ChainID",  # chain identifier
    "Seq_Num",  # residue sequence number (string)
    "InsCode",  # insertion code
    "Coord_X",  # x
    "Coord_Y",  # y
    "Coord_Z",  # z
    "SD_Occp",  # occupancy
    "SD_Temp",  # B-factor / temp factor
    "Element",  # element symbol
    "Charges",  # charge
]

AMINO_ACIDS = _npchar.asarray(
    [
        "Ala",
        "Arg",
        "Asn",
        "Asp",
        "Cys",
        "Gln",
        "Glu",
        "Gly",
        "His",
        "Ile",
        "Leu",
        "Lys",
        "Met",
        "Phe",
        "Pro",
        "Ser",
        "Thr",
        "Trp",
        "Tyr",
        "Val",
    ]
).upper()

# Nucleic-acid residue names to preserve. Upper-cased for comparison.
NA_RESIDUES = {
    "DA",
    "DC",
    "DG",
    "DT",
    "A",
    "C",
    "G",
    "U",
    "ADE",
    "CYT",
    "GUA",
    "THY",
    "URA",
    "PSU",
    "H2U",
    "1MA",
    "5MC",
    "5MU",
    "OMG",
    "I",
}


def _pad80(s: str) -> str:
    """Return a string of at least 80 characters (PDB width)."""
    return (s.rstrip("\n") + " " * 80)[:80]


def _extract_method(lines: Sequence[str]) -> str:
    """Extract EXPDTA method if present; else empty string."""
    for line in lines:
        if line.startswith("EXPDTA"):
            m = re.search(r"EXPDTA\s+(.+\w+)", line)
            return m.group(1) if m else ""
    return ""


def _rk_altloc(c: str) -> int:
    """
    Rank altLoc for tie-breaking:
      blank (' ') = best (-1), then 'A' = 0, 'B' = 1, ...
    """
    return -1 if (c == " " or c == "") else (ord(c) - ord("A"))


# --------------------------- Parsing ---------------------------


def pdb_structure(line: str) -> list[str]:
    """
    Parse one PDB line (ATOM/HETATM/TER) into fixed fields.

    Defensive against short lines by padding to 80 columns.
    """
    s = _pad80(line)

    return [
        s[0:6].strip(),  # 0 Records
        s[6:12].strip(),  # 1 AtomSeq
        s[12:16].strip(),  # 2 AtomTyp
        s[16],  # 3 Alt_Loc
        s[17:20].strip(),  # 4 ResName
        s[21],  # 5 ChainID
        s[22:26].strip(),  # 6 Seq_Num
        s[26],  # 7 InsCode
        s[30:38].strip(),  # 8 Coord_X
        s[38:46].strip(),  # 9 Coord_Y
        s[46:54].strip(),  # 10 Coord_Z
        s[54:60].strip(),  # 11 SD_Occp
        s[60:66].strip(),  # 12 SD_Temp
        s[76:78].strip(),  # 13 Element
        s[78:80].strip(),  # 14 Charges
    ]


def _parse_pdb_lines(lines: Sequence[str]) -> pd.DataFrame:
    """Convert PDB file lines to a DataFrame with columns PDB_ITEMS."""
    data: list[list[str]] = []
    for line in lines:
        if line.startswith(("ATOM", "HETATM", "TER")):
            data.append(pdb_structure(line))
    df = pd.DataFrame(data, columns=PDB_ITEMS)
    return df


def pdb_reader(filename: str) -> tuple[str, pd.DataFrame, pd.DataFrame]:
    """
    Read a PDB file and split into (method, main_pdb_df, ligand_df).

    - 'method' from EXPDTA if present; otherwise "".
    - 'main_pdb_df' includes records up to the last TER (if any), else all.
    - 'ligand_df' includes records after the last TER (if any), else empty.
    """
    with open(filename) as f:
        lines = f.readlines()

    method = _extract_method(lines)
    pdb_df = _parse_pdb_lines(lines)

    ter_idx = pdb_df.index[pdb_df["Records"] == "TER"].tolist()
    if ter_idx:
        last_ter = ter_idx[-1]
        pdb_main = pdb_df.iloc[: last_ter + 1, :].copy()
        ligand = pdb_df.iloc[last_ter + 1 :, :].copy()
    else:
        pdb_main = pdb_df
        ligand = pd.DataFrame(columns=PDB_ITEMS)

    return method, pdb_main, ligand


# --------------------------- Checks & Reports ---------------------------


def find_pdb_files(path: str) -> Iterable[str]:
    """Yield *.pdb files in a directory (non-recursive)."""
    for f in os.listdir(path):
        if f.lower().endswith(".pdb"):
            yield f


def check_ligand(filename: str, ligand: pd.DataFrame) -> tuple[str, np.ndarray] | None:
    """Return (filename, unique ligand ResName) if ligands exist; else None."""
    if ligand is not None and not ligand.empty:
        return (filename, ligand.ResName.unique())
    return None


def check_altloc(filename: str, pdb_df: pd.DataFrame) -> tuple[str, np.ndarray] | None:
    """Return (filename, altLoc set) if any altLoc != ' '; else None."""
    altloc = pdb_df.Alt_Loc.fillna(" ").unique()
    if len(altloc) > 1:
        return (filename, np.sort(altloc))
    return None


def non_std_residues(
    filename: str, pdb_df: pd.DataFrame
) -> tuple[str, np.ndarray] | None:
    """
    Identify non-standard residues (excluding amino acids and nucleic acids).
    Returns (filename, unique_resnames) or None.
    """
    allowed = set(AMINO_ACIDS.tolist()) | NA_RESIDUES
    upper = pdb_df.ResName.astype(str).str.upper()
    mask = ~upper.isin(allowed)
    vals = upper[mask].unique()
    if vals.size:
        return (filename, vals)
    return None


def check_negative_seqnum(filename: str, pdb_df: pd.DataFrame) -> str | None:
    """Report filename if any residue has a negative sequence number."""
    seq_num = pd.to_numeric(pdb_df.Seq_Num, errors="coerce")
    if (seq_num < 0).any():
        return filename
    return None


def check_sequence_gaps(
    filename: str, pdb_df: pd.DataFrame
) -> tuple[str, list[tuple[str, str]]] | None:
    """
    Find sequence gaps (|diff| > 1) over overall order (not per chain).
    Returns (filename, [(A:223, A:237), ...]) or None.
    """
    df = pdb_df[pdb_df["Records"].isin(["ATOM", "HETATM"])]
    if df.empty:
        return None

    seq = pd.to_numeric(df["Seq_Num"], errors="coerce").values
    chain = df["ChainID"].astype(str).values
    seq_diff = np.abs(np.diff(seq))

    if np.any(seq_diff > 1):
        idx = np.where(seq_diff > 1)[0]
        heads = [f"{chain[i]}:{int(seq[i])}" for i in idx]
        tails = [f"{chain[i + 1]}:{int(seq[i + 1])}" for i in idx]
        return (filename, list(zip(heads, tails, strict=False)))
    return None


def check_insertion_code(
    filename: str, pdb_df: pd.DataFrame
) -> tuple[str, np.ndarray] | None:
    """Report insertion codes present (if any besides blank)."""
    insert = pdb_df.InsCode.fillna(" ").unique()
    if len(insert) > 1:
        return (filename, np.sort(insert))
    return None


def check_multiple_chains(
    filename: str, pdb_df: pd.DataFrame
) -> tuple[str, np.ndarray] | None:
    """Report multiple chains if present."""
    chains = pdb_df.ChainID.fillna("").unique()
    if len(chains) > 1:
        return (filename, chains)
    return None


def check_hydrogen(filename: str, pdb_df: pd.DataFrame) -> str | None:
    """
    Report the filename if hydrogen atoms exist.

    Detection: any row with Element == 'H' or AtomTyp starting with 'H'.
    """
    has_h = (pdb_df.Element == "H").any() or (
        pdb_df.AtomTyp.astype(str).str.startswith("H").any()
    )
    return filename if has_h else None


def save_report(
    path: str,
    ligand_info: list,
    altloc_info: list,
    non_std_Res: list,
    hydrogens: list,
    seqGap_info: list,
    insert_info: list,
    multiChains: list,
    negativeSeq: list,
    drawline: str,
    out_name: str = "special_PDB_files.txt",
) -> None:
    """Write a summary report of special conditions found in PDBs."""
    report = out_name
    string = "The files below have"

    with open(os.path.join(path, report), "w") as fw:
        now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        fw.write(f"Summary\t({now})\n")

        def _block(title: str, header: str, rows: Iterable[tuple[str, Iterable[str]]]):
            fw.write(f"{drawline}{string} {title}\n")
            for item in rows:
                if not item:
                    continue
                name, values = item
                values = list(values)
                fw.write(f"{name}\t{header} ")
                fw.write("".join(f"{v:<6s}" for v in values))
                fw.write("\n")

        _block("ligands", "ligands:", ligand_info)
        _block("alternate locations", "alternate locations:", altloc_info)
        _block("non-standard residues", "non-standard residues:", non_std_Res)
        fw.write(f"{drawline}{string} hydrogen atoms\n")
        for h in hydrogens:
            if h:
                fw.write(f"{h}\n")
        fw.write(f"{drawline}{string} sequence gaps\n")
        for g in seqGap_info:
            if g:
                name, gaps = g
                fw.write(f"{name}\tsequence gaps: ")
                fw.write(" ".join(str(x) for x in gaps))
                fw.write("\n")
        _block("insertion code", "insertion code:", insert_info)
        _block("multiple chains", "multiple chains:", multiChains)
        fw.write(f"{drawline}{string} negative sequence number\n")
        for s in negativeSeq:
            if s:
                fw.write(f"{s}\n")


# --------------------------- Strict Writer ---------------------------


def _fmt_float(v: str | float | int, width: int, prec: int) -> str:
    """Right-align float to PDB width; blanks if missing/unparseable."""
    try:
        return f"{float(v):>{width}.{prec}f}"
    except Exception:
        return " " * width


def _atom_field(atom: str, element: str) -> str:
    """
    PDB atom-name alignment:
      - If atom name is 4 chars, no leading space.
      - Else left pad one space so it sits in cols 13–16 correctly.
    """
    a = (atom or "").strip()
    if len(a) == 4:
        return f"{a:>4s}"
    return f" {a:>3s}"


def _format_ter_line(row: dict) -> str:
    """
    TER line fields:
      1-6 'TER', 7-11 serial, 18-20 resName, 22 chain, 23-26 resSeq, 27 iCode.
    """
    serial = str(row.get("AtomSeq") or "").strip()
    res = (row.get("ResName") or "")[:3]
    chain = (row.get("ChainID") or " ")[:1]
    seq = str(row.get("Seq_Num") or "")
    ins = (row.get("InsCode") or " ")[:1]
    return (
        f"{'TER':<6s}"
        f"{serial:>5s}"
        f"{'':6s}"  # 12–17
        f"{res:>3s} "
        f"{chain:1s}"
        f"{seq:>4s}"
        f"{ins:1s}\n"
    )


def save_cleaned_pdb(outputf: str, pdb_df: pd.DataFrame) -> None:
    """
    Write a cleaned DataFrame back to strict PDB fixed-width format.

    Columns (PDB 3.3 style):
      1-6  Record name
      7-11 Atom serial number
      13-16 Atom name (special alignment; see _atom_field)
      17    altLoc
      18-20 resName
      21    blank
      22    chainID
      23-26 resSeq
      27    iCode
      28-30 blanks
      31-38 X (8.3)
      39-46 Y (8.3)
      47-54 Z (8.3)
      55-60 occupancy (6.2)
      61-66 tempFactor (6.2)
      67-76 blanks
      77-78 element
      79-80 charge
    """
    last_res_row = None
    had_ter = False

    with open(outputf, "w") as fw:
        for tup in pdb_df.itertuples(index=False):
            row = {
                "Records": tup[0],
                "AtomSeq": tup[1],
                "AtomTyp": tup[2],
                "Alt_Loc": tup[3],
                "ResName": tup[4],
                "ChainID": tup[5],
                "Seq_Num": tup[6],
                "InsCode": tup[7],
                "Coord_X": tup[8],
                "Coord_Y": tup[9],
                "Coord_Z": tup[10],
                "SD_Occp": tup[11],
                "SD_Temp": tup[12],
                "Element": tup[13],
                "Charges": tup[14],
            }

            rec = (row["Records"] or "").strip().upper()

            if rec == "TER":
                fw.write(_format_ter_line(row))
                had_ter = True
                continue

            if rec not in ("ATOM", "HETATM"):
                # Ignore any other record types (HEADER, REMARK, ...)
                continue

            last_res_row = row

            rec_f = f"{rec:<6s}"
            serial_f = f"{str(row['AtomSeq'] or ''):>5s}"
            atom_f = _atom_field(str(row["AtomTyp"] or ""), str(row["Element"] or ""))
            alt_f = f"{(row['Alt_Loc'] or ' ')[:1]:1s}"
            res_f = f"{(row['ResName'] or '')[:3]:>3s}"
            chain_f = f"{(row['ChainID'] or ' ')[:1]:1s}"
            seq_f = f"{str(row['Seq_Num'] or ''):>4s}"
            ins_f = f"{(row['InsCode'] or ' ')[:1]:1s}"
            x_f = _fmt_float(row["Coord_X"], 8, 3)
            y_f = _fmt_float(row["Coord_Y"], 8, 3)
            z_f = _fmt_float(row["Coord_Z"], 8, 3)
            occ_f = _fmt_float(row["SD_Occp"], 6, 2)
            b_f = _fmt_float(row["SD_Temp"], 6, 2)
            spacer_67_76 = " " * 10
            elem_f = f"{(row['Element'] or '')[:2]:>2s}"
            chg_f = f"{(row['Charges'] or '')[:2]:>2s}"

            line = (
                rec_f
                + serial_f
                + " "
                + atom_f
                + alt_f
                + res_f
                + " "
                + chain_f
                + seq_f
                + ins_f
                + "   "
                + x_f
                + y_f
                + z_f
                + occ_f
                + b_f
                + spacer_67_76
                + elem_f
                + chg_f
                + "\n"
            )
            fw.write(line)

        if not had_ter and last_res_row is not None:
            fw.write(_format_ter_line(last_res_row))

        fw.write("END\n")


# --------------------------- Cleaning Pipeline ---------------------------


def _select_chains(pdb_df: pd.DataFrame, keep: str) -> pd.DataFrame:
    """
    Apply chain policy:
      - "all": keep all chains (default)
      - "one": keep the longest chain (mode of ChainID)
      - "A" or "A,B": keep listed chains
    """
    keep = (keep or "all").strip().upper()
    if keep == "ALL":
        return pdb_df

    if keep == "ONE":
        counts = pdb_df[
            pdb_df["Records"].isin(["ATOM", "HETATM"])
        ].ChainID.value_counts()
        if counts.empty:
            return pdb_df
        top_chain = counts.index[0]
        return pdb_df[pdb_df.ChainID == top_chain]

    wanted = {c.strip() for c in keep.split(",") if c.strip()}
    if wanted:
        return pdb_df[pdb_df.ChainID.str.upper().isin(wanted)]
    return pdb_df


def _resolve_altloc(pdb_df: pd.DataFrame) -> pd.DataFrame:
    """
    Collapse alternative locations to a single record per
    (Seq_Num, ChainID, AtomTyp).
    """
    if "Alt_Loc" not in pdb_df.columns or pdb_df.Alt_Loc.fillna(" ").nunique() <= 1:
        return pdb_df

    occ = pd.to_numeric(pdb_df["SD_Occp"], errors="coerce").fillna(0.0)
    alt = pdb_df["Alt_Loc"].astype(str).fillna(" ")
    rk = alt.map(_rk_altloc)

    df = pdb_df.assign(__occ=occ, __rk=rk)
    df = (
        df.sort_values(
            ["Seq_Num", "ChainID", "AtomTyp", "__occ", "__rk"],
            ascending=[True, True, True, False, True],
        )
        .drop_duplicates(subset=["Seq_Num", "ChainID", "AtomTyp"], keep="first")
        .drop(columns=["__occ", "__rk"])
    )
    return df


def clean_one_file(
    input_path: str,
    keep: str = "all",
    remove_h: bool = False,
    drop_hetatm: bool = False,
) -> str:
    """
    Programmatic single-file cleaner for MAWS (call before running tleap).

    Args:
        input_path: Path to a PDB file.
        keep: "all" | "one" | list of chains like "A" or "A,B".
        remove_h: If True, drop hydrogen atoms (Element == 'H').
        drop_hetatm: If True, drop all HETATM records (ligands/solvents/metals).

    Returns:
        Path to the cleaned PDB: "{stem}_cleaned.pdb" alongside the input.
    """
    input_path = os.path.abspath(input_path)
    path, fname = os.path.dirname(input_path), os.path.basename(input_path)

    _method, pdb_df, _ligand = pdb_reader(input_path)

    # Chain policy
    pdb_df = _select_chains(pdb_df, keep)

    # Optional removals
    if drop_hetatm:
        pdb_df = pdb_df[pdb_df["Records"] != "HETATM"]
    if remove_h:
        pdb_df = pdb_df[pdb_df["Element"].astype(str).str.upper() != "H"]

    # Remove insertion codes (keep only blank)
    if pdb_df.InsCode.fillna(" ").nunique() > 1:
        pdb_df = pdb_df[pdb_df.InsCode.fillna(" ") == " "]

    # Resolve altLocs
    pdb_df = _resolve_altloc(pdb_df)

    # Remove non-standard residues (keep amino acids and nucleic acids)
    ns = non_std_residues(fname, pdb_df)
    if ns:
        bad = set(ns[1].tolist())
        mask = ~pdb_df.ResName.astype(str).str.upper().isin(bad)
        pdb_df = pdb_df[mask]

    # Output
    stem, _ = os.path.splitext(fname)
    out_path = os.path.join(path, f"{stem}_cleaned.pdb")
    save_cleaned_pdb(out_path, pdb_df)
    return out_path


# def batch_clean(
#     dir_path: str,
#     keep: str = "all",
#     remove_h: bool = False,
#     drop_hetatm: bool = False,
#     write_report: bool = True,
# ) -> None:
#     """
#     Clean all *.pdb files in a directory (non-recursive).

#     Writes "{stem}_cleaned.pdb" for each input and an optional summary:
#     "special_PDB_files.txt".
#     """
#     ligand_info = []
#     altloc_info = []
#     non_std_Res = []
#     Hatoms_info = []
#     negativeSeq = []
#     seqGap_info = []
#     insert_info = []
#     multiChains = []

#     drawline = "\n" + "-" * 79 + "\n"
#     time_fmt = "Step time: {:.4f}s"

#     start_all = time.time()
#     for i, f in enumerate(find_pdb_files(dir_path), start=1):
#         t0 = time.time()
#         print(f"{drawline}Check Point:{i:>6}\tPDB: {f}")

#         filename = os.path.join(dir_path, f)
#         _method, pdb_df, ligand = pdb_reader(filename)

#         # Collect info
#         lig = check_ligand(filename, ligand)
#         if lig:
#             ligand_info.append(lig)
#         alt = check_altloc(filename, pdb_df)
#         if alt:
#             altloc_info.append(alt)
#         ns = non_std_residues(filename, pdb_df)
#         if ns:
#             non_std_Res.append(ns)
#         neg = check_negative_seqnum(filename, pdb_df)
#         if neg:
#             negativeSeq.append(neg)
#         gaps = check_sequence_gaps(filename, pdb_df)
#         if gaps:
#             seqGap_info.append(gaps)
#         ins = check_insertion_code(filename, pdb_df)
#         if ins:
#             insert_info.append(ins)
#         chs = check_multiple_chains(filename, pdb_df)
#         if chs:
#             multiChains.append(chs)
#         h = check_hydrogen(filename, pdb_df)
#         if h:
#             Hatoms_info.append(h)

#         # Clean & save
#         out = clean_one_file(
#             filename, keep=keep, remove_h=remove_h, drop_hetatm=drop_hetatm
#         )
#         print(f"Saved cleaned PDB: {out}")
#         print(time_fmt.format(time.time() - t0))

#     if write_report:
#         save_report(
#             dir_path,
#             ligand_info,
#             altloc_info,
#             non_std_Res,
#             Hatoms_info,
#             seqGap_info,
#             insert_info,
#             multiChains,
#             negativeSeq,
#             drawline,
#         )
#         print("Wrote summary report: special_PDB_files.txt")

#     print(f"{drawline}All done in {time.time() - start_all:.3f}s")


# # --------------------------- Optional CLI ---------------------------


# def _build_cli() -> argparse.ArgumentParser:
#     p = argparse.ArgumentParser(
#         description="Clean PDB files programmatically (keeps original attribution)."
#     )
#     sub = p.add_subparsers(dest="cmd", required=True)

#     s1 = sub.add_parser("one", help="Clean a single PDB file")
#     s1.add_argument("input_path", type=str)
#     s1.add_argument(
#         "--keep",
#         type=str,
#         default="all",
#         help='Chain policy: "all", "one", or list "A" / "A,B"',
#     )
#     s1.add_argument("--remove-h", action="store_true", help="Remove hydrogens")
#     s1.add_argument(
#         "--drop-hetatm",
#         action="store_true",
#         help="Remove HETATM records (ligands/solvents)",
#     )

#     s2 = sub.add_parser(
#         "dir", help="Clean all PDB files in a directory (non-recursive)"
#     )
#     s2.add_argument("dir_path", type=str)
#     s2.add_argument("--keep", type=str, default="all")
#     s2.add_argument("--remove-h", action="store_true")
#     s2.add_argument("--drop-hetatm", action="store_true")
#     s2.add_argument(
#         "--no-report", action="store_true", help="Do not write the summary report"
#     )

#     return p


# if __name__ == "__main__":
#     cli = _build_cli()
#     args = cli.parse_args()

#     if args.cmd == "one":
#         out = clean_one_file(
#             args.input_path,
#             keep=args.keep,
#             remove_h=args.remove_h,
#             drop_hetatm=args.drop_hetatm,
#         )
#         print(out)
#     elif args.cmd == "dir":
#         batch_clean(
#             args.dir_path,
#             keep=args.keep,
#             remove_h=args.remove_h,
#             drop_hetatm=args.drop_hetatm,
#             write_report=not args.no_report,
#         )
