"""
maws.forcefield
===============

Which physics to use for a design run, as one value.

MAWS designs an *aptamer*: a short strand of RNA or DNA that folds up against a
target molecule and sticks to it. To score how well a given shape binds, both
molecules need a *force field* — the collection of numbers saying how strongly
each pair of atoms pulls on or pushes away the other. The aptamer and the
target are different kinds of matter and get different force fields, named
strings such as ``"leaprc.RNA.OL3"`` and ``"leaprc.protein.ff19SB"``. Those
names are what ``tleap``, the AmberTools program that builds structures, is
given to load the parameters.

A run also needs a salt concentration. MAWS does not place individual water
molecules; it uses *implicit solvent*, in which the surrounding water is a
smooth continuum. Salt dissolved in that water damps the attraction and
repulsion between charged atoms, an effect called *screening*. A nucleic acid
backbone carries a negative charge on every residue, so with no screening at
all those charges dominate the score and the strand pushes itself apart.

Both force field names and the salt concentration travel together everywhere,
so :class:`ForceField` holds them as one frozen value rather than three
arguments that could be passed inconsistently.

Examples
--------
>>> ff = ForceField.for_target("RNA", "protein")
>>> ff.aptamer_source
'leaprc.RNA.OL3'
>>> ff.alphabet
'GAUC'
>>> ff.parameterized  # the protein force field already covers proteins
True
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import ClassVar, Literal

from maws.errors import ConfigurationError

__all__ = ["AptamerType", "ForceField", "MoleculeType"]

AptamerType = Literal["RNA", "DNA"]
"""Which nucleic acid the aptamer is made of."""

MoleculeType = Literal["protein", "organic", "lipid"]
"""What kind of molecule the aptamer is designed against."""


@dataclass(frozen=True, slots=True)
class ForceField:
    """The complete physics setup for one design run.

    Nothing here can be changed after the object is made. One run therefore
    cannot build a structure under one force field and then score it under
    another: the same value is used at every step, or a different value is
    made and passed explicitly.

    Parameters
    ----------
    aptamer_source : str
        Name ``tleap`` is given to load the nucleic acid parameters, e.g.
        ``"leaprc.RNA.OL3"``. Choosing an RNA or a DNA name here also settles
        which four nucleotides the search may build with, reported by
        :attr:`alphabet`.
    ligand_source : str
        Name ``tleap`` is given to load the target's parameters, e.g.
        ``"leaprc.protein.ff19SB"``. It has to suit what the target actually
        is; a protein scored with the parameters meant for small organic
        molecules gives numbers that mean nothing.
    salt_conc : float, default=0.15
        Concentration of dissolved salt in mol/L, used to damp the pull between
        charged atoms. The default is roughly the saltiness of blood, which is
        the condition an aptamer is normally wanted to work in. ``0.0`` removes
        the damping entirely, so the aptamer's own backbone charges push it
        open. Only singly charged ions such as Na⁺ are modelled this way, so
        raising this number cannot stand in for Mg²⁺.
    parameterized : bool, default=True
        Whether `ligand_source` already describes the target. When True the
        target's PDB file goes straight to ``tleap``. When False, MAWS first
        runs ``antechamber`` and ``parmchk2`` to work out the missing
        parameters for that particular molecule, which takes noticeably longer.

    See Also
    --------
    for_target : Picks both names from what the aptamer and target are.

    Examples
    --------
    >>> ff = ForceField.for_target("DNA", "organic", salt_conc=0.0)
    >>> ff.ligand_source, ff.parameterized
    ('leaprc.gaff2', False)
    """

    aptamer_source: str
    ligand_source: str
    salt_conc: float = 0.15
    parameterized: bool = True

    _APTAMER: ClassVar[dict[str, str]] = {
        "RNA": "leaprc.RNA.OL3",
        "DNA": "leaprc.DNA.OL21",
    }
    _LIGAND: ClassVar[dict[str, tuple[str, bool]]] = {
        # name -> (parameter set to load, does it already cover the target?)
        "protein": ("leaprc.protein.ff19SB", True),
        "organic": ("leaprc.gaff2", False),
        "lipid": ("leaprc.lipid21", False),
    }

    def __post_init__(self) -> None:
        if self.salt_conc < 0:
            raise ConfigurationError(
                f"salt_conc is a concentration and cannot be negative, got "
                f"{self.salt_conc}"
            )

    @classmethod
    def for_target(
        cls,
        aptamer: AptamerType,
        molecule: MoleculeType,
        *,
        salt_conc: float = 0.15,
    ) -> ForceField:
        """Return the usual force field pairing for a kind of design run.

        Saves having to remember the parameter set names: say what the two
        molecules are and this picks a matching pair.

        Parameters
        ----------
        aptamer : {"RNA", "DNA"}
            Which nucleic acid the aptamer is made of.
        molecule : {"protein", "organic", "lipid"}
            What the target is. ``"protein"`` is the only one the parameters
            already cover, so the other two mean ``antechamber`` and
            ``parmchk2`` have to run before the structure can be built.
        salt_conc : float, default=0.15
            Concentration of dissolved salt in mol/L. See :class:`ForceField`.

        Returns
        -------
        ForceField
            The matching setup.

        Raises
        ------
        maws.errors.ConfigurationError
            If either argument is not one of the listed options. The message
            lists the ones that are.

        Examples
        --------
        >>> ForceField.for_target("RNA", "organic").ligand_source
        'leaprc.gaff2'
        >>> ForceField.for_target("DNA", "protein").aptamer_source
        'leaprc.DNA.OL21'
        """
        if aptamer not in cls._APTAMER:
            raise ConfigurationError(
                f"aptamer must be one of {sorted(cls._APTAMER)}, got {aptamer!r}"
            )
        if molecule not in cls._LIGAND:
            raise ConfigurationError(
                f"molecule must be one of {sorted(cls._LIGAND)}, got {molecule!r}"
            )
        ligand_source, parameterized = cls._LIGAND[molecule]
        return cls(
            aptamer_source=cls._APTAMER[aptamer],
            ligand_source=ligand_source,
            salt_conc=salt_conc,
            parameterized=parameterized,
        )

    @property
    def alphabet(self) -> str:
        """str : The four nucleotides the search may grow the aptamer with.

        RNA is built from G, A, U and C; DNA has T where RNA has U. Worked out
        from `aptamer_source` each time it is asked for rather than stored, so
        it cannot end up naming nucleotides the loaded parameters do not
        describe.

        Examples
        --------
        >>> ForceField.for_target("DNA", "protein").alphabet
        'GATC'
        """
        return "GAUC" if "RNA" in self.aptamer_source else "GATC"

    def leap_preamble(self) -> list[str]:
        """Return the two lines that open every ``tleap`` input script.

        A ``tleap`` script starts by loading its parameters, one ``source``
        line per force field, before it mentions any molecule.

        Returns
        -------
        list of str
            One ``source`` line per force field, the aptamer's first.

        Examples
        --------
        >>> ForceField.for_target("RNA", "protein").leap_preamble()
        ['source leaprc.RNA.OL3', 'source leaprc.protein.ff19SB']
        """
        return [
            f"source {self.aptamer_source}",
            f"source {self.ligand_source}",
        ]
