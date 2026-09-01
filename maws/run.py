from __future__ import annotations

import copy
import logging
from collections.abc import Sequence
from dataclasses import dataclass
from operator import attrgetter
from pathlib import Path
from typing import Literal, NamedTuple

import numpy as np
from openmm import app

import maws.space as space
from maws.complex import Complex
from maws.dna_structure import load_dna_structure
from maws.pdb_cleaner import resolve_pdb_path
from maws.rna_structure import load_rna_structure
from maws.scoring import entropy_score

AptamerType = Literal["RNA", "DNA"]
MoleculeType = Literal["protein", "organic", "lipid"]
SamplerMode = Literal["sphere", "surface-following"]
PDBInput = str | Path


class Candidate(NamedTuple):
    """One nucleotide choice a search step scored.

    Parameters
    ----------
    entropy : float
        Score from :func:`maws.scoring.entropy_score` for the step that
        added this nucleotide. Lower is better.
    total : float
        Sum of `entropy` over every step so far, including this one. EFBA's
        entropy is extensive, so this is the score of the whole partial
        aptamer, and it is what candidates are ranked on.
    energy : float
        Lowest total potential energy seen while sampling this candidate,
        in kJ/mol.
    sequence : str
        Aptamer sequence this candidate would give.
    positions : list of openmm.Vec3
        Coordinates of the whole complex at the lowest-energy pose.
    topology : openmm.app.Topology, optional
        Topology matching `positions`. Carried by callers that write a PDB
        for every step; left unset when the sequence is enough to rebuild
        it.

    See Also
    --------
    select_beam : Chooses which of these survive a step.
    """

    entropy: float
    total: float
    energy: float
    sequence: str
    positions: object
    topology: object = None


def select_beam(candidates, width):
    """Return the `width` best candidates, lowest running total first.

    Ranks on :attr:`Candidate.total` rather than the score of the step that
    produced each candidate. Candidates in one step can descend from
    different parents, and those parents scored differently. Ranking on the
    step alone would weigh a strong lineage against a weak one as though
    their histories were equal, letting one lucky step displace a
    consistently better sequence.

    With ``width=1`` every candidate shares a parent, so the common part of
    the total cancels and the order matches the step scores.

    Parameters
    ----------
    candidates : iterable of Candidate
        Every nucleotide choice scored in one search step.
    width : int
        How many to carry into the next step. 1 gives a greedy search.

    Returns
    -------
    list of Candidate
        At most `width` candidates, ordered by running total. Shorter than
        `width` when fewer were scored.

    See Also
    --------
    MawsRunner : Sets `width` from its ``beam`` parameter.

    Examples
    --------
    The second candidate scored better on this step, and still loses.

    >>> from maws.run import Candidate, select_beam
    >>> scored = [
    ...     Candidate(-0.2, -2.0, 0.0, "GG", None),
    ...     Candidate(-0.9, -1.0, 0.0, "AA", None),
    ... ]
    >>> [c.sequence for c in select_beam(scored, width=1)]
    ['GG']
    """
    return sorted(candidates, key=attrgetter("total"))[:width]


@dataclass(frozen=True)
class MawsResult:
    """
    Result of a MAWS aptamer design run.

    Attributes
    ----------
    sequence : str
        Best aptamer sequence found.
    energy : float
        Lowest total potential energy of the whole complex over the poses
        sampled for the winning candidate, in kJ/mol.

        .. warning::
            This is not a binding energy and nothing in the search branches
            on it. There is no unbound reference subtracted, so it is
            dominated by the target's own internal energy, and candidates
            with different atom counts do not produce comparable values.
            Selection is decided entirely by `entropy`. See issue #49 (C4).
    entropy : float
        Score the selection was actually made on, from
        :func:`maws.scoring.entropy_score`. At most 0, and the candidate
        scoring lowest is the one kept.
    pdb_path : str
        Path to the saved result PDB file produced by the run.
    seed : int
        The seed the run actually used. Passing it back as
        ``MawsRunner(seed=...)`` reproduces this result.
    """

    sequence: str
    energy: float
    entropy: float
    pdb_path: str | None = None
    seed: int | None = None


class MawsRunner:
    r"""Design an aptamer against a target molecule.

    Grows a strand of DNA or RNA one nucleotide at a time. Each step scores
    every nucleotide that could be added, using the entropic criterion of
    :func:`maws.scoring.entropy_score`, and carries the best `beam`
    candidates into the next step.

    Parameters
    ----------
    num_nucleotides : int
        Length of the aptamer to design.
    aptamer_type : {"RNA", "DNA"}
        Chemistry of the strand being grown.
    molecule_type : {"protein", "organic", "lipid"}
        Chemistry of the target, which selects its force field.
    beam : int, default=1
        How many candidates to carry from one step into the next. 1 commits
        to the single best nucleotide at every step. Larger values keep
        runners-up alive, so a wrong early choice stays recoverable, at a
        cost in run time that grows linearly.

        .. versionadded:: 0.1
    beta : float, default=0.01
        How sharply lower energies are favoured in the score, in mol/kJ.
        A Lagrange multiplier rather than a physical inverse temperature;
        see :func:`maws.scoring.entropy_score`.
    first_chunk_size : int, default=5000
        Conformations sampled per candidate in the first step.
    second_chunk_size : int, default=5000
        Conformations sampled per candidate in every step after the first.
    clean_pdb : bool, default=False
        If True, repairs the input PDB before LEaP reads it. Use for
        protein targets.
    keep_chains : str, default="all"
        Which chains the cleaner keeps: ``"all"``, ``"one"``, or a comma
        separated list such as ``"A,B"``.
    remove_h : bool, default=False
        If True, the cleaner strips hydrogens.
    drop_hetatm : bool, default=False
        If True, the cleaner drops every HETATM record.
    verbose : bool, default=False
        If True, reports each step at INFO level rather than DEBUG.
    sampler_mode : {"surface-following", "sphere"}
        Shape of the region poses are drawn from. ``"surface-following"``
        keeps poses within `d_max` of the target's surface. ``"sphere"``
        fills a ball around the target, most of which is open solvent.
    reach : float, default=10.0
        How far past the target's furthest atom the sampling region
        extends, in angstrom.
    d_max : float, default=6.0
        How far from the target's surface a pose may sit, in angstrom, for
        ``sampler_mode="surface-following"``.
    site_centre : sequence of float, optional
        Sample around this point rather than the whole target, in the input
        PDB's coordinates. Give it when the binding site is known.
    site_radius : float, optional
        How far the region reaches from `site_centre`, in angstrom.
    probe : float, default=1.4
        Radius in angstrom of the ball rolled over the target to find its
        surface. 1.4 is the size of a water molecule.
    clash_tolerance : float, default=1.0
        How far a placed strand may overlap the target's van der Waals
        spheres before the pose is redrawn, in angstrom.
    salt_conc : float, default=0.15
        Monovalent salt concentration in mol/L, for Debye-Huckel screening
        in the implicit solvent. 0 leaves electrostatics unscreened.
    seed : int, optional
        Seed for every random draw, making the run repeatable. Defaults to
        a fresh seed each run, reported in the log and in
        :attr:`MawsResult.seed`.

    See Also
    --------
    maws.scoring.entropy_score : The criterion each step is decided on.
    select_beam : Chooses which candidates survive a step.

    Notes
    -----
    The scoring criterion and the one-nucleotide-at-a-time growth follow the
    Entropic Fragment-Based Approach of Tseng et al. [1]_, the method MAWS
    was built on [2]_.

    ``beam`` is a deliberate departure from that method, added by Siddharth
    in 2026. EFBA specifies a greedy seed-and-grow search, which commits to
    one nucleotide per step and never revisits it, so a wrong choice early
    constrains every step after it. A beam keeps the runners-up alive and
    gives the search a way back. **The default of 1 reproduces the
    published method exactly**; any value above 1 leaves it.

    References
    ----------
    .. [1] Tseng, C.-Y., Ashrafuzzaman, M., Mane, J. Y., Kapty, J., Mercer,
           J. R., Tuszynski, J. A. (2011). "Entropic Fragment-Based Approach
           to Aptamer Design". Chemical Biology & Drug Design 78(1), 1-13.
    .. [2] Kalinowski, M. et al. (2016). "MAWS - Making Aptamers Without
           SELEX". iGEM Heidelberg.

    Examples
    --------
    >>> runner = MawsRunner(  # doctest: +SKIP
    ...     num_nucleotides=15, aptamer_type="RNA", molecule_type="protein"
    ... )
    >>> result = runner.run(pdb="data/1BRQ.pdb")  # doctest: +SKIP
    """

    def __init__(
        self,
        *,
        num_nucleotides: int,
        aptamer_type: AptamerType,
        molecule_type: MoleculeType,
        beam: int = 1,
        beta: float = 0.01,
        first_chunk_size: int = 5000,
        second_chunk_size: int = 5000,
        clean_pdb: bool = False,
        keep_chains: str = "all",
        remove_h: bool = False,
        drop_hetatm: bool = False,
        verbose: bool = False,
        sampler_mode: SamplerMode = "surface-following",
        reach: float = 10.0,
        d_max: float = 6.0,
        site_centre: Sequence[float] | None = None,
        site_radius: float | None = None,
        probe: float = 1.4,
        clash_tolerance: float = 1.0,
        salt_conc: float = 0.15,
        seed: int | None = None,
    ) -> None:
        if num_nucleotides <= 0:
            raise ValueError(
                f"num_nucleotides must be greater than 0, got {num_nucleotides}"
            )
        if beam < 1:
            raise ValueError(f"beam must be >= 1, got {beam}")
        if first_chunk_size <= 0 or second_chunk_size <= 0:
            raise ValueError("Chunk size must be greater than 0")
        if reach < 0:
            raise ValueError(f"reach must be >= 0, got {reach}")
        if probe < 0:
            raise ValueError(f"probe must be >= 0, got {probe}")
        if clash_tolerance < 0:
            raise ValueError(f"clash_tolerance must be >= 0, got {clash_tolerance}")
        if salt_conc < 0:
            raise ValueError(f"salt_conc must be >= 0, got {salt_conc}")
        if seed is not None and not isinstance(seed, int | np.integer):
            raise TypeError(f"seed must be an int or None, got {type(seed).__name__}")

        self.num_nucleotides = num_nucleotides
        self.aptamer_type = aptamer_type
        self.molecule_type = molecule_type
        self.beam = beam
        self.beta = beta
        self.first_chunk_size = first_chunk_size
        self.second_chunk_size = second_chunk_size
        self.clean_pdb = clean_pdb
        self.keep_chains = keep_chains
        self.remove_h = remove_h
        self.drop_hetatm = drop_hetatm
        self.verbose = verbose
        self.sampler_mode = sampler_mode
        self.reach = reach
        self.d_max = d_max
        self.site_centre = site_centre
        self.site_radius = site_radius
        self.probe = probe
        self.clash_tolerance = clash_tolerance
        self.salt_conc = salt_conc
        self.seed = seed

    def run(
        self,
        *,
        pdb: PDBInput,
        name: str = "MAWS_aptamer",
        output_pdb: str | Path | None = None,
    ) -> MawsResult:
        """
        Run the MAWS algorithm.

        Parameters
        ----------
        pdb : str | Path
            Input ligand PDB file path.
        name : str
            Run name used only for log context and artifact naming.
        output_pdb : str | Path | None
            If provided:
              - if it's an existing directory -> writes `{name}_RESULT.pdb` inside it
              - otherwise treated as the exact output file path (parent dirs created)
            If None, no PDB is written.

        Returns
        -------
        MawsResult
        """
        N_BACKBONE_TORSIONS = (
            4  # MAWS rotates 4 backbone torsions per residue in this implementation
        )
        log = logging.getLogger(__name__)

        # A run with no seed still gets one, so its result can be reproduced
        # from the log afterwards.
        seed = np.random.SeedSequence().entropy if self.seed is None else self.seed
        rng = np.random.default_rng(seed)
        log.info("Random seed: %s", seed)

        if self.verbose:
            log.info("MAWS run started: name=%s", name)
        log.debug(
            "Config: num_nucleotides=%d aptamer_type=%s molecule_type=%s beta=%s "
            "c1=%d c2=%d clean_pdb=%s keep_chains=%s remove_h=%s drop_hetatm=%s "
            "salt_conc=%s",
            self.num_nucleotides,
            self.aptamer_type,
            self.molecule_type,
            self.beta,
            self.first_chunk_size,
            self.second_chunk_size,
            self.clean_pdb,
            self.keep_chains,
            self.remove_h,
            self.drop_hetatm,
            self.salt_conc,
        )

        # Resolve (and optionally clean) the PDB path before LEaP calls
        pdb_path, original_pdb_path = resolve_pdb_path(
            str(pdb),
            self.molecule_type,
            clean_pdb=self.clean_pdb,
            keep_chains=self.keep_chains,
            remove_h=self.remove_h,
            drop_hetatm=self.drop_hetatm,
            logger=log,
        )
        log.debug("Input PDB original=%s final=%s", original_pdb_path, pdb_path)

        # Choose aptamer FF and residue template
        if self.aptamer_type == "RNA":
            molecule = load_rna_structure()
            nt_list = "GAUC"
            force_field_aptamer = "leaprc.RNA.OL3"
        else:  # DNA
            molecule = load_dna_structure()
            nt_list = "GATC"
            force_field_aptamer = "leaprc.DNA.OL21"

        # Choose ligand FF
        if self.molecule_type == "protein":
            force_field_ligand = "leaprc.protein.ff19SB"
            parameterized = True
        elif self.molecule_type == "organic":
            force_field_ligand = "leaprc.gaff2"
            parameterized = False
        else:
            force_field_ligand = "leaprc.lipid21"
            parameterized = False

        log.debug(
            "Forcefields: aptamer=%s ligand=%s parameterized=%s nt_list=%s",
            force_field_aptamer,
            force_field_ligand,
            parameterized,
            nt_list,
        )

        # Template complex with empty aptamer chain + ligand from PDB
        cpx = Complex(
            force_field_aptamer=force_field_aptamer,
            force_field_ligand=force_field_ligand,
            salt_conc=self.salt_conc,
        )
        cpx.add_chain("", molecule)  # empty aptamer chain
        cpx.add_chain_from_pdb(
            pdb_path=pdb_path,
            force_field_aptamer=force_field_aptamer,
            force_field_ligand=force_field_ligand,
            parameterized=parameterized,
        )

        # Ligand-only complex for COM sampling center
        ligand_only = Complex(
            force_field_aptamer=force_field_aptamer,
            force_field_ligand=force_field_ligand,
            salt_conc=self.salt_conc,
        )
        ligand_only.add_chain_from_pdb(
            pdb_path=pdb_path,
            force_field_aptamer=force_field_aptamer,
            force_field_ligand=force_field_ligand,
            parameterized=parameterized,
        )
        ligand_only.build()

        sampler = space.make_sampler(
            ligand_only,
            mode=self.sampler_mode,
            reach=self.reach,
            d_max=self.d_max,
            site_centre=self.site_centre,
            site_radius=self.site_radius,
            probe=self.probe,
            rng=rng,
        )
        rotations = space.NAngles(N_BACKBONE_TORSIONS, rng=rng)

        if self.verbose:
            log.info("MAWS step 1: selecting first nucleotide")
        else:
            log.debug("Step 1 start")

        # ---- Step 1: choose first nucleotide ----
        scored: list[Candidate] = []
        for ntide in nt_list:
            energies = []
            free_E = None
            position = None

            cx = copy.deepcopy(cpx)
            aptamer = cx.aptamer_chain()
            aptamer.create_sequence(ntide)
            cx.build()

            positions0 = cx.positions[:]
            clash = space.ClashFilter(
                cx, aptamer.element, tolerance=self.clash_tolerance
            )

            for _ in range(self.first_chunk_size):
                space.draw_clear_conformation(cx, aptamer, sampler, rotations, clash)

                energy = cx.get_energy()[0]
                if free_E is None or energy < free_E:
                    free_E = energy
                    position = cx.positions[:]
                energies.append(energy)

                cx.positions = positions0[:]

            entropy = entropy_score(energies, beta=self.beta)
            log.debug(
                "Step1 candidate=%s entropy=%s best_E=%s",
                aptamer.alias_sequence,
                entropy,
                free_E,
            )

            scored.append(Candidate(entropy, entropy, free_E, ntide, position[:]))

        beam = select_beam(scored, self.beam)
        log.debug("After step1 beam=%s", [(c.sequence, c.entropy) for c in beam])

        # ---- Steps 2..N: grow sequence (append or prepend) ----
        if self.verbose:
            log.info("MAWS steps 2..N: growing sequence")
        for i in range(1, self.num_nucleotides):
            scored = []
            log.debug("Step%d starting from %s", i + 1, [c.sequence for c in beam])

            for parent in beam:
                best_old_sequence = parent.sequence
                best_old_positions = parent.positions[:]

                for ntide, append in ((n, a) for n in nt_list for a in (True, False)):
                    energies = []
                    free_E = None
                    position = None

                    cx = copy.deepcopy(cpx)
                    aptamer = cx.aptamer_chain()
                    aptamer.create_sequence(best_old_sequence)

                    cx.build()  # cached
                    cx.positions = best_old_positions[:]

                    if append:
                        aptamer.append_sequence(ntide)
                    else:
                        aptamer.prepend_sequence(ntide)

                    cx.rebuild()
                    # Only the aptamer may move: the docking target is rigid.
                    cx.pert_min(
                        size=0.5,
                        atoms=range(aptamer.element[0], aptamer.element[2]),
                        rng=rng,
                    )

                    positions0 = cx.positions[:]
                    clash = space.ClashFilter(
                        cx,
                        space.new_residue_element(aptamer, append=append),
                        tolerance=self.clash_tolerance,
                    )

                    for _ in range(self.second_chunk_size):
                        space.draw_clear_torsions(
                            cx, aptamer, rotations, clash, append=append
                        )

                        energy = cx.get_energy()[0]
                        if free_E is None or energy < free_E:
                            free_E = energy
                            position = cx.positions[:]
                        energies.append(energy)

                        cx.positions = positions0[:]

                    entropy = entropy_score(energies, beta=self.beta)
                    log.debug(
                        "Step%d candidate=%s (%s %s) entropy=%s best_E=%s",
                        i + 1,
                        aptamer.alias_sequence,
                        "append" if append else "prepend",
                        ntide,
                        entropy,
                        free_E,
                    )

                    scored.append(
                        Candidate(
                            entropy,
                            parent.total + entropy,
                            free_E,
                            aptamer.alias_sequence,
                            position[:],
                        )
                    )

            beam = select_beam(scored, self.beam)

        best = beam[0]
        best_entropy = best.entropy
        best_energy = best.energy
        best_sequence = best.sequence
        best_positions = best.positions

        # Optional final PDB artifact
        written_pdb = None
        if output_pdb is not None:
            out = Path(output_pdb)
            if out.exists() and out.is_dir():
                out = out / f"{name}_RESULT.pdb"
            else:
                out.parent.mkdir(parents=True, exist_ok=True)

            result_complex = copy.deepcopy(cpx)
            aptamer = result_complex.aptamer_chain()
            aptamer.create_sequence(best_sequence)
            result_complex.build()
            result_complex.positions = best_positions[:]

            with open(out, "w") as f:
                app.PDBFile.writeModel(
                    result_complex.topology,
                    result_complex.positions,
                    file=f,
                )
            written_pdb = str(out)
            log.debug("Wrote PDB artifact: %s", written_pdb)

        if self.verbose:
            log.info("MAWS run finished: name=%s sequence=%s", name, best_sequence)
        else:
            log.debug("MAWS run finished: name=%s sequence=%s", name, best_sequence)

        return MawsResult(
            sequence=str(best_sequence),
            energy=float(best_energy) if best_energy is not None else float("nan"),
            entropy=float(best_entropy) if best_entropy is not None else float("nan"),
            pdb_path=written_pdb,
            seed=int(seed),
        )
