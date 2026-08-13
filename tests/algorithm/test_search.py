"""
Tests for :mod:`maws.search`, the growth algorithm itself.

Every test here runs the real algorithm. What makes that affordable is that
the structures come from :class:`~maws.build.FakeBuilder` and the energies from
:class:`~maws.energy.StubEnergy`, so a whole search finishes in milliseconds
with nothing installed. The arrangement of atoms is not chemistry, but the
control flow, the bookkeeping and the choices are exactly the real ones.
"""

from __future__ import annotations

import numpy as np
import pytest

from maws.build import FakeBuilder
from maws.errors import ConfigurationError, SamplingError
from maws.forcefield import ForceField
from maws.libraries import rna
from maws.regrow import grow_chain
from maws.sampling import Placement, SurfaceSampler
from maws.search import (
    Candidate,
    CandidateScored,
    SearchFinished,
    StepCompleted,
    StepStarted,
    _growth_torsions,
    grow_aptamer,
    stub_energy,
)
from maws.topology import Assembly

ALPHABET_SIZE = 4
"""How many nucleotides the search may choose from at each step."""

DIRECTIONS = 2
"""How many ends of the strand a nucleotide can be added to."""


def run_search(system, builder, **overrides):
    """Run a whole search and return its events as a list.

    Parameters
    ----------
    system : maws.topology.BuiltSystem
        A target with an empty strand beside it.
    builder : maws.build.Builder
        What rebuilds the structure as the strand grows.
    **overrides
        Passed through to :func:`~maws.search.grow_aptamer`, replacing the
        small defaults used here.

    Returns
    -------
    list of maws.search.StepEvent
        Every event the search produced, in order.
    """
    settings = {
        "energy": stub_energy(),
        "builder": builder,
        "n_nucleotides": 3,
        "first_samples": 4,
        "samples": 4,
        "rng": np.random.default_rng(0),
    }
    settings.update(overrides)
    return list(grow_aptamer(system, **settings))


class TestEventStream:
    """What the search reports, and in what order."""

    def test_the_search_ends_with_exactly_one_finished_event(
        self, empty_system, builder
    ):
        """One event marks the end, and it is the last thing produced."""
        events = run_search(empty_system, builder)
        assert isinstance(events[-1], SearchFinished)
        assert sum(isinstance(e, SearchFinished) for e in events) == 1

    def test_every_step_is_announced_before_its_candidates(self, empty_system, builder):
        """A step starts, then its candidates are scored, then it completes."""
        events = run_search(empty_system, builder, n_nucleotides=2)
        kinds = [type(event).__name__ for event in events]
        assert kinds[0] == "StepStarted"
        assert kinds.index("StepCompleted") > kinds.index("CandidateScored")

    def test_one_step_is_announced_per_nucleotide(self, empty_system, builder):
        """A run of n nucleotides takes n steps."""
        events = run_search(empty_system, builder, n_nucleotides=3)
        assert sum(isinstance(e, StepStarted) for e in events) == 3
        assert sum(isinstance(e, StepCompleted) for e in events) == 3

    def test_the_first_step_tries_one_candidate_per_nucleotide(
        self, empty_system, builder
    ):
        """With nothing to extend, only the four nucleotides are tried."""
        events = run_search(empty_system, builder)
        first = [e for e in events if isinstance(e, CandidateScored) and e.step == 1]
        assert len(first) == ALPHABET_SIZE

    def test_later_steps_try_every_nucleotide_at_both_ends(self, empty_system, builder):
        """A strand can be extended at either end, so there are twice as many."""
        events = run_search(empty_system, builder)
        second = [e for e in events if isinstance(e, CandidateScored) and e.step == 2]
        assert len(second) == ALPHABET_SIZE * DIRECTIONS

    def test_the_first_step_starts_from_an_empty_strand(self, empty_system, builder):
        """There is nothing to grow from until the first nucleotide is chosen."""
        events = run_search(empty_system, builder)
        started = next(e for e in events if isinstance(e, StepStarted))
        assert len(started.sequence) == 0

    def test_each_step_starts_from_the_previous_winner(self, empty_system, builder):
        """The strand a step grows from is the one the last step settled on."""
        events = run_search(empty_system, builder)
        completed = [e for e in events if isinstance(e, StepCompleted)]
        started = [e for e in events if isinstance(e, StepStarted)]
        assert started[1].sequence == completed[0].winner.sequence

    def test_every_event_knows_which_step_it_belongs_to(self, empty_system, builder):
        """Step numbers count from one and never skip."""
        events = run_search(empty_system, builder)
        numbers = [e.step for e in events if isinstance(e, StepStarted)]
        assert numbers == [1, 2, 3]

    def test_the_finished_event_repeats_the_last_winner(self, empty_system, builder):
        """Reading only the final event is enough to get the answer."""
        events = run_search(empty_system, builder)
        last_completed = [e for e in events if isinstance(e, StepCompleted)][-1]
        assert events[-1].winner is last_completed.winner


class TestGrowth:
    """How the strand itself changes from step to step."""

    def test_the_strand_gains_exactly_one_nucleotide_per_step(
        self, empty_system, builder
    ):
        """Growth is one residue at a time, which is what makes it tractable."""
        events = run_search(empty_system, builder)
        lengths = [
            len(e.winner.sequence) for e in events if isinstance(e, StepCompleted)
        ]
        assert lengths == [1, 2, 3]

    def test_the_finished_strand_is_as_long_as_asked_for(self, empty_system, builder):
        """A run of n nucleotides produces a strand of n nucleotides."""
        events = run_search(empty_system, builder, n_nucleotides=4)
        assert len(events[-1].winner.sequence) == 4

    def test_a_candidate_keeps_the_nucleotide_it_added(self, empty_system, builder):
        """The candidate says what it added and where, not just the result."""
        events = run_search(empty_system, builder)
        scored = [e for e in events if isinstance(e, CandidateScored) and e.step == 2]
        added = {(e.candidate.token, e.candidate.direction) for e in scored}
        assert added == {
            (token, direction) for token in "GAUC" for direction in ("3prime", "5prime")
        }

    def test_growing_at_the_five_prime_end_puts_the_nucleotide_first(
        self, one_residue_system, builder
    ):
        """Which end grows decides where the nucleotide appears in the strand."""
        grown = grow_chain(
            one_residue_system,
            one_residue_system.pose,
            role="aptamer",
            token="A",
            direction="5prime",
            builder=builder,
        )
        assert str(grown.system.chain("aptamer").sequence) == "A G"

    def test_the_target_never_moves_while_the_strand_grows(
        self, one_residue_system, builder
    ):
        """Rebuilding the strand must not disturb what it is binding to."""
        before = one_residue_system.pose.atoms(one_residue_system.chain("ligand").span)
        grown = grow_chain(
            one_residue_system,
            one_residue_system.pose,
            role="aptamer",
            token="A",
            direction="3prime",
            builder=builder,
        )
        after = grown.pose.atoms(grown.system.chain("ligand").span)
        np.testing.assert_array_equal(before, after)


class TestChoosing:
    """How the search picks a winner from the candidates it scored."""

    def test_the_winner_is_the_candidate_with_the_lowest_score(
        self, empty_system, builder
    ):
        """Lower is better, so the winner is the minimum, never the maximum."""
        events = run_search(empty_system, builder)
        for step in (1, 2, 3):
            scored = [
                e.candidate
                for e in events
                if isinstance(e, CandidateScored) and e.step == step
            ]
            winner = next(
                e.winner
                for e in events
                if isinstance(e, StepCompleted) and e.step == step
            )
            assert winner.score == min(c.score for c in scored)

    def test_the_winner_is_one_of_the_candidates_that_was_scored(
        self, empty_system, builder
    ):
        """The search reports every candidate it considered, including the winner."""
        events = run_search(empty_system, builder)
        scored = [
            e.candidate
            for e in events
            if isinstance(e, CandidateScored) and e.step == 1
        ]
        winner = next(
            e.winner for e in events if isinstance(e, StepCompleted) and e.step == 1
        )
        assert winner in scored

    def test_a_candidate_carries_the_positions_that_gave_its_energy(
        self, empty_system, builder
    ):
        """The reported energy can be reproduced from the reported positions."""
        events = run_search(empty_system, builder)
        candidate = events[-1].winner
        model = stub_energy()(candidate.system)
        assert model.evaluate(candidate.pose) == pytest.approx(candidate.energy)

    def test_a_candidates_positions_fit_its_own_structure(self, empty_system, builder):
        """Each candidate is a different structure, with positions to match."""
        events = run_search(empty_system, builder)
        for event in events:
            if isinstance(event, CandidateScored):
                assert event.candidate.pose.n_atoms == event.candidate.system.n_atoms


class TestRepeatability:
    """A run with a fixed seed can be reproduced exactly."""

    def test_the_same_seed_gives_the_same_strand(self, empty_system, builder):
        """Two runs with the same seed agree, so a result can be checked."""
        first = run_search(empty_system, builder, rng=np.random.default_rng(7))
        second = run_search(empty_system, builder, rng=np.random.default_rng(7))
        assert str(first[-1].winner.sequence) == str(second[-1].winner.sequence)

    def test_the_same_seed_gives_the_same_scores(self, empty_system, builder):
        """Repeatability covers the numbers, not only the answer."""
        first = run_search(empty_system, builder, rng=np.random.default_rng(7))
        second = run_search(empty_system, builder, rng=np.random.default_rng(7))
        assert first[-1].winner.score == second[-1].winner.score

    def test_different_seeds_explore_differently(self, empty_system, builder):
        """Without a fixed seed the search is genuinely random."""
        first = run_search(empty_system, builder, rng=np.random.default_rng(1))
        second = run_search(empty_system, builder, rng=np.random.default_rng(2))
        assert first[-1].winner.energy != second[-1].winner.energy


class TestSettings:
    """Arguments that change what the search does."""

    def test_the_alphabet_limits_which_nucleotides_are_tried(
        self, empty_system, builder
    ):
        """Restricting the alphabet restricts what can appear in the strand."""
        events = run_search(empty_system, builder, alphabet="G", n_nucleotides=2)
        assert set(str(events[-1].winner.sequence).split()) == {"G"}

    def test_a_dna_design_uses_the_dna_nucleotides(self, builder):
        """The nucleotides come from the force field unless overridden."""
        forcefield = ForceField.for_target("DNA", "protein")
        from maws.libraries import dna

        system = builder.build(
            Assembly().with_aptamer(dna()).with_ligand_stub(20), forcefield
        )
        events = run_search(system, builder, n_nucleotides=1)
        assert set(str(events[-1].winner.sequence).split()) <= set("GATC")

    def test_the_first_step_can_sample_more_than_later_ones(
        self, empty_system, builder
    ):
        """The first step also searches over where to put the strand."""
        events = run_search(
            empty_system, builder, n_nucleotides=2, first_samples=8, samples=2
        )
        assert isinstance(events[-1], SearchFinished)


class TestStopping:
    """Reading the events as a stream, rather than waiting for the end."""

    def test_a_search_can_be_abandoned_part_way(self, empty_system, builder):
        """Because progress is a stream, a caller can simply stop reading."""
        seen = []
        for event in grow_aptamer(
            empty_system,
            energy=stub_energy(),
            builder=builder,
            n_nucleotides=10,
            first_samples=4,
            samples=4,
            rng=np.random.default_rng(0),
        ):
            seen.append(event)
            if isinstance(event, StepCompleted) and event.step == 2:
                break

        assert not any(isinstance(e, SearchFinished) for e in seen)
        assert len(seen[-1].winner.sequence) == 2


class TestShapesThatBuryTheStrandAreNeverScored:
    """A shape with the strand inside the target is not a shape at all.

    The sampler proposes where the middle of the strand goes, and rejects a
    middle that lands inside the target. That is not the same question as
    whether the strand clears it: a nucleotide is about ten ångström across, so
    it can have its middle in clear space and its far end buried. The search
    therefore asks the sampler again once the atoms are actually somewhere.
    """

    class RefuseEverything:
        """A sampler that proposes freely and then rejects every shape.

        Stands in for a target the strand can find no room beside, without
        having to construct one.
        """

        def sample(self):
            """Return a placement at the origin, turned not at all."""
            return Placement(
                position=np.zeros(3), axis=np.array([0.0, 0.0, 1.0]), angle=0.0
            )

        def accepts(self, points):
            """Reject every shape, whatever the atoms did.

            Parameters
            ----------
            points : numpy.ndarray
                Shape ``(N, 3)``, in ångström. Ignored.

            Returns
            -------
            bool
                Always False.
            """
            return False

        def __next__(self):
            return self.sample()

        def __iter__(self):
            return self

    def test_a_run_with_no_room_at_all_says_so(self, empty_system, builder):
        """Rather than scoring nothing, or scoring the buried shapes anyway."""
        with pytest.raises(SamplingError, match="clear of the target"):
            run_search(empty_system, builder, sampler=self.RefuseEverything())

    def test_the_message_says_how_many_shapes_were_tried(self, empty_system, builder):
        """So the reader can tell "no room" from "asked for too many"."""
        with pytest.raises(SamplingError, match="out of 80 tried"):
            run_search(
                empty_system, builder, sampler=self.RefuseEverything(), first_samples=4
            )

    def test_every_shape_scored_keeps_the_strand_clear_of_the_target(
        self, empty_system, builder, rng
    ):
        """The check the search applies, applied again from the outside.

        The winner's positions are the best of the shapes that were scored, so
        if any buried shape had been let through this is where it would show.
        """
        sampler = SurfaceSampler.around(empty_system, rng=rng)
        events = run_search(empty_system, builder, sampler=sampler)
        winner = events[-1].winner
        strand = winner.system.chain("aptamer").span
        assert sampler.excluder.all_clear(winner.pose.atoms(strand))


class TestDegreesOfFreedom:
    """Which bonds a step is allowed to turn, once the strand has grown.

    The private helper is used directly because this is the one place where
    the two ends of the strand can be compared side by side. Going through a
    whole search would only show the consequence — one end winning more often
    than it should — with no way to say why.
    """

    def _grown(self, builder, direction):
        """Return a three-nucleotide strand grown once more at one end.

        Parameters
        ----------
        builder : maws.build.Builder
            What builds the structures.
        direction : {"3prime", "5prime"}
            Which end to add the fourth nucleotide at.

        Returns
        -------
        maws.pose.ChainView
            The strand after growing, ready to be asked for its bonds.
        """
        system = builder.build(
            Assembly().with_aptamer(rna(), "G A U"),
            ForceField.for_target("RNA", "protein"),
        )
        grown = grow_chain(
            system,
            system.pose,
            role="aptamer",
            token="C",
            direction=direction,
            builder=builder,
        )
        return grown.system.chain("aptamer")

    @pytest.mark.parametrize("direction", ["3prime", "5prime"])
    def test_a_step_gets_as_many_bonds_as_it_asked_for(self, builder, direction):
        """Asking for four bonds to turn gives four, at either end."""
        chain = self._grown(builder, direction)
        assert len(_growth_torsions(chain, direction, 4)) == 4

    @pytest.mark.parametrize("direction", ["3prime", "5prime"])
    def test_every_bond_a_step_turns_moves_the_new_residue(self, builder, direction):
        """No bond may leave the new residue alone.

        The residue just added is the only thing this step is placing. A bond
        that does not move any of its atoms spends one of the step's random
        angles on shuffling the part of the strand that was already settled
        against the target.
        """
        chain = self._grown(builder, direction)
        newest = chain.residue(-1 if direction == "3prime" else 0)
        atoms = set(newest.span)
        for torsion in _growth_torsions(chain, direction, 4):
            assert atoms & set(torsion.moving), (
                f"bond {torsion.pivot}-{torsion.bond} moves nothing of the "
                f"residue just added at the {direction} end"
            )

    def test_both_ends_are_offered_the_same_number_of_bonds(self, builder):
        """Neither end may be handed fewer bonds to turn than the other.

        The search decides which end to grow at by comparing the scores the two
        produce. An end offered fewer bonds would try a smaller set of shapes,
        find a worse best one, and lose that comparison for a reason that has
        nothing to do with the molecule.
        """
        three = len(_growth_torsions(self._grown(builder, "3prime"), "3prime", 4))
        five = len(_growth_torsions(self._grown(builder, "5prime"), "5prime", 4))
        assert three == five


class TestRejectedArguments:
    """Arguments the search refuses, and what it says about them."""

    def test_a_strand_of_no_nucleotides_is_rejected(self, empty_system, builder):
        """There is no such thing as an aptamer of length zero."""
        with pytest.raises(ConfigurationError, match="at least 1"):
            run_search(empty_system, builder, n_nucleotides=0)

    def test_a_strand_that_already_has_residues_is_rejected(
        self, one_residue_system, builder
    ):
        """The search builds a strand from nothing; it does not extend one."""
        with pytest.raises(ConfigurationError, match="already holds"):
            run_search(one_residue_system, builder)

    def test_an_empty_alphabet_is_rejected(self, empty_system, builder):
        """With no nucleotides to choose from there is nothing to score."""
        with pytest.raises(ConfigurationError, match="alphabet is empty"):
            run_search(empty_system, builder, alphabet="")


class TestCandidateValue:
    """The record of one scored way of extending the strand."""

    def test_a_candidate_cannot_be_changed_after_it_is_made(self):
        """Candidates are kept and compared, so they are frozen."""
        builder = FakeBuilder()
        forcefield = ForceField.for_target("RNA", "protein")
        system = builder.build(
            Assembly().with_aptamer(rna(), "G").with_ligand_stub(5), forcefield
        )
        candidate = Candidate(
            sequence=system.chain("aptamer").sequence,
            token="G",
            direction="3prime",
            score=-0.1,
            energy=-10.0,
            system=system,
            pose=system.pose,
        )
        with pytest.raises(AttributeError):
            candidate.score = 0.0
