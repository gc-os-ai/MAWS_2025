"""
Tests for :mod:`maws.sampling`.

Every sampler here draws random numbers, so every test that depends on the
values passes a generator built with a fixed seed. That a fixed seed gives a
fixed answer is itself one of the things pinned down: a design run that cannot
be repeated cannot be checked.

The distances are worked out from the van der Waals radii the module lists.
Carbon's is 1.70 Å, and with the default 1.4 Å probe a position is acceptable
only when it is more than 3.10 Å from a carbon atom.
"""

from __future__ import annotations

import warnings
from itertools import islice

import numpy as np
import pytest

from maws.errors import ConfigurationError, SamplingError
from maws.sampling import (
    DEFAULT_VDW,
    Excluder,
    FixedSampler,
    Placement,
    Sphere,
    SurfaceFollowingSampler,
    SurfaceSampler,
    TorsionAngles,
    make_sampler,
)

CARBON_CLEARANCE = 1.70 + 1.4
"""How far a position must be from a lone carbon atom to be acceptable, in Å."""


def lone_carbon(**kwargs) -> Excluder:
    """Return the surface test for a single carbon atom at the origin.

    Parameters
    ----------
    **kwargs
        Passed on to :class:`~maws.sampling.Excluder`, e.g. ``probe``.

    Returns
    -------
    maws.sampling.Excluder
        A test with one atom in it.
    """
    return Excluder(np.zeros((1, 3)), ["C"], **kwargs)


def along_x(distance: float) -> np.ndarray:
    """Return the position `distance` ångström from the origin along x.

    Parameters
    ----------
    distance : float
        How far out, in ångström.

    Returns
    -------
    numpy.ndarray
        Shape ``(3,)``.
    """
    return np.array([distance, 0.0, 0.0])


class TestTorsionAngles:
    """Drawing one angle for each bond that is to be turned."""

    @pytest.mark.parametrize("n", [0, 1, 7])
    def test_it_draws_one_angle_per_bond(self, n, rng):
        """The number asked for is the number handed back."""
        assert TorsionAngles(n, rng=rng).sample().shape == (n,)

    def test_every_angle_is_somewhere_in_a_full_turn(self, rng):
        """Angles are in radians, from zero up to but not including 2π."""
        drawn = TorsionAngles(200, rng=rng).sample()
        assert (drawn >= 0.0).all()
        assert (drawn < 2 * np.pi).all()

    def test_the_same_seed_gives_the_same_angles(self):
        """This is why a generator is passed in rather than made inside."""
        first = TorsionAngles(5, rng=np.random.default_rng(7)).sample()
        second = TorsionAngles(5, rng=np.random.default_rng(7)).sample()
        assert first == pytest.approx(second)

    def test_different_seeds_give_different_angles(self):
        """Two runs that were meant to differ do differ."""
        first = TorsionAngles(5, rng=np.random.default_rng(7)).sample()
        second = TorsionAngles(5, rng=np.random.default_rng(8)).sample()
        assert not np.allclose(first, second)

    def test_a_negative_count_is_rejected(self):
        """There is no such thing as minus three bonds to turn."""
        with pytest.raises(ConfigurationError, match="must not be negative"):
            TorsionAngles(-1)


class TestSphere:
    """Drawing positions from a ball of space."""

    def test_every_position_drawn_is_inside_the_ball(self, rng):
        """A proposal outside the region would defeat the point of having one."""
        ball = Sphere(10.0, np.zeros(3), rng=rng)
        for _ in range(200):
            assert np.linalg.norm(ball.sample().position) <= 10.0

    def test_positions_are_measured_from_the_centre_given(self, rng):
        """A ball centred on the target is where proposals should come from."""
        centre = np.array([100.0, -50.0, 7.0])
        ball = Sphere(3.0, centre, rng=rng)
        assert np.linalg.norm(ball.sample().position - centre) <= 3.0

    def test_the_turn_axis_is_a_direction_of_length_one(self, rng):
        """A rotation is about a line, so its axis carries no length of its
        own."""
        drawn = Sphere(10.0, np.zeros(3), rng=rng).sample()
        assert np.linalg.norm(drawn.axis) == pytest.approx(1.0)

    def test_the_same_seed_gives_the_same_draw(self):
        """A whole run can be repeated by fixing this one seed."""
        first = Sphere(10.0, np.zeros(3), rng=np.random.default_rng(3)).sample()
        second = Sphere(10.0, np.zeros(3), rng=np.random.default_rng(3)).sample()
        assert first.position == pytest.approx(second.position)
        assert first.angle == pytest.approx(second.angle)


class TestExcluder:
    """Deciding whether a proposed position would sit inside the target."""

    def test_a_position_well_away_from_the_target_is_acceptable(self):
        """Nothing is nearby, so there is nothing to overlap."""
        assert lone_carbon().is_clear(along_x(20.0))

    def test_a_position_closer_than_the_atom_plus_the_probe_is_not(self):
        """Carbon's 1.70 Å and the 1.4 Å probe come to 3.10 Å, and 3.0 Å is
        inside that."""
        assert not lone_carbon().is_clear(along_x(3.0))

    def test_a_position_just_beyond_that_is(self):
        """3.2 Å clears the 3.10 Å the two radii add up to."""
        assert lone_carbon().is_clear(along_x(3.2))

    def test_a_larger_probe_keeps_proposals_further_out(self):
        """Rolling a bigger ball over the target thickens the layer it
        excludes."""
        assert lone_carbon().is_clear(along_x(4.0))
        assert not lone_carbon(probe=3.0).is_clear(along_x(4.0))

    def test_positions_and_element_symbols_must_line_up(self):
        """One symbol per atom, or some atom has no size to check against."""
        with pytest.raises(ConfigurationError, match="element symbols"):
            Excluder(np.zeros((3, 3)), ["C", "N"])

    def test_positions_must_be_three_coordinates_per_atom(self):
        """Anything else is not a set of positions."""
        with pytest.raises(ConfigurationError, match=r"shaped \(N, 3\)"):
            Excluder(np.zeros((3, 2)), ["C", "N", "O"])

    def test_a_negative_probe_is_rejected(self):
        """A radius below zero has no meaning."""
        with pytest.raises(ConfigurationError, match="probe must not be negative"):
            Excluder(np.zeros((1, 3)), ["C"], probe=-0.5)


class TestUnknownElements:
    """What happens when an atom's symbol is not one of the listed ones."""

    def unlisted(self, symbol: str) -> None:
        """Forget that `symbol` has already been complained about.

        The complaint is made once per symbol for as long as the program runs,
        so a test that wants to see it has to start from a clean slate.

        Parameters
        ----------
        symbol : str
            The symbol to forget.
        """
        Excluder._warned.discard(symbol)

    def test_an_unlisted_symbol_is_assumed_to_be_the_size_of_carbon(self):
        """It takes up 1.70 Å, so with the default probe it blocks out 3.10 Å
        exactly as a carbon atom would."""
        self.unlisted("Xx")
        with pytest.warns(RuntimeWarning):
            excluder = Excluder(np.zeros((1, 3)), ["Xx"])
        assert DEFAULT_VDW == 1.70
        assert not excluder.is_clear(along_x(CARBON_CLEARANCE - 0.1))
        assert excluder.is_clear(along_x(CARBON_CLEARANCE + 0.1))

    def test_an_unlisted_symbol_is_complained_about(self):
        """A whole target of unrecognised atoms is a mistake worth noticing."""
        self.unlisted("Yy")
        with pytest.warns(RuntimeWarning, match="unknown element"):
            Excluder(np.zeros((4, 3)), ["Yy"] * 4)

    def test_it_is_complained_about_only_once(self):
        """A large target would otherwise bury everything else in warnings."""
        self.unlisted("Zz")
        with pytest.warns(RuntimeWarning):
            Excluder(np.zeros((1, 3)), ["Zz"])
        with warnings.catch_warnings(record=True) as seen:
            warnings.simplefilter("always")
            Excluder(np.zeros((1, 3)), ["Zz"])
        assert seen == []


class TestSurfaceSampler:
    """Proposing positions around a target and throwing away the bad ones."""

    def test_every_placement_it_yields_clears_the_target(self, rng):
        """That is the whole job: a position inside the target is impossible."""
        excluder = lone_carbon()
        sampler = SurfaceSampler(Sphere(20.0, np.zeros(3), rng=rng), excluder)
        for placement in islice(sampler, 50):
            assert excluder.is_clear(placement.position)

    def test_it_can_be_asked_for_one_placement_at_a_time(self, rng):
        """``next(sampler)`` is how a search takes a single guess."""
        sampler = SurfaceSampler(Sphere(20.0, np.zeros(3), rng=rng), lone_carbon())
        assert isinstance(next(sampler), Placement)

    def test_it_can_be_looped_over_like_a_stream(self, rng):
        """A search asks for thousands of guesses with ``islice``."""
        sampler = SurfaceSampler(Sphere(20.0, np.zeros(3), rng=rng), lone_carbon())
        assert len(list(islice(sampler, 5))) == 5

    def test_a_region_that_lies_entirely_inside_the_target_gives_up(self, rng):
        """The ball has radius 5 Å around the origin, and a 100 Å probe on the
        atom at the origin blocks everything within 101.7 Å of it, so no
        proposal can ever pass and the sampler stops instead of looping."""
        blocked = SurfaceSampler(
            Sphere(5.0, np.zeros(3), rng=rng),
            lone_carbon(probe=100.0),
            max_rejections=5,
        )
        with pytest.raises(SamplingError, match="no position clear of the target"):
            blocked.sample()

    def test_around_a_target_draws_from_a_region_that_encloses_it(
        self, empty_system, rng
    ):
        """Every one of the target's own atoms is inside the ball proposals come
        from, so a position on any side of it can be reached."""
        sampler = SurfaceSampler.around(empty_system, rng=rng)
        span = empty_system.chain("ligand").span
        distances = np.linalg.norm(
            empty_system.pose.atoms(span) - sampler.envelope.centre, axis=1
        )
        assert distances.max() < sampler.envelope.radius

    def test_a_larger_reach_gives_a_larger_region(self, empty_system, rng):
        """Reach is how far past the furthest atom the region extends."""
        near = SurfaceSampler.around(empty_system, reach=5.0, rng=rng)
        far = SurfaceSampler.around(empty_system, reach=50.0, rng=rng)
        assert far.envelope.radius == pytest.approx(near.envelope.radius + 45.0)


class TestFixedSampler:
    """Handing out placements from a list written by hand."""

    def placements(self, n: int) -> list[Placement]:
        """Return `n` placements that can be told apart by their position.

        Parameters
        ----------
        n : int
            How many to make.

        Returns
        -------
        list of maws.sampling.Placement
            Placements at 0, 1, ... ångström along x.
        """
        return [
            Placement(along_x(float(index)), np.array([0.0, 0.0, 1.0]), 0.0)
            for index in range(n)
        ]

    def test_it_hands_them_out_in_the_order_they_were_written(self):
        """A test using this knows exactly where the strand will be put."""
        sampler = FixedSampler(self.placements(3))
        handed_out = [float(next(sampler).position[0]) for _ in range(3)]
        assert handed_out == [0.0, 1.0, 2.0]

    def test_it_starts_again_from_the_beginning_by_default(self):
        """The stream never ends, so a search cannot run out of guesses."""
        sampler = FixedSampler(self.placements(2))
        handed_out = [float(next(sampler).position[0]) for _ in range(5)]
        assert handed_out == [0.0, 1.0, 0.0, 1.0, 0.0]

    def test_with_cycling_off_it_stops_when_the_list_runs_out(self):
        """This catches a test that asked for more guesses than it meant to."""
        sampler = FixedSampler(self.placements(2), cycle=False)
        next(sampler)
        next(sampler)
        with pytest.raises(SamplingError, match="have been used"):
            next(sampler)

    def test_an_empty_list_is_rejected(self):
        """A sampler with nothing to hand out could never produce a guess."""
        with pytest.raises(ConfigurationError, match="at least one placement"):
            FixedSampler([])


class TestMakeSampler:
    """Choosing a sampler by name."""

    @pytest.mark.parametrize(
        ("mode", "expected"),
        [("sphere", SurfaceSampler), ("surface-following", SurfaceFollowingSampler)],
    )
    def test_each_name_builds_its_own_sampler(self, empty_system, rng, mode, expected):
        """The two names are the two ways of proposing positions."""
        assert isinstance(make_sampler(empty_system, mode=mode, rng=rng), expected)

    def test_the_default_is_the_ball_around_the_target(self, empty_system, rng):
        """Simple and fast, and right for most runs."""
        assert isinstance(make_sampler(empty_system, rng=rng), SurfaceSampler)

    def test_an_unrecognised_name_is_rejected(self, empty_system, rng):
        """The message lists the two names that work."""
        with pytest.raises(ConfigurationError, match="mode must be"):
            make_sampler(empty_system, mode="spherical", rng=rng)

    def test_a_negative_reach_is_rejected(self, empty_system, rng):
        """A region cannot stop short of the target it is meant to surround."""
        with pytest.raises(ConfigurationError, match="reach must not be negative"):
            make_sampler(empty_system, reach=-1.0, rng=rng)

    def test_a_negative_probe_is_rejected(self, empty_system, rng):
        """The ball rolled over the target has no negative size."""
        with pytest.raises(ConfigurationError, match="probe must not be negative"):
            make_sampler(empty_system, probe=-1.0, rng=rng)
