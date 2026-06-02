"""
Unit tests for maws.complex module.

These tests verify the Complex class's pure Python functionality:
- initialization
- chain management (add_chain, get_chain, aptamer_chain, ligand_chain)
- No external dependencies (AmberTools/OpenMM not required for these tests)

Note: Tests for build(), rebuild(), add_chain_from_pdb(), and rotation
operations require AmberTools/OpenMM and are in test_chain_complex.py.
"""

import pytest

from maws.complex import Complex
from maws.structure import Structure


@pytest.fixture
def simple_structure():
    """Create a simple Structure for testing."""
    residues = ["A", "B", "C"]
    lengths = [10, 15, 20]
    alias = [
        ["A", "A", "A", "A", "A"],
        ["B", "B", "B", "B", "B"],
        ["C", "C", "C", "C", "C"],
    ]
    return Structure(residues, residue_length=lengths, alias=alias)


class TestComplexInit:
    """Tests for Complex initialization."""

    def test_complex_default_init(self):
        """Complex initializes with default force fields."""
        cpx = Complex()
        assert cpx.chains == []
        assert cpx.positions is None
        assert cpx.topology is None
        assert "RNA.OL3" in cpx.build_string
        assert "ff19SB" in cpx.build_string

    def test_complex_custom_force_fields(self):
        """Complex accepts custom force field specifications."""
        cpx = Complex(
            force_field_aptamer="leaprc.DNA.OL21",
            force_field_ligand="leaprc.gaff2",
        )
        assert "DNA.OL21" in cpx.build_string
        assert "gaff2" in cpx.build_string

    def test_complex_default_salt_conc(self):
        """Complex defaults to physiological monovalent salt (0.15 mol/L)."""
        cpx = Complex()
        assert cpx.salt_conc == 0.15

    def test_complex_custom_salt_conc(self):
        """Complex stores a custom salt concentration."""
        cpx = Complex(salt_conc=0.3)
        assert cpx.salt_conc == 0.3

    def test_complex_zero_salt_conc(self):
        """salt_conc=0.0 (documented unscreened mode) is a valid, stored value."""
        cpx = Complex(salt_conc=0.0)
        assert cpx.salt_conc == 0.0


class _FakePrmtop:
    """Records the kwargs passed to createSystem (no OpenMM build required)."""

    def __init__(self):
        self.create_system_kwargs = None

    def createSystem(self, **kwargs):  # noqa: N802 (mirrors OpenMM API)
        self.create_system_kwargs = kwargs
        return "SYSTEM_SENTINEL"


class TestComplexMakeSystem:
    """Tests that salt concentration is threaded into createSystem."""

    def test_make_system_passes_salt_conc(self):
        """_make_system forwards salt_conc as a Debye-Hückel screening term."""
        from openmm import unit

        cpx = Complex(salt_conc=0.2)
        cpx.prmtop = _FakePrmtop()

        system = cpx._make_system()

        assert system == "SYSTEM_SENTINEL"
        assert (
            cpx.prmtop.create_system_kwargs["implicitSolventSaltConc"]
            == 0.2 * unit.molar
        )

    def test_make_system_uses_default_salt_conc(self):
        """Default Complex feeds 0.15 mol/L into createSystem."""
        from openmm import unit

        cpx = Complex()
        cpx.prmtop = _FakePrmtop()

        cpx._make_system()

        assert (
            cpx.prmtop.create_system_kwargs["implicitSolventSaltConc"]
            == 0.15 * unit.molar
        )

    def test_make_system_zero_salt_conc_is_unscreened(self):
        """salt_conc=0.0 forwards 0 mol/L (reproduces the old unscreened build)."""
        from openmm import unit

        cpx = Complex(salt_conc=0.0)
        cpx.prmtop = _FakePrmtop()

        cpx._make_system()

        assert (
            cpx.prmtop.create_system_kwargs["implicitSolventSaltConc"]
            == 0.0 * unit.molar
        )


class TestComplexAddChain:
    """Tests for Complex.add_chain() method."""

    def test_add_first_chain(self, simple_structure):
        """First chain starts at index 0."""
        cpx = Complex()
        cpx.add_chain("A B", simple_structure)

        assert len(cpx.chains) == 1
        assert cpx.chains[0].start == 0
        assert cpx.chains[0].id == 0
        assert cpx.chains[0].length == 25  # A(10) + B(15)

    def test_add_second_chain(self, simple_structure):
        """Second chain starts after first chain's atoms."""
        cpx = Complex()
        cpx.add_chain("A", simple_structure)  # length = 10
        cpx.add_chain("B", simple_structure)  # should start at 10

        assert len(cpx.chains) == 2
        assert cpx.chains[0].start == 0
        assert cpx.chains[1].start == 10
        assert cpx.chains[1].id == 1

    def test_add_empty_chain(self, simple_structure):
        """Empty chain can be added."""
        cpx = Complex()
        cpx.add_chain("", simple_structure)

        assert len(cpx.chains) == 1
        assert cpx.chains[0].length == 0


class TestComplexGetChain:
    """Tests for Complex chain accessor methods."""

    def test_get_chain_by_index(self, simple_structure):
        """get_chain() returns chain by positive index."""
        cpx = Complex()
        cpx.add_chain("A", simple_structure)
        cpx.add_chain("B", simple_structure)

        chain = cpx.get_chain(0)
        assert chain is cpx.chains[0]

        chain = cpx.get_chain(1)
        assert chain is cpx.chains[1]

    def test_get_chain_negative_index(self, simple_structure):
        """get_chain() supports negative indices."""
        cpx = Complex()
        cpx.add_chain("A", simple_structure)
        cpx.add_chain("B", simple_structure)

        chain = cpx.get_chain(-1)
        assert chain is cpx.chains[-1]

    def test_get_chain_no_chains(self):
        """get_chain() raises IndexError when no chains exist."""
        cpx = Complex()

        with pytest.raises(IndexError, match="no chains"):
            cpx.get_chain(0)

    def test_get_chain_out_of_bounds(self, simple_structure):
        """get_chain() raises IndexError for out-of-bounds index."""
        cpx = Complex()
        cpx.add_chain("A", simple_structure)

        with pytest.raises(IndexError):
            cpx.get_chain(5)


class TestComplexConvenienceMethods:
    """Tests for aptamer_chain() and ligand_chain() convenience methods."""

    def test_aptamer_chain(self, simple_structure):
        """aptamer_chain() returns first chain."""
        cpx = Complex()
        cpx.add_chain("A", simple_structure)
        cpx.add_chain("B", simple_structure)

        aptamer = cpx.aptamer_chain()
        assert aptamer is cpx.chains[0]

    def test_aptamer_chain_no_chains(self):
        """aptamer_chain() raises IndexError when no chains."""
        cpx = Complex()

        with pytest.raises(IndexError):
            cpx.aptamer_chain()

    def test_ligand_chain(self, simple_structure):
        """ligand_chain() returns second chain."""
        cpx = Complex()
        cpx.add_chain("A", simple_structure)
        cpx.add_chain("B", simple_structure)

        ligand = cpx.ligand_chain()
        assert ligand is cpx.chains[1]

    def test_ligand_chain_only_one_chain(self, simple_structure):
        """ligand_chain() raises IndexError when only one chain."""
        cpx = Complex()
        cpx.add_chain("A", simple_structure)

        with pytest.raises(IndexError):
            cpx.ligand_chain()


class TestComplexBuildValidation:
    """Tests for Complex.build() validation (without actually running LEaP)."""

    def test_build_no_chains_raises(self):
        """build() raises ValueError when no chains."""
        cpx = Complex()

        with pytest.raises(ValueError, match="Empty Complex"):
            cpx.build()


class TestComplexCacheKey:
    """Tests for _build_cache_key() method."""

    def test_cache_key_deterministic(self, simple_structure):
        """Same inputs produce same cache key."""
        cpx1 = Complex()
        cpx1.add_chain("A B", simple_structure)

        cpx2 = Complex()
        cpx2.add_chain("A B", simple_structure)

        key1 = cpx1._build_cache_key()
        key2 = cpx2._build_cache_key()

        assert key1 == key2

    def test_cache_key_differs_for_different_sequences(self, simple_structure):
        """Different sequences produce different cache keys."""
        cpx1 = Complex()
        cpx1.add_chain("A B", simple_structure)

        cpx2 = Complex()
        cpx2.add_chain("A C", simple_structure)

        key1 = cpx1._build_cache_key()
        key2 = cpx2._build_cache_key()

        assert key1 != key2

    def test_cache_key_differs_for_different_lib_contents(self, tmp_path):
        """Two different ligands sharing the same LIG.lib path must not collide.

        Regression test for the build-cache collision: ``make_lib`` hardcodes the
        residue name ``LIG``, so every protein writes ``LIG.lib`` to the same path
        and yields the same canonical sequence ``"LIG"``. The cache key must depend
        on the .lib *contents*, otherwise the second protein silently reuses the
        first protein's cached topology (dropping/replacing its residues).
        """
        lib = tmp_path / "LIG.lib"

        # First ligand parameterization writes LIG.lib
        lib.write_text('!!index array str\n "LIG"\n# protein A library contents\n')
        struct_a = Structure(["LIG"], residue_length=[100], residue_path=str(tmp_path))
        cpx_a = Complex()
        cpx_a.add_chain("LIG", struct_a)
        key_a = cpx_a._build_cache_key()

        # Second ligand overwrites LIG.lib at the same path with different contents
        lib.write_text('!!index array str\n "LIG"\n# protein B library contents\n')
        struct_b = Structure(["LIG"], residue_length=[200], residue_path=str(tmp_path))
        cpx_b = Complex()
        cpx_b.add_chain("LIG", struct_b)
        key_b = cpx_b._build_cache_key()

        assert key_a != key_b
