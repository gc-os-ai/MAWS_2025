import pathlib
import pytest


from MAWS.src.Structure import Structure, RotTriplet


def test_translate_single():
    s = Structure(["ALA"])
    assert s.translate("ALA") == "ALA"


# def test_translate_multi_default():
#     s = Structure(["ALA", "GLY"])
#     expected = "ALA- -x- -ALA"
#     assert s.translate("ALA GLY ALA") == expected
def test_translate_multi_default():
    s = Structure(["ALA", "GLY"])
    # No alias decorations by default
    assert s.translate("ALA GLY ALA") == "ALA GLY ALA"

@pytest.mark.parametrize(
    "triplet, length, expected",
    [
        ((-1, 4, None), 10, (9, 4, None)),
        ((0, -5, -1), 12, (0, 7, 11)),
    ],
)
def test_normalise_triple(triplet: RotTriplet, length: int, expected: RotTriplet):
    from MAWS.src.Structure import _normalise
    out = tuple(_normalise(i, length) if i is not None else None for i in triplet)
    assert out == expected
