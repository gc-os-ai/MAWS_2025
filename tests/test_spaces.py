import numpy as np
from MAWS.src.Spaces import Cube, Sphere, NAngles

def test_cube_contains():
    cube = Cube(2.0, centre=(0.0, 0.0, 0.0))
    p = cube.generator()[:3]  # position part
    assert np.all(np.abs(p) <= 1.0)

def test_sphere_radius():
    sph = Sphere(3.0)
    p = sph.generator()[:3]
    assert np.linalg.norm(p) <= 3.0

def test_angles_range():
    angles = NAngles(5).generator()
    assert np.all((0.0 <= angles) & (angles <= 2*np.pi))
