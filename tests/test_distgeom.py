import numpy as np
from molSimplify.Scripts.distgeom import GetCMDists, CosRule


def test_GetCMDists():
    # Dummy points in 3d
    xyzs = np.array([[0., 0., 0.],
                     [1., 0., 0.],
                     [0., 1., 0.],
                     [0., 1., 1.],
                     [0., 0., 0.]])

    dist_mat = np.sqrt(np.sum(
        (xyzs[np.newaxis, :, :] - xyzs[:, np.newaxis, :])**2, axis=-1))
    # Calculate geometric center
    center = np.mean(xyzs, axis=0)
    # Distances to the center:
    dists_ref = np.sqrt(np.sum((xyzs - center)**2, axis=-1))

    dists, status = GetCMDists(dist_mat, len(xyzs))
    np.testing.assert_allclose(dists, dists_ref)
    # What is status even used for? Always returns True!
    assert status


def test_CosRule(atol=1e-6):
    # Test simple "unit" triangle
    assert abs(CosRule(1., 1., 90) - np.sqrt(2)) < atol
    assert abs(CosRule(1, np.sqrt(2), 45) - 1.) < atol
    assert abs(CosRule(np.sqrt(2), 1., 45) - 1.) < atol
    # Test Pythagoras
    assert abs(CosRule(3., 4., 90) - 5.) < atol
    # Test 30-60-90 triangle
    a = 2.345
    assert abs(CosRule(a, 0.5*a*np.sqrt(3), 30) - 0.5*a) < atol
    assert abs(CosRule(a, 0.5*a, 60) - 0.5*a*np.sqrt(3)) < atol
    assert abs(CosRule(0.5*a, 0.5*a*np.sqrt(3), 90) - a) < atol
    # Test symmetry on 30-60-90
    assert abs(CosRule(0.5*a*np.sqrt(3), a, 30) - 0.5*a) < atol
    assert abs(CosRule(0.5*a, a, 60) - 0.5*a*np.sqrt(3)) < atol
    assert abs(CosRule(0.5*a*np.sqrt(3), 0.5*a, 90) - a) < atol
