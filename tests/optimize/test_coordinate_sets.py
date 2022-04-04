import pytest
import numpy as np
import geometric.internal
from utils import g2_molecules
from molSimplify.Scripts.rmsd import kabsch_rmsd
from molSimplify.optimize.primitives import (Distance, Angle,
                                             LinearAngle, Dihedral)
from molSimplify.optimize.coordinate_sets import (InternalCoordinates,
                                                  DelocalizedCoordinates)


@pytest.mark.parametrize('name', g2_molecules.keys())
def test_redundant_internals(name):
    if name == 'Si2H6':
        # Skip Si2H6 because of a 0 != 2 pi error
        return
    atoms = g2_molecules[name]['atoms']
    mol = g2_molecules[name]['mol']
    coords_ref = geometric.internal.PrimitiveInternalCoordinates(
        mol, connect=True)

    primitives = []
    for ic in coords_ref.Internals:
        if isinstance(ic, geometric.internal.Distance):
            primitives.append(Distance(ic.a, ic.b))
        elif isinstance(ic, geometric.internal.Angle):
            primitives.append(Angle(ic.a, ic.b, ic.c))
        elif isinstance(ic, geometric.internal.Dihedral):
            primitives.append(Dihedral(ic.a, ic.b, ic.c, ic.d))
        elif isinstance(ic, geometric.internal.OutOfPlane):
            primitives.append(Dihedral(ic.a, ic.b, ic.c, ic.d))
        elif isinstance(ic, geometric.internal.LinearAngle):
            primitives.append(LinearAngle(ic.a, ic.b, ic.c, ic.axis))
        else:
            raise NotImplementedError(f'Internal {type(ic)} not implemented')

    coords = InternalCoordinates(primitives)

    xyzs = atoms.get_positions()
    # Test transformation to internal coordinates
    np.testing.assert_allclose(coords.to_internals(xyzs),
                               coords_ref.calculate(xyzs),
                               atol=1e-8)

    np.testing.assert_allclose(coords.B(xyzs),
                               coords_ref.wilsonB(xyzs),
                               atol=1e-8)

    np.testing.assert_allclose(coords.Ginv(xyzs),
                               coords_ref.GInverse(xyzs),
                               atol=1e-8)

    # Test transformation back to cartesian
    # Try to reconstruct the geometry from a distorted reference
    np.random.seed(4321)
    xyzs_dist = xyzs + 0.1*np.random.randn(*xyzs.shape)
    np.testing.assert_allclose(coords.to_internals(xyzs_dist),
                               coords_ref.calculate(xyzs_dist),
                               atol=1e-8)

    dq_ref = coords_ref.calcDiff(xyzs, xyzs_dist)
    dq = coords.diff_internals(xyzs, xyzs_dist)
    np.testing.assert_allclose(dq, dq_ref, atol=1e-8)
    xyzs2_ref = coords_ref.newCartesian(
        xyzs_dist.flatten(), dq.flatten()).reshape(-1, 3)

    if coords_ref.bork:  # Meaning that the transformation failed:
        with pytest.raises(RuntimeError):
            coords.to_cartesians(dq, xyzs_dist, maxstep=np.inf, tol=1e-6)
        return
    xyzs2 = coords.to_cartesians(dq, xyzs_dist, maxstep=0.1, tol=1e-6)

    assert kabsch_rmsd(xyzs2, xyzs2_ref, translate=True) < 1e-5
    # Test that final geometry is close to original
    assert kabsch_rmsd(xyzs2, xyzs, translate=True) < 1e-5


@pytest.mark.parametrize('name', g2_molecules.keys())
def test_delocalized_internals(name):
    if name == 'Si2H6':
        # Skip Si2H6 because of a 0 != 2 pi error
        return
    atoms = g2_molecules[name]['atoms']
    mol = g2_molecules[name]['mol']
    coords_ref = geometric.internal.DelocalizedInternalCoordinates(
        mol, connect=True, build=True)

    primitives = []
    for ic in coords_ref.Prims.Internals:
        if isinstance(ic, geometric.internal.Distance):
            primitives.append(Distance(ic.a, ic.b))
        elif isinstance(ic, geometric.internal.Angle):
            primitives.append(Angle(ic.a, ic.b, ic.c))
        elif isinstance(ic, geometric.internal.Dihedral):
            primitives.append(Dihedral(ic.a, ic.b, ic.c, ic.d))
        elif isinstance(ic, geometric.internal.OutOfPlane):
            primitives.append(Dihedral(ic.a, ic.b, ic.c, ic.d))
        elif isinstance(ic, geometric.internal.LinearAngle):
            primitives.append(LinearAngle(ic.a, ic.b, ic.c, ic.axis))
        else:
            raise NotImplementedError(f'Internal {type(ic)} not implemented')

    xyzs = atoms.get_positions()
    coords = DelocalizedCoordinates(primitives, xyzs, threshold=1e-6)

    # Obtaining the same eigenvectors is almost impossible due to possible
    # degeneracies. TODO: figure out way to test this nevertheless.

    # From here on out geomeTRIC is forced to use our results.
    coords_ref.Vecs = coords.U

    q = coords.to_internals(xyzs)
    q_ref = coords_ref.calculate(xyzs)
    np.testing.assert_allclose(q, q_ref, atol=1e-8)

    B = coords.B(xyzs)
    B_ref = coords_ref.wilsonB(xyzs)
    np.testing.assert_allclose(B, B_ref, atol=1e-8)
    # Try to reconstruct the geometry from a distorted reference
    np.random.seed(4321)
    xyzs_dist = xyzs + 0.1*np.random.randn(*xyzs.shape)

    dq = coords.diff_internals(xyzs, xyzs_dist)
    dq_ref = coords_ref.calcDiff(xyzs, xyzs_dist)
    np.testing.assert_allclose(dq, dq_ref, atol=1e-8)

    xyzs2 = coords.to_cartesians(dq, xyzs_dist, maxstep=0.1, tol=1e-6)
    # Test that final geometry is close to original
    assert kabsch_rmsd(xyzs2, xyzs, translate=True) < 1e-5
