import os
import pytest
import numpy as np
import ase.build
import ase.collections
import geometric.molecule
import geometric.internal
from molSimplify.Scripts.rmsd import kabsch_rmsd
from molSimplify.optimize.coordinates import (Distance, Angle,
                                              LinearAngle, Dihedral,
                                              InternalCoordinates)


@pytest.mark.parametrize('name', ase.collections.g2.names)
def test_redundant_internals(tmpdir, name):
    atoms = ase.build.molecule(name)
    if len(atoms) == 1 or name == 'Si2H6':
        # Skip single atom systems and Si2H6 because of a 0 != 2 pi error
        return
    path = os.path.join(tmpdir, 'tmp.xyz')
    ase.io.write(path, atoms, plain=True)
    mol = geometric.molecule.Molecule(path)
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





