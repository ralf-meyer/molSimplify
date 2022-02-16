import os
import pytest
import numpy as np
import ase.build
import ase.io
import ase.collections
import geometric.molecule
import geometric.internal
from molSimplify.optimize.connectivity import (find_connectivity,
                                               find_primitives)


@pytest.mark.parametrize('name', ase.collections.g2.names)
def test_connectivity(tmpdir, name):
    # For some reason the current verion of geomeTRIC uses a covalent radius
    # of zero for Na. Therefore, Na containing molecules have to be skipped.
    if 'Na' in name:
        return
    atoms = ase.build.molecule(name)
    path = os.path.join(tmpdir, 'tmp.xyz')
    ase.io.write(path, atoms, plain=True)
    # geomeTRIC uses a threshold of 1.2 on the unsquared distances.
    # This correspondes to using 1.2^2 in the Billeter et al. alogrithm.
    bonds = find_connectivity(atoms, threshold=1.2**2, connect_fragments=False)
    mol = geometric.molecule.Molecule(path)
    mol.build_topology()
    bonds_ref = list(mol.topology.edges())
    assert bonds == bonds_ref


@pytest.mark.parametrize('name', ase.collections.g2.names)
def test_find_primitives(tmpdir, name):
    atoms = ase.build.molecule(name)
    path = os.path.join(tmpdir, 'tmp.xyz')
    ase.io.write(path, atoms, plain=True)
    # geomeTRIC uses a threshold of 1.2 on the unsquared distances.
    # This correspondes to using 1.2^2 in the Billeter et al. alogrithm.
    bonds = find_connectivity(atoms, threshold=1.2**2, connect_fragments=True)
    # geoemTRIC uses a threshold of cos(theta) = 0.95 for linear angles.
    linear_threshold = np.arccos(0.95) * 180 / np.pi
    bends, linear_bends, torsions, planars = find_primitives(
        atoms.get_positions(), bonds, linear_threshold=linear_threshold,
        planar_threshold=0.95)

    mol = geometric.molecule.Molecule(path)
    mol.build_topology()
    coords_ref = geometric.internal.PrimitiveInternalCoordinates(
        mol, connect=True)
    # Compare bonds
    bonds_ref = [(ic.a, ic.b) for ic in coords_ref.Internals
                 if isinstance(ic, geometric.internal.Distance)]
    assert bonds == bonds_ref
    # Compare dihedrals
    if name != '2-butyne':  # Linear chains not yet implemented
        torsions_ref = [(ic.a, ic.b, ic.c, ic.d) for ic in coords_ref.Internals
                        if isinstance(ic, geometric.internal.Dihedral)]
        assert len(torsions) == len(torsions_ref)
        for t in torsions:
            assert t in torsions_ref or t[::-1] in torsions_ref
    # Compare linear bends. Note that geomeTRIC does not use the
    # "connected to exactly two atoms" rule by Billeter et al for linear bends
    if name != 'ClF3':  # Central atom bound to three neighbors
        linear_bends_ref = [(ic.a, ic.b, ic.c) for ic in coords_ref.Internals
                            if isinstance(ic, geometric.internal.LinearAngle)]
        # Every linear bend appears twice in the reference
        assert linear_bends == linear_bends_ref[::2]
    # Compare planars
    planars_ref = [(ic.a, ic.b, ic.c, ic.d) for ic in coords_ref.Internals
                   if isinstance(ic, geometric.internal.OutOfPlane)]
    assert planars == planars_ref
    # Compare angles. This is where it gets difficult
    bends_ref = [(ic.a, ic.b, ic.c) for ic in coords_ref.Internals
                 if isinstance(ic, geometric.internal.Angle)]
    assert sorted(bends) == sorted(bends_ref)
