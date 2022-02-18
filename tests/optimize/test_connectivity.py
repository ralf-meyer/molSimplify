import pytest
import numpy as np
import geometric.internal
from utils import g2_molecules
from molSimplify.optimize.connectivity import (find_connectivity,
                                               find_primitives)


@pytest.mark.parametrize('system', g2_molecules)
def test_connectivity(system):
    # For some reason the current verion of geomeTRIC uses a covalent radius
    # of zero for Na. Therefore, Na containing molecules have to be skipped.
    name = system['name']
    if 'Na' in name:
        return
    atoms = system['atoms']
    mol = system['mol']
    mol.build_topology()
    # geomeTRIC uses a threshold of 1.2 on the unsquared distances.
    # This correspondes to using 1.2^2 in the Billeter et al. alogrithm.
    bonds = find_connectivity(atoms, threshold=1.2**2, connect_fragments=False)

    bonds_ref = list(mol.topology.edges())
    assert bonds == bonds_ref


@pytest.mark.parametrize('system', g2_molecules)
def test_find_primitives(system):
    name = system['name']
    atoms = system['atoms']
    mol = system['mol']
    mol.build_topology()
    # geomeTRIC uses a threshold of 1.2 on the unsquared distances.
    # This correspondes to using 1.2^2 in the Billeter et al. alogrithm.
    bonds = find_connectivity(atoms, threshold=1.2**2, connect_fragments=True)
    # geoemTRIC uses a threshold of cos(theta) = 0.95 for linear angles.
    linear_threshold = np.arccos(0.95) * 180 / np.pi
    bends, linear_bends, torsions, planars = find_primitives(
        atoms.get_positions(), bonds, linear_threshold=linear_threshold,
        planar_threshold=0.95)

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
    # Compare planars. For 'C3H4_C2v the methods do not agree on the ordering
    # and for 'methylenecyclopropane' geomeTRIC does not find the second OOP.
    # Second list of exclusions is for a bug / implementation detail how
    # non-planar systems are recognized in geomeTRIC (dot product of the normal
    # vectors of angles a-i-j and i-j-k instead of a-i-j and a-j-k)
    if (name not in ['methylenecyclopropane', 'C3H4_C2v'] and
            name not in ['C3H9C', 'ClF3', 'CH2NHCH2', 'CH2OCH2', 'CH3CONH2',
                         'C3H7', 'C5H8']):
        planars_ref = [(ic.a, ic.b, ic.c, ic.d) for ic in coords_ref.Internals
                       if isinstance(ic, geometric.internal.OutOfPlane)]
        assert planars == planars_ref
        # Compare angles.
        bends_ref = [(ic.a, ic.b, ic.c) for ic in coords_ref.Internals
                     if isinstance(ic, geometric.internal.Angle)]
        assert sorted(bends) == sorted(bends_ref)
