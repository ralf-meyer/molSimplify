import os
import pytest
import ase.build
import ase.io
import ase.collections
import geometric.molecule
from molSimplify.optimize.coordinates import find_connectivity


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
    print(atoms.get_positions())
    assert bonds == bonds_ref
