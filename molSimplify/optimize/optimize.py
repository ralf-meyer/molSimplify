import sys
import ase.io
import ase.optimize
import ase.constraints
import geometric
import logging
import numpy as np
import pathlib
from molSimplify.optimize.calculators import get_calculator
from molSimplify.optimize.params import parse_args
from molSimplify.optimize.hessians import compute_guess_hessian
from molSimplify.optimize.coordinates import find_connectivity
from molSimplify.Classes.globalvars import metalslist


logger = logging.getLogger(__name__)


def run_preoptimization(atoms, method, optimizer=ase.optimize.LBFGS):
    method = method.lower()
    atoms.calc = get_calculator(method)
    opt = optimizer(atoms)
    opt.run(fmax=0.01)


def read_molecule(terachem_file):
    # Convert to path. This allows to finds the path of the xyz file later.
    terachem_file = pathlib.Path(terachem_file)
    with open(terachem_file, 'r') as fin:
        lines = fin.readlines()
    for line in lines:
        if line.startswith('coordinates'):
            xyz_path = terachem_file.parent.joinpath(line.split()[1])
            charge = 0.0
            spin = 1.0
        elif line.startswith('charge'):
            charge = float(line.split()[1])
        elif line.startswith('spinmult'):
            spin = float(line.split()[1])
    logger.debug(f'Read terachem file with xyz: {xyz_path}, '
                 f'charge: {charge}, spin: {spin}')
    atoms = ase.io.read(xyz_path)
    # Assign spin and charge to first atom
    q = np.zeros_like(atoms.get_initial_charges())
    q[0] = charge
    atoms.set_initial_charges(q)
    s = np.zeros_like(atoms.get_initial_magnetic_moments())
    # Number of unpaired electrons is multiplicity - 1
    s[0] = spin - 1
    atoms.set_initial_magnetic_moments(s)
    return atoms, xyz_path


def run_preprocessing(args):
    atoms, xyz_path = read_molecule(args.get('input'))
    # Find metal indices
    metal_ind = [i for i in range(len(atoms))
                 if atoms[i].get('symbol') in metalslist]
    # raise error if more than one metal
    if len(metal_ind) > 1:
        raise NotImplementedError('Currently only systems with a single metal '
                                  'atom can be handled.')
    logger.info(f'Metal indices {metal_ind}')
    # Find graph
    connectivity = find_connectivity(atoms)
    preopt = args.get('preopt', False)
    if preopt:
        # Collect indices of frozen atoms starting with a copy of metal_ind
        frozen_inds = [i for i in metal_ind]
        # Freeze connecting atoms
        for bond in connectivity:
            if bond[0] in metal_ind:
                frozen_inds.append(bond[1])
            elif bond[1] in metal_ind:
                frozen_inds.append(bond[0])
        logger.info('Freezing the following atom positions during '
                    f'preoptimization: {frozen_inds}')
        atoms.set_constraint(ase.constraints.FixAtoms(indices=frozen_inds))
        run_preoptimization(atoms, preopt)
        ase.io.write(xyz_path, atoms, plain=True)
        # Remove constraints
        atoms.set_constraint()
    guess_hessian = args.get('guess_hessian', False)
    if guess_hessian:
        H = compute_guess_hessian(atoms, guess_hessian)
        np.savetxt('hessian.txt', H)


def main():
    args, unknown_args = parse_args(sys.argv[1:])
    logging.basicConfig(level=logging.DEBUG)
    logger.info(f'sys.argv[1:]: {sys.argv[1:]}')
    logger.info(f'unknown_args: {unknown_args}')
    # Parse geometric arguments
    geometric_args = geometric.params.parse_optimizer_args(unknown_args)
    # Copy the terachem input file to our args
    args['input'] = geometric_args['input']

    run_preprocessing(args)
    if args.get('guess_hessian', False):
        geometric_args['hessian'] = 'file:./hessian.txt'
    # Call external program
    geometric.optimize.run_optimizer(**geometric_args)


if __name__ == '__main__':
    main()
