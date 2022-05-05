import sys
import ase.io
import ase.optimize
import ase.constraints
import ase.units
import logging
from molSimplify.optimize.io import read_terachem_input
from molSimplify.optimize.calculators import get_calculator
from molSimplify.optimize.params import parse_args
from molSimplify.optimize.hessian_guess import (compute_hessian_guess,
                                                filter_hessian)
from molSimplify.optimize.connectivity import (find_connectivity)
from molSimplify.optimize.coordinate_sets import get_coordinate_set
from molSimplify.optimize.optimizers import TerachemConvergence, RFO
from molSimplify.Classes.globalvars import metalslist


logger = logging.getLogger(__name__)


def run_preoptimization(atoms, method, optimizer=ase.optimize.LBFGS):
    method = method.lower()
    symbols = atoms.get_chemical_symbols()
    if method in ['mmff94']:
        new_symbols = symbols.copy()
        for i, sym in enumerate(symbols):
            if sym in metalslist:
                new_symbols[i] = 'Si'
        atoms.set_chemical_symbols(new_symbols)
    original_calc = atoms.calc
    atoms.calc = get_calculator(method)
    opt = optimizer(atoms)
    opt.run(fmax=0.01)
    # Reset atom types and calculator
    atoms.set_chemical_symbols(symbols)
    atoms.calc = original_calc


def run_preprocessing(atoms, preopt='xtb', name='molsimp'):
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
    if preopt is not None:
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
        ase.io.write(f'{name}_preopt.xyz', atoms, plain=True)
        # Remove constraints
        atoms.set_constraint()


def run_optimization(atoms, coords='cart', hessian_guess='trivial',
                     hessian_thresh=None, name='molsimp'):

    # Build optimizer with terachem convergence criteria
    class MolSimplifyOpt(TerachemConvergence, RFO):
        pass

    coord_set = get_coordinate_set(atoms, coords)
    H0 = None
    if hessian_guess is not None:
        H0 = compute_hessian_guess(atoms, hessian_guess)
        if hessian_thresh is not None:
            H0 = filter_hessian(H0, hessian_thresh)

    opt = MolSimplifyOpt(atoms, coordinate_set=coord_set, H0=H0,
                         trajectory=f'{name}_optim.traj')
    opt.run()


def main():
    args = parse_args(sys.argv[1:])
    logging.basicConfig(level=logging.DEBUG)
    logger.info(f'sys.argv[1:]: {sys.argv[1:]}')

    atoms = read_terachem_input(args['input'])

    run_preprocessing(args)

    run_optimization(atoms,
                     coords=args['coords'],
                     hessian_guess=args['hessian_guess'],
                     hessian_thresh=args['hessian_thresh'])


if __name__ == '__main__':
    main()
