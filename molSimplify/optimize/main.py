import sys
import os
import time
import pkg_resources
import ase.io
import ase.optimize
import ase.constraints
import ase.units
import geometric
import logging
import numpy as np
import pathlib
from molSimplify.optimize.calculators import get_calculator
from molSimplify.optimize.params import parse_args
from molSimplify.optimize.hessian_guess import (compute_hessian_guess,
                                                filter_hessian)
from molSimplify.optimize.connectivity import (find_connectivity,
                                               find_primitives)
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
    atoms.calc = get_calculator(method)
    opt = optimizer(atoms)
    opt.run(fmax=0.01)
    # Reset atom types
    atoms.set_chemical_symbols(symbols)


def read_molecule(terachem_file):
    # Convert to path. This allows to finds the path of the xyz file later.
    terachem_file = pathlib.Path(terachem_file)
    with open(terachem_file, 'r') as fin:
        lines = fin.readlines()
    # Set defaults
    charge = 0.0
    spin = 1.0
    for line in lines:
        if line.startswith('coordinates'):
            xyz_path = terachem_file.parent.joinpath(line.split()[1])
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
    hessian_guess = args.get('hessian_guess', False)
    if hessian_guess:
        H = compute_hessian_guess(atoms, hessian_guess)
        # Transform to Hartree/bohr^2 for geometric
        H = H * ase.units.Bohr**2/ase.units.Hartree
        # Filter small and negative eigenvalues
        H = filter_hessian(H, thresh=args.get('hessian_thresh', 1e-4))
        np.savetxt('hessian.txt', H)


def run_geometric(**kwargs):
    """
    Run geometry optimization, constrained optimization, or
    constrained scan job given arguments from command line.
    """
    if kwargs.get('logIni') is None:
        import geometric.optimize
        logIni = pkg_resources.resource_filename(
            geometric.optimize.__name__, 'config/log.ini')
    else:
        logIni = kwargs.get('logIni')
    logfilename = kwargs.get('prefix')
    # Input file for optimization; QC input file or OpenMM .xml file
    inputf = kwargs.get('input')
    verbose = kwargs.get('verbose', False)
    # Get calculation prefix and temporary directory name
    arg_prefix = kwargs.get('prefix', None)  # prefix for output file
    prefix = (arg_prefix if arg_prefix is not None
              else os.path.splitext(inputf)[0])
    logfilename = prefix + ".log"
    # Create a backup if the log file already exists
    backed_up = geometric.nifty.bak(logfilename)
    import logging.config
    logging.config.fileConfig(logIni, defaults={'logfilename': logfilename},
                              disable_existing_loggers=False)
    # ==============================#
    # | End log file configuration |#
    # ==============================#

    import geometric
    logger.info('geometric-optimize called with the following command line:\n')
    logger.info(' '.join(sys.argv)+'\n')
    geometric.info.print_logo(logger)
    logger.info('-=# \x1b[1;94m geomeTRIC started. '
                'Version: %s \x1b[0m #=-\n' % geometric.__version__)
    if backed_up:
        logger.info('Backed up existing log file: %s -> %s\n' % (
            logfilename, os.path.basename(backed_up)))

    t0 = time.time()

    # Create the params object, containing data to be passed into the optimizer
    params = geometric.params.OptParams(**kwargs)
    params.printInfo()

    # Create "dirname" folder for writing
    dirname = prefix+".tmp"
    if not os.path.exists(dirname):
        os.makedirs(dirname)
    kwargs['dirname'] = dirname

    # Get the Molecule and engine objects needed for optimization
    M, engine = geometric.prepare.get_molecule_engine(**kwargs)

    # Create Work Queue object
    if kwargs.get('port', 0):
        logger.info("Creating Work Queue object for "
                    "distributed Hessian calculation\n")
        geometric.nifty.createWorkQueue(kwargs['port'], debug=verbose > 1)

    # Get initial coordinates in bohr
    coords = M.xyzs[0].flatten() * geometric.nifty.ang2bohr

    # Read in the constraints
    constraints = kwargs.get('constraints', None)  # Constraint input file

    if constraints is not None:
        Cons, CVals = geometric.prepare.parse_constraints(
            M, open(constraints).read())
    else:
        Cons = None
        CVals = None

    # =========================================#
    # | Set up the internal coordinate system |#
    # =========================================#
    # First item in tuple: The class to be initialized
    # Second item in tuple: Whether to connect nonbonded fragments
    # Third item in tuple: Whether to throw in all Cartesians
    #                      (no effect if second item is True)
    CoordSysDict = {
        'cart': (geometric.internal.CartesianCoordinates,
                 False, False),
        'prim': (geometric.internal.PrimitiveInternalCoordinates,
                 True, False),
        'dlc': (geometric.internal.DelocalizedInternalCoordinates,
                True, False),
        'dlc-new': (geometric.internal.DelocalizedInternalCoordinates,
                    True, False),
        'hdlc': (geometric.internal.DelocalizedInternalCoordinates,
                 False, True),
        'tric-p': (geometric.internal.PrimitiveInternalCoordinates,
                   False, False),
        'tric': (geometric.internal.DelocalizedInternalCoordinates,
                 False, False)}
    coordsys = kwargs.get('coordsys', 'tric')
    CoordClass, connect, addcart = CoordSysDict[coordsys.lower()]

    if coordsys == 'dlc-new':
        atoms = ase.io.read(engine.tcin['coordinates'])
        # geomeTRIC uses a threshold of 1.2 on the unsquared distances.
        # This correspondes to using 1.2^2 in the Billeter et al. alogrithm.
        bonds = find_connectivity(atoms, threshold=1.2**2)
        bends, linear_bends, torsions, planars = find_primitives(
            atoms.get_positions(), bonds, planar_method='molsimplify')
        prims = geometric.internal.PrimitiveInternalCoordinates(M)
        prims.Internals = (
            [geometric.internal.Distance(*b) for b in bonds]
            + [geometric.internal.Angle(*a) for a in bends]
            + [geometric.internal.LinearAngle(*a, axis=axis)
               for a in linear_bends for axis in (0, 1)]
            + [geometric.internal.Dihedral(*d) for d in torsions]
            + [geometric.internal.OutOfPlane(*p) for p in planars])
        IC = CoordClass(M, build=False, connect=connect, addcart=addcart,
                        constraints=Cons,
                        cvals=CVals[0] if CVals is not None else None,
                        conmethod=params.conmethod)
        IC.Prims = prims
        xyz = M.xyzs[0].flatten() * geometric.nifty.ang2bohr
        IC.build_dlc(xyz)
    else:
        IC = CoordClass(M, build=True, connect=connect, addcart=addcart,
                        constraints=Cons,
                        cvals=CVals[0] if CVals is not None else None,
                        conmethod=params.conmethod)
    # ========================================#
    # | End internal coordinate system setup |#
    # ========================================#

    # Auxiliary functions (will not do optimization):
    displace = kwargs.get('displace', False)  # Write out the displacements.
    if displace:
        geometric.ic_tools.write_displacements(coords, M, IC, dirname, verbose)
        return

    # Check internal coordinate gradients using finite difference..
    fdcheck = kwargs.get('fdcheck', False)
    if fdcheck:
        IC.Prims.checkFiniteDifferenceGrad(coords)
        IC.Prims.checkFiniteDifferenceHess(coords)
        geometric.ic_tools.check_internal_grad(
            coords, M, IC.Prims, engine, dirname, verbose)
        geometric.ic_tools.check_internal_hess(
            coords, M, IC.Prims, engine, dirname, verbose)
        return

    # Print out information about the coordinate system
    if isinstance(IC, geometric.internal.CartesianCoordinates):
        logger.info("%i Cartesian coordinates being used\n" % (3*M.na))
    else:
        logger.info(f"{len(IC.Internals)} internal coordinates being used"
                    f" (instead of {3*M.na} Cartesians)\n")
    logger.info(IC)
    logger.info("\n")

    if Cons is None:
        # Run a standard geometry optimization
        params.xyzout = prefix+"_optim.xyz"
        progress = geometric.optimize.Optimize(
            coords, M, IC, engine, dirname, params)
    else:
        # Run a single constrained geometry optimization or scan over a grid
        if isinstance(IC, (geometric.internal.CartesianCoordinates,
                           geometric.internal.PrimitiveInternalCoordinates)):
            raise RuntimeError("Constraints only work with delocalized"
                               " internal coordinates")
        Mfinal = None
        for ic, CVal in enumerate(CVals):
            if len(CVals) > 1:
                logger.info(
                    "---=== Scan %i/%i : Constrained Optimization ===---\n" % (
                        ic+1, len(CVals)))
            IC = CoordClass(M, build=True, connect=connect, addcart=addcart,
                            constraints=Cons, cvals=CVal,
                            conmethod=params.conmethod)
            IC.printConstraints(coords, thre=-1)
            if len(CVals) > 1:
                params.xyzout = prefix+"_scan-%03i.xyz" % (ic+1)
                # In the special case of a constraint scan, we write out
                # multiple qdata.txt files
                if params.qdata is not None:
                    params.qdata = 'qdata_scan-%03i.txt' % (ic+1)
            else:
                params.xyzout = prefix+"_optim.xyz"
            if ic == 0:
                progress = geometric.optimize.Optimize(
                    coords, M, IC, engine, dirname, params)
            else:
                progress += geometric.optimize.Optimize(
                    coords, M, IC, engine, dirname, params)
            # update the structure for next optimization in SCAN (by CNH)
            M.xyzs[0] = progress.xyzs[-1]
            coords = progress.xyzs[-1].flatten() * geometric.nifty.ang2bohr
            if Mfinal:
                Mfinal += progress[-1]
            else:
                Mfinal = progress[-1]
            cNames = IC.getConstraintNames()
            cVals = IC.getConstraintTargetVals()
            comment = ', '.join(["%s = %.2f" % (cName, cVal)
                                 for cName, cVal in zip(cNames, cVals)])
            Mfinal.comms[-1] = "Scan Cycle %i/%i ; %s ; %s" % (
                ic+1, len(CVals), comment, progress.comms[-1])
        if len(CVals) > 1:
            Mfinal.write('scan-final.xyz')
            if params.qdata is not None:
                Mfinal.write('qdata-final.txt')
    geometric.optimize.print_citation(logger)
    logger.info("Time elapsed since start of run_optimizer: %.3f seconds\n" % (
        time.time()-t0))
    return progress


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
    if args.get('hessian_guess', False):
        geometric_args['hessian'] = 'file:./hessian.txt'

    opt_method = args.get('optimizer', 'geometric')
    if opt_method == 'geometric':
        # Call external program
        run_geometric(**geometric_args)
    else:
        raise NotImplementedError(f'Optimization method {opt_method} not '
                                  'implemented.')


if __name__ == '__main__':
    main()
