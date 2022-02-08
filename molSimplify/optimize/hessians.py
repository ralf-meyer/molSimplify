import numdifftools as nd
from molSimplify.optimize.calculators import (_available_calculators,
                                              get_calculator)


def compute_guess_hessian(atoms, guess_hessian):
    if guess_hessian.lower() in _available_calculators:
        atoms.calc = get_calculator(guess_hessian.lower())
        return numerical_hessian(atoms)


def numerical_hessian(atoms, step=1e-5):
    x0 = atoms.get_positions()

    def fun(x):
        atoms.set_positions(x.reshape(-1, 3))
        return -atoms.get_forces().flatten()
    H = nd.Jacobian(fun, step=step)(x0.flatten())
    atoms.set_positions(x0)
    return 0.5*(H + H.T)
