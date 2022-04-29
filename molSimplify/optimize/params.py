import argparse


def parse_args(*args):

    parser = argparse.ArgumentParser()
    preproc = parser.add_argument_group('preprocessing',
                                        'Options for preprocessing steps.')
    preproc.add_argument('--preopt', type=str, help='')
    preproc.add_argument('--hessian_guess', type=str, help='')
    preproc.add_argument('--hessian_thresh', type=float, help='')
    opt = parser.add_argument_group('optimization',
                                    'Options for the main geometry '
                                    'optimization procedure.')
    opt.add_argument('--optimizer', type=str, help='')

    optimize_args, unknown_args = parser.parse_known_args(*args)
    # Convert to dict following
    # https://docs.python.org/3/library/argparse.html#the-namespace-object
    args_dict = {}
    for key, value in vars(optimize_args).items():
        if value is not None:
            args_dict[key] = value

    return args_dict, unknown_args