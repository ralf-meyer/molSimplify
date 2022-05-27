import os
from molSimplify.Scripts.inparse import parseinputfile
from argparse import Namespace
from pkg_resources import resource_filename, Requirement


def test_parseinputfile_empty():
    defaults = {'skipANN': False, 'oldANN': False,
                'dbvdent': False, 'dbvconns': False,
                'dbvhyb': False, 'dbvlinks': False,
                'rprompt': False, 'rundir': f'{os.getcwd()}/Runs'}

    args = Namespace()
    parseinputfile(args, inputfile_str=' ')

    # Assert defaults are set
    assert args.__dict__ == defaults


def test_parseinputfile_inputfile_kwarg():
    """Test that both methods of calling parseinputfile
    (with an args.i file or inputfile_str kwarg) yield the same result"""
    infile = resource_filename(Requirement.parse(
        "molSimplify"), "tests/inputs/example_1.in")

    args1 = Namespace(i=infile)
    parseinputfile(args1)

    with open(infile, 'r') as fin:
        lines = fin.read()
    args2 = Namespace()
    parseinputfile(args2, inputfile_str=lines)
    # Add 'i' argument for comparison
    args2.i = infile

    assert args1 == args2
