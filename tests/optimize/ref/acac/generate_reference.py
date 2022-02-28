import subprocess
import shutil
import os
from openbabel import openbabel


def ob_minimize(path, method, frozen_inds):
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats('xyz', 'xyz')

    if method in ['MMFF94']:
        with open(path, 'r') as fin:
            lines = fin.read()
        lines = lines.replace('Fe', 'Si')
        with open(path, 'w') as fout:
            fout.write(lines)

    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, path)

    forcefield = openbabel.OBForceField.FindForceField(method)

    # initialize constraints
    constr = openbabel.OBFFConstraints()
    for i in frozen_inds:
        constr.AddAtomConstraint(i+1)

    print('Forcefield setup successful:', forcefield.Setup(mol, constr))
    print('Energy prior to optimization',
          forcefield.Energy(), forcefield.GetUnit())
    forcefield.SteepestDescent(10000, 1e-8)
    forcefield.ConjugateGradients(10000, 1e-8)
    print('Success in copying:', forcefield.GetCoordinates(mol))
    print('Energy after optimization',
          forcefield.Energy(), forcefield.GetUnit())
    print('Success in writing:', obConversion.WriteFile(mol, path))

    if method in ['MMFF94']:
        with open(path, 'r') as fin:
            lines = fin.read()
        lines = lines.replace('Si', 'Fe')
        with open(path, 'w') as fout:
            fout.write(lines)


def main():
    initial_path = '../../inputs/acac/fe_oct_2_acac_3_s_5_conf_1.xyz'
    frozen_inds = [0, 1, 6, 15, 20, 29, 34]
    with open('xtb.inp', 'w') as fout:
        fout.write('$opt\n')
        fout.write('engine=inertial\n')
        fout.write('$fix\n')
        fout.write('atoms: ')
        fout.write(','.join([str(i+1) for i in frozen_inds]))
        fout.write('\n$end\n')
    subprocess.run(['xtb', '--chrg', '-1', '--uhf', '4', '--opt', 'vtight',
                    '--input', 'xtb.inp', initial_path])
    shutil.move('xtbopt.xyz', 'xtb_opt.xyz')
    for f in ['.xtboptok', 'charges', 'wbo', 'xtbopt.log',
              'xtbrestart', 'xtbtopo.mol', 'xtb.inp']:
        os.remove(f)

    for ff in ['uff', 'mmff94', 'gaff']:
        shutil.copy(initial_path, f'./{ff}_opt.xyz')
        ob_minimize(f'./{ff}_opt.xyz', ff.upper(), frozen_inds)


if __name__ == '__main__':
    main()
