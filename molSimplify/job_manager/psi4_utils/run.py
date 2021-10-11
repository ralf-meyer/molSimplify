import psi4
import argparse
import os
import numpy as np
import shutil
import subprocess
import iodata
import json
import glob
from os.path import expanduser
from molSimplify.job_manager.psi4_utils.molden2psi4wfn import tcmolden2psi4wfn_ao_mapping
from molSimplify.job_manager.psi4_utils.molden2psi4wfn_spherical import tcmolden2psi4wfn_ao_mapping_spherical
from pkg_resources import resource_filename, Requirement


def lacvps(mol, role):
    '''
    Define the LACVP* basis set.
    lanl2dz for metals and 6-31g* for others.
    '''
    basstrings = {}
    mol.set_basis_all_atoms("6-31g*", role=role)
    mol.set_basis_by_symbol("fe", "lanl2dz", role=role)
    mol.set_basis_by_symbol("co", "lanl2dz", role=role)
    mol.set_basis_by_symbol("cr", "lanl2dz", role=role)
    mol.set_basis_by_symbol("mn", "lanl2dz", role=role)
    mol.set_basis_by_symbol("mo", "lanl2dz", role=role)
    mol.set_basis_by_symbol("tc", "lanl2dz", role=role)
    mol.set_basis_by_symbol("ru", "lanl2dz", role=role)
    mol.set_basis_by_symbol("rh", "lanl2dz", role=role)
    mol.set_basis_by_symbol("I", "lanl2dz", role=role)
    mol.set_basis_by_symbol("Br", "lanl2dz", role=role)
    mol.set_basis_by_symbol("hf", "lanl2dz", role=role)
    return basstrings


def get_molecule(xyzfile, charge, spin, sym='c1'):
    '''
    Assemble a molecule object from xyzfile, charge and spin.
    '''
    wholetext = "%s %s\n" % (charge, spin)
    if os.path.isfile(xyzfile):
        with open(xyzfile, "r") as fo:
            natoms = int(fo.readline().split()[0])
            fo.readline()
            for ii in range(natoms):
                wholetext += fo.readline()
    wholetext += "\nsymmetry %s\nnoreorient\nnocom\n"%sym
    mol = psi4.geometry("""%s""" % wholetext)
    return mol


def setup_dft_parameters(psi4_config):
    psi4.set_memory(psi4_config["memory"])
    psi4.set_num_threads(psi4_config["num_threads"])
    if psi4_config["basis"] == "lacvps":
        psi4.qcdb.libmintsbasisset.basishorde['LACVPS'] = lacvps
        psi4.set_options({"puream": False})
    elif psi4_config["basis"] == "def2-tzvp":
        psi4.set_options({"puream": True})
    else:
        psi4.set_options({"puream": True})
        # pass
        # raise ValueError("Only lacvps is supported!")
    psi4.set_options({
        'reference': psi4_config["ref"],
        "DF_SCF_GUESS": False,
        "scf_type": "df",
        "dft_pruning_scheme": "robust",
        "basis": psi4_config["basis"],
        "DFT_BASIS_TOLERANCE": 1e-10,
        "INTS_TOLERANCE": 1e-10,
        "PRINT_MOS": False,
        "dft_spherical_points": 590,
        "dft_radial_points": 99,
        # "SOSCF": True,
        # "DAMPING_PERCENTAGE": 20,
        # "BASIS_GUESS": True,
        "guess": "read", })


def ensure_dir(dirpath):
    if not os.path.isdir(dirpath):
        os.makedirs(dirpath)


def check_sucess():
    success = False
    with open("output.dat", "r") as fo:
        txt = "".join(fo.readlines())
    if 'Computation Completed' in txt:
        success = True
    return success


def b3lyp_hfx():
    b3lyp_d = {}
    for hfx in [0, 5, 10, 15, 20, 25, 30, 35, 40]:
        hfx_func = {
            "name": "hfx_func",
            "x_functionals": {
                "GGA_X_B88": {"alpha": 0.9*(1-hfx*0.01)},
                "LDA_X": {"alpha": 0.1*(1-hfx*0.01)}
                    },
            "x_hf": {"alpha": hfx*0.01},
            "c_functionals": {
                "GGA_C_LYP": {"alpha": 0.81 },
                "LDA_C_VWN_RPA": {"alpha": 0.19 }
            }
        }
        b3lyp_d["b3lyp_" + str(hfx)] = hfx_func
    return b3lyp_d


def run_b3lyp(psi4_config, rundir="./b3lyp"):
    b3lyp_d = b3lyp_hfx()
    home = expanduser("~")
    psi4_scr = './'
    filename = "output"
    basedir = os.getcwd()
    d = json.load(open(psi4_config["charge-spin-info"], "r"))
    psi4_config.update(d)
    ensure_dir(rundir)
    shutil.copyfile("geo.xyz", rundir + '/geo.xyz')
    sym = 'c1' if 'sym' not in psi4_config else psi4_config['sym']
    mol = get_molecule(psi4_config["xyzfile"], psi4_config["charge"], psi4_config["spin"], sym)
    setup_dft_parameters(psi4_config)
    pid = str(os.getpid())
    if os.path.isfile(psi4_config["moldenfile"]):
        shutil.copyfile(psi4_config["moldenfile"], rundir + '/'+ psi4_config["moldenfile"])
        os.chdir(rundir)
        psi4.core.set_output_file(filename + '.dat', False)
        ## 1-step SCF
        psi4.set_options({
            "maxiter": 5,
            "D_CONVERGENCE": 1e5,
            "E_CONVERGENCE": 1e5,
            "fail_on_maxiter": False})
        if psi4_config["basis"] == "def2-tzvp":
            psi4.set_options({"basis": "def2-sv(p)"})
        if "b3lyp_hfx" in psi4_config and psi4_config["b3lyp_hfx"] != 20:
            print("customized b3lyp with different HFX: ", psi4_config["b3lyp_hfx"])
            functional = "b3lyp_" + str(psi4_config["b3lyp_hfx"])
            e, wfn = psi4.energy("scf", dft_functional=b3lyp_d[functional],  molecule=mol, return_wfn=True)
        else:
            e, wfn = psi4.energy('b3lyp', molecule=mol, return_wfn=True)
        wfn.to_file("wfn-1step.180")
        ## Get converged WFN
        d_molden = iodata.molden.load_molden(psi4_config["moldenfile"])
        restricted = True if any(x in psi4_config["ref"] for x in ["r", "R"]) else False
        if not psi4_config["basis"] == "def2-tzvp":
            Ca, Cb, mapping = tcmolden2psi4wfn_ao_mapping(d_molden, restricted=restricted)
        else:
            Ca, Cb, mapping = tcmolden2psi4wfn_ao_mapping_spherical(d_molden, restricted=restricted)
        wfn_minimal_np = np.load("wfn-1step.180.npy", allow_pickle=True)
        wfn_minimal_np[()]['matrix']["Ca"] = Ca
        if not restricted:
            wfn_minimal_np[()]['matrix']["Cb"] = Cb
        else:
            wfn_minimal_np[()]['matrix']["Cb"] = Ca
        np.save("wfn-1step-tc.180.npy", wfn_minimal_np)
        ## Copy wfn file to the right place with a right name
        pid = str(os.getpid())
        targetfile = psi4_scr + filename + '.default.' + pid + '.180.npy'
        shutil.copyfile("wfn-1step-tc.180.npy", targetfile)
        # Final scf---
        psi4.set_options({
            "maxiter": 50,
            "D_CONVERGENCE": 3e-6,
            "E_CONVERGENCE": 3e-6,
            "fail_on_maxiter": True})
    else:
        os.chdir(rundir)
        psi4.core.set_output_file(filename + '.dat', False)
        print("Warning: no Molden file is used to initialize this calculation!")
        psi4.set_options({
            "maxiter": 250,
            "D_CONVERGENCE": 3e-5,
            "E_CONVERGENCE": 3e-5,
            "fail_on_maxiter": True})
        if psi4_config["basis"] == "def2-tzvp":
            psi4.set_options({"basis": "def2-sv(p)"})
    sucess = False
    try:
        if "b3lyp_hfx" in psi4_config and psi4_config["b3lyp_hfx"] != 20:
            print("customized b3lyp with different HFX: ", psi4_config["b3lyp_hfx"])
            functional = "b3lyp_" + str(psi4_config["b3lyp_hfx"])
            e, wfn = psi4.energy("scf", dft_functional=b3lyp_d[functional],  molecule=mol, return_wfn=True)
        else:
            e, wfn = psi4.energy('b3lyp', molecule=mol, return_wfn=True)
        wfn.to_file("wfn.180")
        sucess = True
    except:
        print("This calculation does not converge.")
    if psi4_config["basis"] == "def2-tzvp" and sucess:
        psi4.set_options({"basis": "def2-tzvp", "maxiter": 200,})
        try:
            if "b3lyp_hfx" in psi4_config and psi4_config["b3lyp_hfx"] != 20:
                print("customized b3lyp with different HFX: ", psi4_config["b3lyp_hfx"])
                functional = "b3lyp_" + str(psi4_config["b3lyp_hfx"])
                e, wfn = psi4.energy("scf", dft_functional=b3lyp_d[functional],  molecule=mol, return_wfn=True)
            else:
                e, wfn = psi4.energy('b3lyp', molecule=mol, return_wfn=True)
            wfn.to_file("wfn.180")
        except:
            print("This calculation does not converge.")
    success = check_sucess()
    for filename in os.listdir("./"):
        if ("psi." in filename) or ("default" in filename):
            print("removing: :", filename)
            os.remove(filename)
    os.chdir(basedir)
    return success


def run_general(psi4_config, functional):
    b3lyp_d = b3lyp_hfx()
    psi4_scr = './'
    filename = "output"
    basedir = os.getcwd()
    rundir = "./" + functional
    d = json.load(open(psi4_config["charge-spin-info"], "r"))
    psi4_config.update(d)
    shutil.copyfile("geo.xyz", functional + '/geo.xyz')
    ensure_dir(rundir)
    os.chdir(rundir)
    psi4.core.set_output_file(filename + '.dat', False)
    sym = 'c1' if 'sym' not in psi4_config else psi4_config['sym']
    mol = get_molecule(psi4_config["xyzfile"], psi4_config["charge"], psi4_config["spin"], sym)
    setup_dft_parameters(psi4_config)
    ## Copy wfn file to the right place with a right name---
    pid = str(os.getpid())
    targetfile = psi4_scr + filename + '.default.' + pid + '.180.npy'
    shutil.copyfile(psi4_config["wfnfile"], targetfile)
    ## Final scf---
    psi4.set_options({
        "maxiter": 50,
        "D_CONVERGENCE": 3e-5,
        "E_CONVERGENCE": 3e-5,
        "fail_on_maxiter": True})
    if True:
        if (not functional in b3lyp_d) and (not "hfx_" in functional):
            e, wfn = psi4.energy(functional, molecule=mol, return_wfn=True)
        elif "hfx_" in functional:
            basefunc, hfx = functional.split("_")[0], int(functional.split("_")[-1])
            print("HFX sampling: ", basefunc, hfx)
            e, wfn = psi4.energy("scf", dft_functional=get_hfx_functional(basefunc, hfx),  molecule=mol, return_wfn=True)
        else:
            print("customized b3lyp with different HFX: ", functional)
            e, wfn = psi4.energy("scf", dft_functional=b3lyp_d[functional],  molecule=mol, return_wfn=True)
        wfn.to_file("wfn.180")
    else:
        print("This calculation does not converge.")
    success = check_sucess()
    for filename in os.listdir("./"):
        if ("psi." in filename) or ("default" in filename):
            print("removing: :", filename)
            os.remove(filename)
    os.chdir(basedir)
    return success


def run_general_hfx(psi4_config, functional, hfx, wfn):
    psi4_scr = './'
    filename = "output"
    basedir = os.getcwd()
    rundir = "./" + functional + "-%d"% hfx
    d = json.load(open(psi4_config["charge-spin-info"], "r"))
    psi4_config.update(d)
    shutil.copyfile("geo.xyz", rundir + '/geo.xyz')
    ensure_dir(rundir)
    os.chdir(rundir)
    psi4.core.set_output_file(filename + '.dat', False)
    sym = 'c1' if 'sym' not in psi4_config else psi4_config['sym']
    mol = get_molecule(psi4_config["xyzfile"], psi4_config["charge"], psi4_config["spin"], sym)
    setup_dft_parameters(psi4_config)
    ## Copy wfn file to the right place with a right name---
    pid = str(os.getpid())
    targetfile = psi4_scr + filename + '.default.' + pid + '.180.npy'
    if not os.path.isfile(wfn):
        print("Previous calculation failed... This one is skipped.")
        return False
    shutil.copyfile(wfn, targetfile)
    ## Final scf---
    psi4.set_options({
        "maxiter": 50,
        "D_CONVERGENCE": 3e-5,
        "E_CONVERGENCE": 3e-5,
        "fail_on_maxiter": True})
    try:
        e, wfn = psi4.energy("scf", molecule=mol, return_wfn=True, dft_functional=get_hfx_functional(functional, hfx))
        wfn.to_file("wfn.180")
    except:
        print("This calculation does not converge.")
    success = check_sucess()
    for filename in os.listdir("./"):
        if ("psi." in filename) or ("default" in filename):
            print("removing: :", filename)
            os.remove(filename)
    os.chdir(basedir)
    return success


def get_hfx_functional(functional, hfx):
    fmap = {"tpss": "TPSS", "scan": "SCAN", "m06-l": "M06_L", "mn15-l": "MN15_L"}
    if functional == "bp86":
        hfx_func = {
            "name": "hfx_func",
            "x_functionals": {"GGA_X_B88": {"alpha": 1-hfx*0.01}},
            "x_hf": {"alpha": hfx*0.01},
            "c_functionals": {"GGA_C_P86": {}}
        }
    elif functional == "blyp":
        hfx_func = {
            "name": "hfx_func",
            "x_functionals": {"GGA_X_B88": {"alpha": 1-hfx*0.01}},
            "x_hf": {"alpha": hfx*0.01},
            "c_functionals": {"GGA_C_LYP": {}}
        }
    elif functional == "b3lyp":
        hfx_func = {
            "name": "hfx_func",
            "x_functionals": {
                "GGA_X_B88": {"alpha": 0.9*(1-hfx*0.01)},
                "LDA_X": {"alpha": 0.1*(1-hfx*0.01)}
                    },
            "x_hf": {"alpha": hfx*0.01},
            "c_functionals": {
                "GGA_C_LYP": {"alpha": 0.81 },
                "LDA_C_VWN_RPA": {"alpha": 0.19 }
            }
        }
    elif functional == "pbe":
        hfx_func = {
            "name": "hfx_func",
            "x_functionals": {"GGA_X_PBE": {"alpha": 1-hfx*0.01}},
            "x_hf": {"alpha": hfx*0.01},
            "c_functionals": {"GGA_C_PBE": {}}
        }
    elif functional in ["m06-l", "mn15-l", "scan", "tpss"]:
        mega = "" if "PBE" in functional else "M"
        hfx_func = {
            "name": "hfx_func",
            "x_functionals": {"%sGGA_X_%s"%(mega, fmap[functional]): {"alpha": 1-hfx*0.01}},
            "x_hf": {"alpha": hfx*0.01},
            "c_functionals": {"%sGGA_C_%s"%(mega, fmap[functional]): {}}
        }
    else:
        raise ValueError("This functional has not been implemented with HFX resampling yet: ", functional)
    return hfx_func


def write_jobscript(psi4_config):
    if not "cluster" in psi4_config:
        mem = int(psi4_config['memory'].split(" ")[0])/1000
        with open("./jobscript.sh", "w") as fo:
            fo.write("#$ -S /bin/bash\n")
            fo.write("#$ -N psi4_dft\n")
            fo.write("#$ -R y\n")
            fo.write("#$ -cwd\n")
            fo.write("#$ -l h_rt=240:00:00\n")
            fo.write("#$ -l h_rss=%dG\n"% (mem))
            fo.write("#$ -q cpus\n")
            fo.write("#$ -l cpus=1\n")
            fo.write("#$ -pe smp %d\n"%psi4_config['num_threads'])
            fo.write("# -fin *.py\n")
            fo.write("# -fin b3lyp\n")
            fo.write("# -fin *.xyz\n")
            fo.write("# -fin *.molden\n")
            fo.write("# -fin *.json\n")
            fo.write("# -fin *\n")

            fo.write("source /home/crduan/.bashrc\n")
            fo.write("conda activate /home/crduan/miniconda/envs/mols_py36\n")
            fo.write("export PSI_SCRATCH='./'\n")
            fo.write("echo 'psi4 scr: ' $PSI_SCRATCH\n")
            fo.write("python -u loop_run.py  > $SGE_O_WORKDIR/nohup1.out\n")
            fo.write("python -u loop_run.py  > $SGE_O_WORKDIR/nohup2.out\n")
            fo.write("python -u loop_run.py  > $SGE_O_WORKDIR/nohup3.out\n")
            if "hfx_rescue" in psi4_config and psi4_config["hfx_rescue"]:
                fo.write("echo rescuing...\n")
                fo.write("python -u loop_rescue.py > $SGE_O_WORKDIR/rescue_nohup1.out\n")
                fo.write("python -u loop_rescue.py > $SGE_O_WORKDIR/rescue_nohup2.out\n")
                fo.write("python -u loop_rescue.py > $SGE_O_WORKDIR/rescue_nohup3.out\n")
            fo.write("mv * $SGE_O_WORKDIR\n")
    elif psi4_config["cluster"] == "expanse":
        mem = int(psi4_config['memory'].split(" ")[0])/psi4_config['num_threads']/1000
        with open("./jobscript.sh", "w") as fo:
            fo.write("#!/bin/sh\n")
            fo.write("#SBATCH -A mit136\n")
            fo.write("#SBATCH --job-name=psi4_multiDFA\n")
            fo.write("#SBATCH --partition=shared\n")
            fo.write("#SBATCH -t 48:00:00\n")
            fo.write("#SBATCH --nodes=1\n")
            fo.write("#SBATCH --ntasks-per-node=16\n")
            fo.write("#SBATCH --error=job.%J.err\n")
            fo.write("#SBATCH --output=job.%J.out\n")
            fo.write("#SBATCH --export=ALL\n")
            fo.write("#SBATCH --mem=64G\n")

            fo.write("source /home/crduan/.bashrc\n")
            fo.write("conda activate mols_psi4\n")
            fo.write("export PSI_SCRATCH='./'\n")
            fo.write("echo 'psi4 scr: ' $PSI_SCRATCH\n")
            fo.write("python -u loop_run.py  > nohup1.out\n")
            fo.write("python -u loop_run.py  > nohup2.out\n")
            fo.write("python -u loop_run.py  > nohup3.out\n")
            if "hfx_rescue" in psi4_config and psi4_config["hfx_rescue"]:
                fo.write("echo rescuing...\n")
                fo.write("python -u loop_rescue.py > rescue_nohup1.out\n")
                fo.write("python -u loop_rescue.py > rescue_nohup2.out\n")
                fo.write("python -u loop_rescue.py > rescue_nohup3.out\n")


def run_bash(cmd, basedir, rundir):
    os.chdir(rundir)
    infile = resource_filename(Requirement.parse("molSimplify"), "molSimplify/job_manager/psi4_utils/loop_run.py")
    shutil.copy(infile, "./")
    infile_rescue = resource_filename(Requirement.parse("molSimplify"), "molSimplify/job_manager/psi4_utils/loop_rescue.py")
    shutil.copy(infile_rescue, "./")
    print("Executing: ", cmd, "at: ", rundir)
    subprocess.call(cmd, shell=True)
    os.chdir(basedir)
