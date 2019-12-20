import subprocess
import copy
import os
import glob
import numpy as np
import shutil
from molSimplify.Classes.mol3D import *
from tools import *
from molSimplifyAD.utils.pymongo_tools import connect2db, query_lowestE_converged


def isCSD(job):
    iscsd = True
    for ii in range(6):
        if not job[ii].isupper():
            iscsd = False
            break
    return iscsd


def call_molsimplify(geodir, job, jobname):
    liglist = ",".join(job["ligstr"].split("_"))
    tmp_name = str(np.random.randint(10 ** 18))  # assign a temporary name so that the results are findable
    bash_command = " ".join(["molsimplify ", '-core ' + job["metal"],
                             '-lig ' + str(liglist), '-ligocc 1,1,1,1,1,1',
                             '-keepHs yes,yes,yes,yes,yes,yes',
                             '-ligloc ' + 'yes', '-calccharge yes',
                             '-spin ' + str(job["spin"]), '-oxstate ' + str(job["ox"]),
                             "-ffoption " + "b", ' -ff UFF',
                             "-name", tmp_name])
    print("call: ", bash_command)
    bash_command = bash_command.split()
    subprocess.call(bash_command)

    file_name = os.path.join(os.path.expanduser('~'), 'Runs', tmp_name)
    inner_folder_path = glob.glob(os.path.join(file_name, '*'))[0]
    xyz_path = glob.glob(os.path.join(inner_folder_path, '*.xyz'))[0]
    xyz_name = os.path.split(xyz_path)[-1]
    charge = False
    with open(inner_folder_path + '/terachem_input', "r") as fo:
        for line in fo:
            if "charge" == line[:len("charge")]:
                charge = int(line.split()[-1])
    if not charge:
        raise ValueError("No charge is extracted from terachem input: ", inner_folder_path + '/terachem_input')
    shutil.copyfile(xyz_path, geodir + '/' + jobname + '.xyz')
    print(xyz_path, geodir + '/' + jobname + '.xyz')
    return charge


def write_xyz_from_db(geodir, jobname, optgeo):
    natoms = len(optgeo.split('\n')) - 1
    with open(geodir + '/' + jobname + ".xyz", "w") as fo:
        fo.write("%d\n" % natoms)
        fo.write("====Geometry adopted from the database====\n")
        fo.write(optgeo)


def generate_fake_results_from_db(rundir, jobname, tmcdoc):
    scrdir = rundir + '/scr'
    if not os.path.isdir(scrdir):
        os.makedirs(scrdir)
    _ = write_xyz_from_db(scrdir, 'optim', tmcdoc["opt_geo"])
    if tmcdoc['wavefunction']:
        if int(tmcdoc['spin']) == 1:
            try:
                shutil.copy(tmcdoc['wavefunction']['c0'], scrdir + '/c0')
            except:
                pass
        else:
            try:
                shutil.copy(tmcdoc['wavefunction']['ca0'], scrdir + '/ca0')
                shutil.copy(tmcdoc['wavefunction']['cb0'], scrdir + '/cb0')
            except:
                pass
    inpath = rundir + '/' + jobname + '.in'
    with open(inpath, "w") as fo:
        fo.write("=======This outfile is FAKE and generated artificialy======\n")
        fo.write("charge %d\n" % int(tmcdoc["charge"]))
        fo.write("spinmult %d\n" % int(tmcdoc["spin"]))
        fo.write("method %s\n" % str(tmcdoc["functional"]))
        fo.write("basis %s\n" % str(tmcdoc["basis"]))
        fo.write("coordinates %s.xyz\n" % jobname)
    outpath = rundir + '/' + jobname + '.out'
    with open(outpath, "w") as fo:
        fo.write("=======This outfile is FAKE and generated artificialy======\n")
        fo.write('Startfile from command line: %s.in\n' % jobname)
        fo.write('Hartree-Fock exact exchange:          %.2f\n' % (float(tmcdoc['alpha']) * 0.01))
        fo.write('DFT Functional requested: %s\n' % str(tmcdoc['functional']))
        fo.write('*                    TeraChem %s            *\n' % str(tmcdoc['terachem_version']))
        fo.write('Alpha level shift: 0.25\n')
        fo.write('Beta level shift: 0.25\n')
        fo.write('Using basis set: %s\n' % str(tmcdoc["basis"]))
        fo.write("Total charge:    %d\n" % int(tmcdoc["charge"]))
        fo.write("Spin multiplicity: %d\n" % int(tmcdoc["spin"]))
        fo.write("SPIN S-SQUARED: %f (exact: %f)\n" % (float(tmcdoc["ss_act"]), float(tmcdoc["ss_target"])))
        fo.write("-=#=-      Optimization Cycle     1   -=#=-\n")
        fo.write("FINAL ENERGY: %.10f a.u.\n" % float(tmcdoc["energy"]))
        fo.write("-=#=-     Optimization Converged.     -=#=-\n")
        fo.write("Total processing time: 0.00 sec\n")
        fo.write("Job finished: A time of mystery\n")
    return outpath


def populate_single_job(basedir, job, db):
    geodir = basedir + "/initial_geometry/"
    if not os.path.isdir(geodir):
        os.makedirs(geodir)
    iscsd = isCSD(job['ligstr'])
    query_constraints = {"metal": job['metal'], "spin": job["spin"], "ligstr": job["ligstr"], "alpha": 20}
    if not iscsd:
        query_constraints.update({"ox": job["ox"]})
        jobname = "_".join([job['metal'], str(job['ox']), str(job['spin']), job['ligstr']])
    else:
        jobname = "_".join([job['ligstr'], job['metal'], str(job['spin'])])
    tmcdoc, recover = None, True
    if not db == None:
        tmcdoc = query_lowestE_converged(db, collection='oct', constraints=query_constraints)
        if not tmcdoc == None:
            print("Bingo! Optimized geometry found in db: ", query_constraints)
            try:
                charge = int(tmcdoc["charge"])
                spin = int(tmcdoc["spin"])
                energy = float(tmcdoc["energy"])
                wfn = tmcdoc['wavefunction']
                ss_act, ss_target = float(tmcdoc["ss_act"]), float(tmcdoc["ss_target"])
                write_xyz_from_db(geodir, jobname, tmcdoc["opt_geo"])
            except:
                recover = False
        else:
            print("generate initial geometry from molsimplify...")
            charge = call_molsimplify(geodir, job, jobname)
    else:
        if not iscsd:
            print("NO db connection! Generate initial geometry from molsimplify...")
            charge = call_molsimplify(basedir, job, jobname)
        else:
            raise ValueError("Cannot generate initial geometry for CSD complex...")
    rundir = basedir + '/' + jobname
    if not os.path.isdir(rundir) and recover:
        os.makedirs(rundir)
        shutil.copyfile(geodir + '/' + jobname + '.xyz', rundir + "/" + jobname + ".xyz")
        os.chdir(rundir)
        # Add fake files etc for a smooth carry-on in job manager for further dependent jobs.
        if not tmcdoc == None:
            outpath = generate_fake_results_from_db(rundir, jobname, tmcdoc)
        else:
            tools.write_input(jobname, charge, int(job["spin"]), run_type='minimize', solvent=False)
            tools.write_jobscript(jobname)
    elif os.path.isdir(rundir):
        print("folder exist.")
    else:
        print('WARNING: cannot recover %s' % jobname)

    os.chdir(basedir)


def populate_list_of_jobs(basedir, jobs, db_communicate=True):
    if db_communicate:
        db = connect2db(user="readonly_user", pwd="readonly", host="localhost",
                        port=27017, database="tmc", auth=True)
    else:
        db = None
    for job in jobs:
        populate_single_job(basedir, job, db)
