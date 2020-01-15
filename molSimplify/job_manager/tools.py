import os
import glob
import copy
import numpy as np
import subprocess
import pandas as pd
import shutil
import time
from molSimplify.job_manager.classes import resub_history, textfile
from ast import literal_eval


def ensure_dir(dirpath):
    if not os.path.isdir(dirpath):
        os.makedirs(dirpath)


def check_valid_outfile(path):
    # The nohup.out file gets caught in the find statement
    # use this function so that we only get TeraChem.outs
    endpath = os.path.split(path)[-1]
    if 'nohup.out' in endpath or endpath.startswith('.'):
        return False
    else:
        return True


def try_float(obj):
    # Converts an object to a floating point if possible
    try:
        floating_point = float(obj)
    except:
        floating_point = obj
    return floating_point


def invert_dictionary(dictionary):
    new_dict = dict()
    for key in list(dictionary.keys()):
        if type(dictionary[key]) == list:
            for entry in dictionary[key]:
                if entry in list(new_dict.keys()):
                    raise Exception('Dictionary inversion failed, values do not serve as unique keys')
                new_dict[entry] = key
        else:
            if dictionary[key] in list(new_dict.keys()):
                raise Exception('Dictionary inversion failed, values do not serve as unique keys')
            new_dict[dictionary[key]] = key
    return new_dict


def call_bash(string, error=False, version=1):
    if version == 1:
        p = subprocess.Popen(string.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    elif version == 2:
        p = subprocess.Popen(string, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    out, err = p.communicate()

    out = out.split('\n')
    if out[-1] == '':
        out = out[:-1]

    if error:
        return out, err
    else:
        return out


def convert_to_absolute_path(path):
    if path[0] != '/':
        path = os.path.join(os.getcwd(), path)

    return path


def list_active_jobs(ids=False, home_directory=False, parse_bundles=False):
    #  @return A list of active jobs for the current user. By job name

    if (ids and parse_bundles) or (parse_bundles and not home_directory):
        raise Exception('Incompatible options passed to list_active_jobs()')
    if home_directory == 'in place':
        home_directory = os.getcwd()

    job_report = textfile()
    try:
        job_report.lines = call_bash("qstat -r")
    except:
        job_report.lines = []

    names = job_report.wordgrab('jobname:', 2)[0]
    names = [i for i in names if i]  # filters out NoneTypes

    if ids:
        job_ids = []
        line_indices_of_jobnames = job_report.wordgrab('jobname:', 2, matching_index=True)[0]
        line_indices_of_jobnames = [i for i in line_indices_of_jobnames if i]  # filters out NoneTypes
        for line_index in line_indices_of_jobnames:
            job_ids.append(int(job_report.lines[line_index - 1].split()[0]))
        if len(names) != len(job_ids):
            print(len(names))
            print(len(job_ids))
            raise Exception('An error has occurred in listing active jobs!')
        return names, job_ids

    if parse_bundles and os.path.isfile(os.path.join(home_directory,'bundle','bundle_id')):

        fil = open(os.path.join(home_directory,'bundle','bundle_id'),'r')
        identifier = fil.readlines()[0]
        fil.close()  

        bundles = [i for i in names if i.startswith('bundle_')]
        bundles = [i.rsplit('_',1)[0] for i in names if i.endswith(identifier)]
        names = [i for i in names if i not in bundles]

        for bundle in bundles:
            info_path = glob.glob(os.path.join(home_directory, 'bundle', bundle, '*_info'))[0]
            fil = open(info_path, 'r')
            lines = fil.readlines()
            lines = [i[:-1] if i.endswith('\n') else i for i in lines]
            fil.close()
            names.extend(lines)

    return names

def get_number_active():
    active_jobs = list_active_jobs()
    def check_active(path, active_jobs=active_jobs):
    # Given a path, checks if it's in the queue currently:
        name = os.path.split(path)[-1]
        name = name.rsplit('.', 1)[0]
        if name in active_jobs:
            return True
        else:
            return False

    outfiles = find('*.out')
    outfiles = filter(check_valid_outfile,outfiles)

    active_non_bundles = [i for i in outfiles if check_active(i)]

    if os.path.isdir('bundle'):
        fil = open(os.path.join('bundle','bundle_id'))
        identifier = fil.readlines()[0]
        if identifier.endswith('\n'):
            identifier = identifier[:-1]
        fil.close()

        active_bundles = [i for i in active_jobs if i.startswith('bundle_')]
        active_bundles = [i for i in active_jobs if i.endswith(identifier)]
        active_jobs = active_bundles+active_non_bundles

        return len(active_jobs)
    else:
        return len(active_non_bundles)

def check_completeness(directory='in place', max_resub=5, configure_dict=False):
    ## Takes a directory, returns lists of finished, failed, and in-progress jobs
    outfiles = find('*.out', directory)
    outfiles = list(filter(check_valid_outfile, outfiles))

    results_tmp = [read_outfile(outfile, short_ouput=True) for outfile in outfiles]
    results_tmp = list(zip(outfiles, results_tmp))
    results_dict = dict()
    for outfile, tmp in results_tmp:
        results_dict[outfile] = tmp

    active_jobs = list_active_jobs(home_directory=directory,parse_bundles=True)

    def check_finished(path, results_dict=results_dict):
        # Return True if the outfile corresponds to a complete job, False otherwise
        if results_dict[path]['finished']:
            return True
        else:
            return False

    def check_active(path, active_jobs=active_jobs):
        # Given a path, checks if it's in the queue currently:
        name = os.path.split(path)[-1]
        name = name.rsplit('.', 1)[0]
        if name in active_jobs:
            return True
        else:
            return False

    def check_needs_resub(path):
        if os.path.isfile(path.rsplit('.', 1)[0] + '.pickle'):
            history = resub_history()
            history.read(path)
            return history.needs_resub
        else:
            return False

    def check_waiting(path):
        if os.path.isfile(path.rsplit('.', 1)[0] + '.pickle'):
            history = resub_history()
            history.read(path)
            if history.waiting:
                return True
        return False

    def grab_waiting(path):
        if os.path.isfile(path.rsplit('.', 1)[0] + '.pickle'):
            history = resub_history()
            history.read(path)
            return history.waiting
        raise Exception('Attempting to grab a "waiting" criteria that does not exist')

    def check_chronic_failure(path):
        if os.path.isfile(path.rsplit('.', 1)[0] + '.pickle'):
            history = resub_history()
            history.read(path)
            if history.resub_number >= max_resub:
                return True
        else:
            return False

    def check_spin_contaminated(path, results_dict=results_dict):
        results = results_dict[path]
        if configure_dict and "ss_cutoff" in configure_dict:
            ss_cutoff = configure_dict['ss_cutoff']
        else:
            ss_cutoff = 1.0
        if results['finished']:
            if type(results['s_squared_ideal']) == float:
                if abs(results['s_squared'] - results['s_squared_ideal']) > ss_cutoff:
                    return True
        return False

    def check_scf_error(path, results_dict=results_dict):
        results = results_dict[path]
        if results['scf_error']:
            return True
        else:
            return False

    def check_thermo_grad_error(path, results_dict=results_dict):
        results = results_dict[path]
        if results['thermo_grad_error']:
            return True
        else:
            return False

    active_jobs = list(filter(check_active, outfiles))
    finished = list(filter(check_finished, outfiles))
    needs_resub = list(filter(check_needs_resub, outfiles))
    waiting = list(filter(check_waiting, outfiles))
    spin_contaminated = list(filter(check_spin_contaminated, outfiles))
    all_scf_errors = list(filter(check_scf_error, outfiles))
    thermo_grad_errors = list(filter(check_thermo_grad_error, outfiles))
    chronic_errors = list(filter(check_chronic_failure, outfiles))
    errors = list(set(outfiles) - set(active_jobs) - set(finished))
    scf_errors = list(filter(check_scf_error, errors))

    # Look for additional active jobs that haven't yet generated outfiles
    jobscript_list = find('*_jobscript', directory)
    jobscript_list = [i.rsplit('_', 1)[0] + '.out' for i in jobscript_list]
    extra_active_jobs = list(filter(check_active, jobscript_list))
    active_jobs.extend(extra_active_jobs)

    # Sort out conflicts in order of reverse priority
    # A job only gets labelled as finished if it's in no other category
    # A job always gets labelled as active if it fits that criteria, even if it's in every other category too
    finished = list(set(finished) - set(needs_resub) - set(spin_contaminated) - set(errors) - set(scf_errors) - set(
        thermo_grad_errors) - set(waiting) - set(chronic_errors) - set(active_jobs))
    needs_resub = list(
        set(needs_resub) - set(spin_contaminated) - set(errors) - set(scf_errors) - set(thermo_grad_errors) - set(
            waiting) - set(chronic_errors) - set(active_jobs))
    spin_contaminated = list(
        set(spin_contaminated) - set(errors) - set(scf_errors) - set(thermo_grad_errors) - set(waiting) - set(
            chronic_errors) - set(active_jobs))
    errors = list(
        set(errors) - set(scf_errors) - set(thermo_grad_errors) - set(waiting) - set(chronic_errors) - set(active_jobs))
    scf_errors = list(set(scf_errors) - set(thermo_grad_errors) - set(waiting) - set(chronic_errors) - set(active_jobs))
    thermo_grad_errors = list(set(thermo_grad_errors) - set(waiting) - set(chronic_errors) - set(active_jobs))
    waiting = list(set(waiting) - set(chronic_errors) - set(active_jobs))
    chronic_errors = list(set(chronic_errors) - set(active_jobs))
    active_jobs = list(set(active_jobs))

    results = {'Finished': finished, 'Active': active_jobs, 'Error': errors, 'Resub': needs_resub,
               'Spin_contaminated': spin_contaminated, 'Chronic_error': chronic_errors,
               'Thermo_grad_error': thermo_grad_errors, 'Waiting': waiting, 'SCF_Error': scf_errors}

    # There are two special categories which operate a bit differently: waiting "SCF_Errors_Including_Active"
    # inverted_results = invert_dictionary(results)
    waiting = [{i: grab_waiting(i)} for i in waiting]
    results['Waiting'] = waiting
    results['SCF_Errors_Including_Active'] = all_scf_errors

    return results


def find(key, directory='in place'):
    ## Looks for all files with a matching key in their name within directory
    #  @return A list of paths
    if directory == 'in place':
        directory = os.getcwd()

    bash = 'find ' + directory + ' -name ' + key

    return call_bash(bash)


def qsub(jobscript_list):
    ## Takes a list of paths to jobscripts (like output by find) and submits them with qsub
    # Handles cases where a single jobscript is submitted instead of a list
    if type(jobscript_list) != list:
        jobscript_list = [jobscript_list]
    jobscript_list = [convert_to_absolute_path(i) for i in jobscript_list]

    stdouts = []
    for i in jobscript_list:
        home = os.getcwd()
        os.chdir(os.path.split(i)[0])
        jobscript = os.path.split(i)[1]
        stdout, stderr = call_bash('qsub ' + jobscript, error=True)
        stdouts.append(stdout)
        if len(stderr) > 0:
            raise Exception(stderr)
        os.chdir(home)
        time.sleep(1)

    return stdouts


def check_original(job):
    # Returns true if this job is an original job
    # Returns false if this job is already a derivative job

    name = os.path.split(job)[-1]
    name = name.rsplit('.', 1)[0]
    name = name.split('_')

    dependent_jobs = ['solvent', 'vertEA', 'vertIP', 'thermo', 'kp', 'rm', 'ultratight', 'HFXresampling',
                      'functional']
    if any(j in name for j in dependent_jobs):
        return False
    else:
        return True


def check_short_single_point(job):
    name = os.path.split(job)[-1]
    name = name.rsplit('.', 1)[0]
    name = name.split('_')

    short_jobs = ['solvent', 'kp', 'rm', 'functional', 'vertEA', 'vertIP']
    if any(j in name for j in short_jobs):
        return True
    else:
        return False


def extract_optimized_geo(PATH, custom_name=False):
    # Given the path to an optim.xyz file, this will extract optimized.xyz, which contains only the last frame
    # The file is written to the same directory as contained optim.xyz

    optim = open(PATH, 'r')
    lines = optim.readlines()
    optim.close()
    lines.reverse()
    if len(lines) == 0:
        lines = []
        print(('optim.xyz is empty for: ' + PATH))
    else:
        for i in range(len(lines)):
            if 'frame' in lines[i].split():
                break
        lines = lines[:i + 2]
        lines.reverse()
    homedir = os.path.split(PATH)[0]

    if custom_name:
        basedir = os.path.split(homedir)[0]
        xyz = glob.glob(os.path.join(basedir, '*.xyz'))[0]
        xyz = os.path.split(xyz)[-1]
        name = xyz[:-4] + '_optimized.xyz'
    else:
        name = 'optimized.xyz'

    optimized = open(os.path.join(homedir, name), 'w')
    for i in lines:
        optimized.write(i)
    optimized.close()

    return lines


def pull_optimized_geos(PATHs=[]):
    # Given a PATH or list of PATHs to optim.xyz files, these will grab the optimized geometries and move them into a local folder called optimized_geos
    if type(PATHs) != list:
        PATHs = [PATHs]
    if len(PATHs) == 0:
        PATHs = find('optim.xyz')  # Default behaviour will pull all optims in the current directory

    if os.path.isdir('opimized_geos'):
        raise Exception('optimized_geos already exists!')

    os.mkdir('optimized_geos')
    home = os.getcwd()
    for Path in PATHs:
        extract_optimized_geo(Path)
        scr_path = os.path.split(Path)[0]
        homedir = os.path.split(scr_path)[0]
        initial_xyz = glob.glob(os.path.join(homedir, '*.xyz'))
        if len(initial_xyz) != 1:
            print('Name could not be identified for: ' + Path)
            print('Naming with a random 6 digit number')
            name = str(np.random.randint(999999)) + '_optimized.xyz'
        else:
            initial_xyz = initial_xyz[0]
            name = os.path.split(initial_xyz)[-1]
            name = name.rsplit('.', 1)[0]
            name += '_optimized.xyz'

        shutil.move(os.path.join(os.path.split(Path)[0], 'optimized.xyz'), os.path.join(home, 'optimized_geos', name))


def read_outfile(outfile_path,short_ouput=False):
    ## Reads TeraChem and ORCA outfiles
    #  @param outfile_path complete path to the outfile to be read, as a string
    #  @return A dictionary with keys finalenergy,s_squared,s_squared_ideal,time
    
    output = textfile(outfile_path)
    output_type = output.wordgrab(['TeraChem','ORCA'],'whole_line')
    for counter,match in enumerate(output_type):
        if match[0]:
            break
        if counter == 1:
            if 'nohup' in outfile_path:
                print('Warning! Nohup file caught in outfile processing')
                print(outfile_path)
                counter = 0
            else:
                raise ValueError('.out file type not recognized for file: '+outfile_path)
    output_type = ['TeraChem','ORCA'][counter]
    if output_type == 'ORCA':
        raise Exception(outfile_path)
    
    name = None
    finished = False
    charge = None
    finalenergy = None
    min_energy = None
    s_squared = None
    s_squared_ideal = None
    scf_error = False
    time = None
    thermo_grad_error = False
    implicit_solvation_energy = None
    geo_opt_cycles = None
    thermo_vib = None
    thermo_vib_f = None
    thermo_suspect = None

    name = os.path.split(outfile_path)[-1]
    name = name.rsplit('.',1)[0]
    if output_type == 'TeraChem':
        
        charge = output.wordgrab(['charge:'],[2],first_line=True)[0]
        if charge:
            charge = int(charge)
        if not short_ouput:
            (finalenergy,s_squared,s_squared_ideal,time,thermo_grad_error,
             implicit_solvation_energy,geo_opt_cycles,
             thermo_vib,thermo_vib_f,thermo_suspect) = output.wordgrab(['FINAL','S-SQUARED:','S-SQUARED:','processing',
                                                                        'Maximum component of gradient is too large',
                                                                        'C-PCM contribution to final energy:',
                                                                        'Optimization Cycle','Thermal vibrational energy',
                                                                        'Thermal vibrational free energy',
                                                                        'Thermochemical Analysis is Suspect'],
                                                                        [2,2,4,3,0,4,3,7,10,0],last_line=True)
        if short_ouput:
            s_squared,s_squared_ideal,thermo_grad_error = output.wordgrab(['S-SQUARED:','S-SQUARED:','Maximum component of gradient is too large'],
                                                                           [2,4,0],last_line=True)

        if thermo_grad_error:
            thermo_grad_error = True
        else:
            thermo_grad_error = False
        if thermo_suspect:
            thermo_suspect = True
        else:
            thermo_suspect = False

        if s_squared_ideal:
            s_squared_ideal = float(s_squared_ideal.strip(')'))
        if implicit_solvation_energy:
            implicit_solvation_energy = try_float(implicit_solvation_energy.split(':')[-1])
            
        min_energy = output.wordgrab('FINAL',2,min_value = True)[0]
        
        is_finished = output.wordgrab(['finished:'],'whole_line',last_line=True)[0]
        if is_finished:
            if is_finished[0] == 'Job' and is_finished[1] == 'finished:':
                finished = True
        
        is_scf_error = output.wordgrab('DIIS',5,matching_index=True)[0]
        if is_scf_error[0]:
            is_scf_error = [output.lines[i].split() for i in is_scf_error]
        else:
            is_scf_error = []
        if type(is_scf_error) == list and len(is_scf_error) > 0:
            for scf in is_scf_error:
                if ('failed' in scf) and ('converge' in scf) and ('iterations,' in scf) and ('ADIIS' in scf):
                    scf = scf[5]
                    scf = int(scf.split('+')[0])
                    if scf > 5000:
                        scf_error = [True,scf]
                    
    if output_type == 'ORCA':
        finalenergy,s_squared,s_squared_ideal = output.wordgrab(['FINAL','<S**2>','S*(S+1)'],[8,13,12],last_line=True)
        timekey = 'TIME:'
        if type(output.wordgrab(timekey,'whole_line')) == list: 
            time = (float(output.wordgrab(timekey,3),last_line=True)*24*60*60
                   +float(output.wordgrab(timekey,5),last_line=True)*60*60
                   +float(output.wordgrab(timekey,7,last_line=True))*60
                   +float(output.wordgrab(timekey,9,last_line=True))
                   +float(output.wordgrab(timekey,11,last_line=True))*0.001)
        
    return_dict = {}
    return_dict['name'] = name
    return_dict['charge'] = try_float(charge)
    return_dict['finalenergy'] = try_float(finalenergy)
    return_dict['time'] = try_float(time)
    return_dict['s_squared'] = try_float(s_squared)
    return_dict['s_squared_ideal'] = try_float(s_squared_ideal)
    return_dict['finished'] = finished
    return_dict['min_energy'] = try_float(min_energy)
    return_dict['scf_error'] = scf_error
    return_dict['thermo_grad_error'] = thermo_grad_error
    return_dict['solvation_energy'] = implicit_solvation_energy
    return_dict['optimization_cycles'] = geo_opt_cycles
    return_dict['thermo_vib_energy'] = try_float(thermo_vib)
    return_dict['thermo_vib_free_energy'] = try_float(thermo_vib_f)
    return_dict['thermo_suspect'] = thermo_suspect
    return return_dict


def read_infile(outfile_path):
    # Takes the path to either the outfile or the infile of a job
    # Returns a dictionary of the job settings included in that infile

    root = outfile_path.rsplit('.', 1)[0]
    inp = textfile(root + '.in')
    charge, spinmult, solvent, run_type, levelshifta, levelshiftb, method, hfx, basis, dispersion, coordinates, guess = inp.wordgrab(
        ['charge ', 'spinmult ', 'epsilon ',
         'run ', 'levelshiftvala ',
         'levelshiftvalb ', 'method ',
         'HFX ', 'basis ', 'dispersion ',
         'coordinates ', 'guess '],
        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
        last_line=True)
    charge, spinmult = int(charge), int(spinmult)
    if guess:
        guess = True
    else:
        guess = False
    if method[0] == 'u':
        method = method[1:]

    convergence_thresholds = inp.wordgrab(
        ['min_converge_gmax ', 'min_converge_grms ', 'min_converge_dmax ', 'min_converge_drms ', 'min_converge_e ',
         'convthre '],
        [1] * 6, last_line=True)

    multibasis = inp.wordgrab(['$multibasis', '$end'], [0, 0], last_line=True, matching_index=True)
    if not multibasis[0]:
        multibasis = False
    else:
        multibasis = inp.lines[multibasis[0] + 1:multibasis[1]]

    constraints = inp.wordgrab(['$constraint_freeze', '$end'], [0, 0], last_line=True, matching_index=True)
    if not constraints[0]:
        constraints = False
    else:
        constraints = inp.lines[constraints[0] + 1:constraints[1]]

    if constraints and multibasis:
        raise Exception(
            'The current implementation of tools.read_infile() is known to behave poorly when an infile specifies both a multibasis and constraints')

    return_dict = {}

    for prop, prop_name in zip([charge, spinmult, solvent, run_type, levelshifta, levelshiftb, method, hfx,
                                basis, convergence_thresholds, multibasis, constraints, dispersion, coordinates, guess],
                               ['charge', 'spinmult', 'solvent', 'run_type', 'levelshifta', 'levelshiftb', 'method',
                                'hfx',
                                'basis', 'convergence_thresholds', 'multibasis', 'constraints', 'dispersion',
                                'coordinates', 'guess']):
        return_dict[prop_name] = prop
    return return_dict


# Read the global and local configure files to determine the derivative jobs requested and the settings for job recovery
# The global configure file should be in the same directory where resub() is called
# The local configure file should be in the same directory as the .out file
def read_configure(home_directory, outfile_path):
    def load_configure_file(directory):
        def strip_new_line(string):
            if string[-1] == '\n':
                return string[:-1]
            else:
                return string

        if directory == 'in place':
            directory = os.getcwd()

        configure = os.path.join(directory, 'configure')
        if os.path.isfile(configure):
            f = open(configure, 'r')
            configure = f.readlines()
            f.close()
            configure = list(map(strip_new_line, configure))
            return configure
        else:
            return []

    home_configure = load_configure_file(home_directory)
    if outfile_path:
        local_configure = load_configure_file(os.path.split(outfile_path)[0])
    else:
        local_configure = []

    # Determine which derivative jobs are requested
    solvent, vertEA, vertIP, thermo, dissociation, hfx_resample, functionalsSP = False, False, False, False, False, False, False
    for line in home_configure + local_configure:
        if 'solvent' in line or 'Solvent' in line:
            solvent = [float(p) for p in line.split()[1:]]
        if 'vertEA' in line or 'VertEA' in line:
            vertEA = True
        if 'vertIP' in line or 'VertIP' in line:
            vertIP = True
        if 'functionalsSP' in line or 'FunctionalsSP' in line:
            functionalsSP = [str(p) for p in line.split()[1:]]
        if 'thermo' in line or 'Thermo' in line:
            thermo = True
        if 'dissociation' in line or 'Dissociation' in line:
            dissociation = True
        if 'hfx_resample' in line or 'HFX_resample' in line:
            hfx_resample = True

    # Determine global settings for this run
    max_jobs, max_resub, levela, levelb, method, hfx, geo_check, sleep, job_recovery, dispersion = False, False, False, False, False, False, False, False, [], False
    ss_cutoff = False
    for configure in [home_configure, local_configure]:
        for line in home_configure:
            if 'max_jobs' in line.split(':'):
                max_jobs = int(line.split(':')[-1])
            if 'max_resub' in line.split(':'):
                max_resub = int(line.split(':')[-1])
            if 'levela' in line.split(':'):
                levela = float(line.split(':')[-1])
            if 'levelb' in line.split(':'):
                levelb = float(line.split(':')[-1])
            if 'method' in line.split(':'):
                method = line.split(':')[-1]
            if 'hfx' in line.split(':'):
                hfx = float(line.split(':')[-1])
            if 'geo_check' in line.split(':'):
                geo_check = line.split(':')[-1]
            if 'sleep' in line.split(':'):
                sleep = int(line.split(':')[-1])
            if 'job_recovery' in line.split(':'):
                job_recovery = line.split(':')[-1]
                # convert the string form of a python list to an actual list
                job_recovery = job_recovery[1:-1]
                job_recovery = job_recovery.split(',')
            if 'dispersion' in line.split(':'):
                dispersion = line.split(':')[-1]
            if 'ss_cutoff' in line.split(':'):
                ss_cutoff = float(line.split(':')[-1])
        # If global settings not specified, choose defaults:
        if not max_jobs:
            max_jobs = 50
        if not max_resub:
            max_resub = 5
        if not levela:
            levela = 0.25
        if not levelb:
            levelb = 0.25
        if not method:
            method = 'b3lyp'
        if not hfx:
            hfx = 0.20
        if not sleep:
            sleep = 7200
        if not ss_cutoff:
            ss_cutoff = 1.0
        # Octahedral defaults to True in original variable initiation

    return {'solvent': solvent, 'vertEA': vertEA, 'vertIP': vertIP, 'thermo': thermo, 'dissociation': dissociation,
            'hfx_resample': hfx_resample, 'max_jobs': max_jobs, 'max_resub': max_resub, 'levela': levela,
            'levelb': levelb, 'method': method, 'hfx': hfx, 'geo_check': geo_check, 'sleep': sleep,
            'job_recovery': job_recovery, 'dispersion': dispersion, 'functionalsSP': functionalsSP,
            'ss_cutoff': ss_cutoff}


def read_charges(PATH):
    # Takes the path to either the outfile or the charge_mull.xls and returns the charges
    PATH = convert_to_absolute_path(PATH)
    if len(PATH.rsplit('.', 1)) > 1:
        if PATH.rsplit('.', 1)[1] == 'out':
            PATH = os.path.join(os.path.split(PATH)[0], 'scr', 'charge_mull.xls')
    try:
        charge_mull = textfile(PATH)
        split_lines = [i.split() for i in charge_mull.lines]
        charges = [i[1] + ' ' + i[2] for i in split_lines]
        return charges
    except:
        return []


def read_mullpop(PATH):
    # Takes the path to either the outfile or the mullpop and returns the mullikan populations
    PATH = convert_to_absolute_path(PATH)
    if len(PATH.rsplit('.', 1)) > 1:
        if PATH.rsplit('.', 1)[1] == 'out':
            PATH = os.path.join(os.path.split(PATH)[0], 'scr', 'mullpop')

    mullpop = textfile(PATH)
    ### If multiple frames in mullpop, grab last frame
    total_lines = mullpop.wordgrab(['------------ ---------- ----------'],[1],matching_index=True)[0]
    if len(total_lines) > 1:
        mullpop.lines = mullpop.lines[total_lines[-2]+2:]

    split_lines = [i.split() for i in mullpop.lines]
    if len(split_lines[2]) == 6:
        pops = [i[1] + ' ' + i[5] for i in split_lines[1:-2]]
    else:
        pops = [i[1] + ' ' + i[5] + ' ' + i[9] for i in split_lines[2:-2]]

    return pops


def write_input(input_dictionary=dict(), name=None, charge=None, spinmult=None,
                run_type='energy', method='b3lyp', solvent=False,
                guess=False, custom_line=None, levelshifta=0.25, levelshiftb=0.25,
                convergence_thresholds=None, basis='lacvps_ecp', hfx=None, constraints=None,
                multibasis=False, coordinates=False, dispersion=False):
    # Writes a generic terachem input file
    # The neccessary parameters can be supplied as arguements or as a dictionary. If supplied as both, the dictionary takes priority
    # "Custom line" can be used to add additional lines to the infile and is not treated by an input dictionary
    # Note that the infile dictionary can have an additional key, "name", which is not poulated by read_infile()
    # If name is specified, the coordinates are generated based on the name, rather than based on the coordinates variable

    infile = dict()
    # If the input_dictionary exists,parse it and set the parameters, overwritting other specifications
    for prop, prop_name in zip([charge, spinmult, solvent, run_type, levelshifta, levelshiftb, method, hfx,
                                basis, convergence_thresholds, multibasis, constraints, dispersion, coordinates,
                                guess, custom_line, name],
                               ['charge', 'spinmult', 'solvent', 'run_type', 'levelshifta', 'levelshiftb', 'method',
                                'hfx',
                                'basis', 'convergence_thresholds', 'multibasis', 'constraints', 'dispersion',
                                'coordinates',
                                'guess', 'custom_line', 'name']):
        if prop_name in list(input_dictionary.keys()):
            infile[prop_name] = input_dictionary[prop_name]
        else:
            infile[prop_name] = prop

    if (not infile['charge'] and infile['charge'] != 0) or (not infile['spinmult'] and infile['spinmult'] != 0) or (
            not infile['name'] and not infile['coordinates']):
        print(('Name: ' + infile['name']))
        print(('Charge: ' + str(infile['charge'])))
        print(('Spinmult: ' + str(infile['spinmult'])))
        raise Exception('Minimum parameters not specified for writing infile')
    if type(infile['charge']) != int or type(infile['spinmult']) != int:
        print(('Charge Type: ' + str(type(infile['charge']))))
        print(('Spinmult Type: ' + str(type(infile['spinmult']))))
        raise Exception('Spin and Charge should both be integers!')

    if infile['spinmult'] != 1:
        infile['method'] = 'u' + infile['method']

    if infile['name']:
        infile['coordinates'] = infile['name'] + '.xyz'
    if not infile['name']:
        infile['name'] = os.path.split(infile['coordinates'])[-1].rsplit('.', 1)[0]

    input_file = open(infile['name'] + '.in', 'w')
    text = ['levelshiftvalb ' + str(infile['levelshiftb']) + '\n',
            'levelshiftvala ' + str(infile['levelshifta']) + '\n',
            'nbo yes\n',
            'run ' + infile['run_type'] + '\n',
            'scf diis+a\n',
            'coordinates ' + infile['coordinates'] + '\n',
            'levelshift yes\n',
            'gpus 1\n',
            'spinmult ' + str(infile['spinmult']) + '\n',
            'scrdir ./scr\n',
            'basis ' + infile['basis'] + '\n',
            'timings yes\n',
            'charge ' + str(infile['charge']) + '\n',
            'method ' + infile['method'] + '\n',
            'new_minimizer yes\n',
            'ml_prop yes\n',
            'poptype mulliken\n',
            'bond_order_list yes\n',
            'end']

    if infile['custom_line']:
        text = text[:15] + [infile['custom_line'] + '\n'] + text[15:]
    if infile['guess']:
        if type(infile['guess']) == bool:
            if infile['spinmult'] == 1:
                infile['guess'] = 'c0'
            else:
                infile['guess'] = 'ca0 cb0'
        text = text[:-1] + ['guess ' + infile['guess'] + '\n',
                            'end']
    if infile['run_type'] != 'ts':
        text = text[:-1] + ['maxit 500\n',
                            'end']

    if type(infile['convergence_thresholds']) == list:
        if infile['convergence_thresholds'][0]:
            thresholds = [line if line.endswith('\n') else line + '\n' for line in infile['convergence_thresholds']]
            tight_thresholds = "min_converge_gmax " + thresholds[0] + "min_converge_grms " + thresholds[
                1] + "min_converge_dmax " + thresholds[2] + "min_converge_drms " + thresholds[3] + "min_converge_e " + \
                               thresholds[4] + "convthre " + thresholds[5]
            text = text[:-1] + ['\n', tight_thresholds, 'end']

    if infile['dispersion']:
        text = text[:-1] + ['dispersion ' + infile['dispersion'] + '\n', 'end']

    if infile['multibasis']:
        multibasis = [line if line.endswith('\n') else line + '\n' for line in infile['multibasis']]
        text = text[:-1] + ['\n', '$multibasis\n'] + multibasis + ['$end\n', 'end']

    if infile['constraints']:
        constraints = [line if line.endswith('\n') else line + '\n' for line in infile['constraints']]
        text = text[:-1] + ['\n', '$constraint_freeze\n'] + constraints + ['$end\n', 'end']

    if type(infile['hfx']) == int or type(infile['hfx']) == float:
        text = text[:-1] + ['\n',
                            'HFX ' + str(infile['hfx']) + '\n',
                            'end']

    if infile['solvent']:  # Adds solvent correction, if requested
        text = text[:-1] + ['\n',
                            'pcm cosmo\n',
                            'epsilon %.2f\n' % float(infile['solvent']),
                            'pcm_radii read\n',
                            'pcm_radii_file /home/harperd/pcm_radii\n',
                            'end']
    for lines in text:
        input_file.write(lines)
    input_file.close()


def write_jobscript(name, custom_line=None, time_limit='96:00:00', terachem_line=True):
    # Writes a generic terachem jobscript
    # custom line allows the addition of extra lines just before the export statement

    jobscript = open(name + '_jobscript', 'w')
    text = ['#$ -S /bin/bash\n',
            '#$ -N ' + name + '\n',
            '#$ -cwd\n',
            '#$ -R y\n',
            '#$ -l h_rt=' + time_limit + '\n',
            '#$ -l h_rss=8G\n',
            '#$ -q gpus|gpusnew|gpusnewer\n',
            '#$ -l gpus=1\n',
            '#$ -pe smp 1\n',
            '# -fin ' + name + '.in\n',
            '# -fin ' + name + '.xyz\n',
            '# -fout scr/\n',
            'source /etc/profile.d/modules.sh\n',
            'module load terachem/tip\n',
            'export OMP_NUM_THREADS=1\n',
            'terachem ' + name + '.in ' + '> $SGE_O_WORKDIR/' + name + '.out\n']
    if not terachem_line:
        text = text[:-1]

    if custom_line:
        if type(custom_line) == list:
            text = text[:12] + custom_line + text[12:]
        else:
            text = text[:12] + [custom_line + '\n'] + text[12:]

    for i in text:
        jobscript.write(i)
    jobscript.close()


def bundle_jobscripts(home_directory,jobscript_paths,max_bundle_size = 10):

    number_of_bundles = int(float(len(jobscript_paths))/float(max_bundle_size))

    bundles, i, many_bundles = [], 0, False
    for i in range(number_of_bundles):
        bundles.append(jobscript_paths[max_bundle_size * i:max_bundle_size * (i + 1)])
        many_bundles = True
    if many_bundles:
        total_length = [len(ii) for ii in bundles]
        total_length = np.sum(np.array(total_length))
        if total_length != len(jobscript_paths): #Triggers when the number of jobscripts is not divisable by the bundle size
            bundles.append(jobscript_paths[max_bundle_size * (i + 1):])
    else: #Triggers when there are fewere jobscripts than the bundle size
        bundles = [jobscript_paths]
    print('Bundling '+str(len(jobscript_paths))+' short jobs into '+str(len(bundles))+' jobscript(s)')

    output_jobscripts = []
    for bundle in bundles:
        output_jobscripts.append(sub_bundle_jobscripts(home_directory, bundle))

    return output_jobscripts



def sub_bundle_jobscripts(home_directory,jobscript_paths):
    #Takes a list of jobscript paths, and bundles them into a single jobscript
    #Records information about which jobs were bundled together in the run's home directory
    if not os.path.isdir(os.path.join(home_directory,'bundle')):
        os.mkdir(os.path.join(home_directory,'bundle'))
        fil = open(os.path.join(home_directory,'bundle','bundle_id'),'w')
        fil.write(str(np.random.randint(100000000000)))
        fil.close()
    fil = open(os.path.join(home_directory,'bundle','bundle_id'),'r')
    identifier = fil.readlines()[0]
    fil.close()  


    jobscript_paths = [convert_to_absolute_path(i) for i in jobscript_paths]

    existing_bundles = glob.glob(os.path.join(home_directory, 'bundle', '*'))
    existing_bundles = [i for i in existing_bundles if os.path.isdir(i)]
    if len(existing_bundles) > 0:
        existing_bundle_numbers = [int(os.path.split(i)[-1].split('_')[-1]) for i in existing_bundles]
    else:
        existing_bundle_numbers = [0]

    # Create a directory for this jobscript bundle
    os.mkdir(os.path.join(home_directory, 'bundle', 'bundle_' + str(max(existing_bundle_numbers) + 1)))

    # Record info about how the jobs are being bundled
    fil = open(os.path.join(home_directory, 'bundle', 'bundle_' + str(max(existing_bundle_numbers) + 1),
                            'bundle_' + str(max(existing_bundle_numbers) + 1) + '_info'), 'w')
    for i in jobscript_paths[:-1]:
        fil.write(os.path.split(i)[-1].rsplit('_', 1)[0] + '\n')
    fil.write(os.path.split(jobscript_paths[-1])[-1].rsplit('_', 1)[0])
    fil.close()

    # Write a jobscript for the job bundle
    home = os.getcwd()
    os.chdir(os.path.join(home_directory,'bundle','bundle_'+str(max(existing_bundle_numbers)+1)))
    write_jobscript(str('bundle_'+str(max(existing_bundle_numbers)+1))+'_'+identifier,terachem_line=False,time_limit='12:00:00')
    shutil.move('bundle_'+str(max(existing_bundle_numbers)+1)+'_'+identifier+'_jobscript','bundle_'+str(max(existing_bundle_numbers)+1))
    fil = open('bundle_'+str(max(existing_bundle_numbers)+1),'a')
    for i in jobscript_paths:
        infile = i.rsplit('_', 1)[0] + '.in'
        outfile = i.rsplit('_', 1)[0] + '.out'
        directory = os.path.split(i)[0]
        text = 'cd ' + directory + '\n' + 'terachem ' + infile + ' > ' + outfile
        fil.write(text + '\n')
    fil.close()
    os.chdir(home)

    return os.path.join(home_directory, 'bundle', 'bundle_' + str(max(existing_bundle_numbers) + 1),
                        'bundle_' + str(max(existing_bundle_numbers) + 1))


def prep_vertical_ip(path):
    # Given a path to the outfile of a finished run, this preps the files for a corresponding vertical IP run
    # Returns a list of the PATH(s) to the jobscript(s) to start the vertical IP calculations(s)
    home = os.getcwd()
    path = convert_to_absolute_path(path)

    results = read_outfile(path)
    if not results['finished']:
        raise Exception('This calculation does not appear to be complete! Aborting...')

    infile_dict = read_infile(path)

    if infile_dict['spinmult'] == 1:
        new_spin = [2]
    else:
        new_spin = [infile_dict['spinmult'] - 1, infile_dict['spinmult'] + 1]

    base = os.path.split(path)[0]

    optimxyz = os.path.join(base, 'scr', 'optim.xyz')
    extract_optimized_geo(optimxyz)

    ipname = results['name'] + '_vertIP'
    vertip_base_path = os.path.join(base, ipname)
    if os.path.isdir(vertip_base_path):
        return ['Directory for vertIP single point already exists']
    os.mkdir(vertip_base_path)
    os.chdir(vertip_base_path)

    jobscripts = []
    for calc in new_spin:
        if calc < 7:
            name = results['name'] + '_vertIP_' + str(calc)
            PATH = os.path.join(vertip_base_path, str(calc))
            if os.path.isdir(PATH):
                jobscripts.append('File for vert IP spin ' + str(calc) + 'already exists')
            else:
                os.mkdir(PATH)
                os.chdir(PATH)

                shutil.copyfile(os.path.join(base, 'scr', 'optimized.xyz'), os.path.join(PATH, name + '.xyz'))

                local_infile_dict = copy.copy(infile_dict)
                local_infile_dict['charge'], local_infile_dict['guess'] = infile_dict['charge'] + 1, False
                local_infile_dict['run_type'], local_infile_dict['spinmult'] = 'energy', calc
                local_infile_dict['name'] = name
                local_infile_dict['levelshifta'], local_infile_dict['levelshiftb'] = 0.25, 0.25
                write_input(local_infile_dict)
                write_jobscript(name)

                jobscripts.append(os.path.join(PATH, name + '_jobscript'))

    os.chdir(home)

    return jobscripts


def prep_vertical_ea(path):
    # Given a path to the outfile of a finished run, this preps the files for a corresponding vertical EA run
    # Returns a list of the PATH(s) to the jobscript(s) to start the vertical IP calculations(s)
    home = os.getcwd()
    path = convert_to_absolute_path(path)

    results = read_outfile(path)
    if not results['finished']:
        raise Exception('This calculation does not appear to be complete! Aborting...')

    infile_dict = read_infile(path)

    if infile_dict['spinmult'] == 1:
        new_spin = [2]
    else:
        new_spin = [infile_dict['spinmult'] - 1, infile_dict['spinmult'] + 1]

    base = os.path.split(path)[0]

    optimxyz = os.path.join(base, 'scr', 'optim.xyz')
    extract_optimized_geo(optimxyz)

    eaname = results['name'] + '_vertEA'
    vertea_base_path = os.path.join(base, eaname)
    if os.path.isdir(vertea_base_path):
        return ['Directory for vertEA single point already exists']
    os.mkdir(vertea_base_path)
    os.chdir(vertea_base_path)

    jobscripts = []
    for calc in new_spin:
        if calc < 7:
            name = results['name'] + '_vertEA_' + str(calc)
            PATH = os.path.join(vertea_base_path, str(calc))
            if os.path.isdir(PATH):
                jobscripts.append('File for vert EA spin ' + str(calc) + 'already exists')
            else:
                os.mkdir(PATH)
                os.chdir(PATH)
                shutil.copyfile(os.path.join(base, 'scr', 'optimized.xyz'), os.path.join(PATH, name + '.xyz'))

                local_infile_dict = copy.copy(infile_dict)
                local_infile_dict['charge'], local_infile_dict['guess'] = infile_dict['charge'] - 1, False
                local_infile_dict['run_type'], local_infile_dict['spinmult'] = 'energy', calc
                local_infile_dict['name'] = name
                local_infile_dict['levelshifta'], local_infile_dict['levelshiftb'] = 0.25, 0.25

                write_input(local_infile_dict)
                write_jobscript(name)
                jobscripts.append(os.path.join(PATH, name + '_jobscript'))
    os.chdir(home)
    return jobscripts


def prep_solvent_sp(path, solvents=[78.9]):
    # Given a path to the outfile of a finished run, this preps the files for a single point solvent run
    # Uses the wavefunction from the gas phase calculation as an initial guess
    # Returns a list of the PATH(s) to the jobscript(s) to start the solvent sp calculations(s)
    home = os.getcwd()
    path = convert_to_absolute_path(path)

    results = read_outfile(path)
    if not results['finished']:
        raise Exception('This calculation does not appear to be complete! Aborting...')

    infile_dict = read_infile(path)

    base = os.path.split(path)[0]

    optimxyz = os.path.join(base, 'scr', 'optim.xyz')
    extract_optimized_geo(optimxyz)

    # Now, start generating the new directory
    solname = results['name'] + '_solvent'
    solvent_base_path = os.path.join(base, solname)
    if os.path.isdir(solvent_base_path):
        return ['Directory for solvent single point already exists']
    os.mkdir(solvent_base_path)
    os.chdir(solvent_base_path)

    jobscripts = []
    for sol_val in solvents:
        PATH = os.path.join(solvent_base_path, str(sol_val).replace('.', '_'))
        ensure_dir(PATH)
        name = results['name'] + "_solvent_" + str(sol_val)
        shutil.copyfile(os.path.join(base, 'scr', 'optimized.xyz'), os.path.join(PATH, name + '.xyz'))
        guess = False
        os.chdir(PATH)
        if infile_dict['spinmult'] == 1:
            if os.path.isfile(os.path.join(base, 'scr', 'c0')):
                shutil.copyfile(os.path.join(base, 'scr', 'c0'), os.path.join(PATH, 'c0'))
                write_jobscript(name, custom_line='# -fin c0')
                guess = True
        else:
            if os.path.isfile(os.path.join(base, 'scr', 'ca0')) and os.path.isfile(os.path.join(base, 'scr', 'cb0')):
                shutil.copyfile(os.path.join(base, 'scr', 'ca0'), os.path.join(PATH, 'ca0'))
                shutil.copyfile(os.path.join(base, 'scr', 'cb0'), os.path.join(PATH, 'cb0'))
                write_jobscript(name, custom_line=['# -fin ca0\n', '# -fin cb0\n'])
                guess = True

        local_infile_dict = copy.copy(infile_dict)
        local_infile_dict['solvent'], local_infile_dict['guess'] = sol_val, guess
        local_infile_dict['run_type'] = 'energy'
        local_infile_dict['name'] = name
        local_infile_dict['levelshifta'], local_infile_dict['levelshiftb'] = 0.25, 0.25

        write_input(local_infile_dict)
        os.chdir(home)
        jobscripts.append(os.path.join(PATH, name + '_jobscript'))
    return jobscripts


def prep_functionals_sp(path, functionalsSP):
    home = os.getcwd()
    path = convert_to_absolute_path(path)
    results = read_outfile(path)
    if not results['finished']:
        raise Exception('This calculation does not appear to be complete! Aborting...')

    infile_dict = read_infile(path)
    base = os.path.split(path)[0]
    optimxyz = os.path.join(base, 'scr', 'optim.xyz')
    extract_optimized_geo(optimxyz)

    # Now, start generating the new directory
    funcname = results['name'] + '_functionalsSP'
    functional_base_path = os.path.join(base, funcname)
    if os.path.isdir(functional_base_path):
        return ['Directory for functional single point already exists']
    os.mkdir(functional_base_path)
    os.chdir(functional_base_path)

    jobscripts = []
    for func in functionalsSP:
        PATH = os.path.join(functional_base_path, str(func))
        ensure_dir(PATH)
        name = results['name'] + "_functional_" + str(func)
        shutil.copyfile(os.path.join(base, 'scr', 'optimized.xyz'), os.path.join(PATH, name + '.xyz'))
        guess = False
        os.chdir(PATH)
        if infile_dict['spinmult'] == 1:
            if os.path.isfile(os.path.join(base, 'scr', 'c0')):
                shutil.copyfile(os.path.join(base, 'scr', 'c0'), os.path.join(PATH, 'c0'))
                write_jobscript(name, custom_line='# -fin c0')
                guess = True
        else:
            if os.path.isfile(os.path.join(base, 'scr', 'ca0')) and os.path.isfile(os.path.join(base, 'scr', 'cb0')):
                shutil.copyfile(os.path.join(base, 'scr', 'ca0'), os.path.join(PATH, 'ca0'))
                shutil.copyfile(os.path.join(base, 'scr', 'cb0'), os.path.join(PATH, 'cb0'))
                write_jobscript(name, custom_line=['# -fin ca0\n', '# -fin cb0\n'])
                guess = True

        local_infile_dict = copy.copy(infile_dict)
        local_infile_dict['solvent'], local_infile_dict['guess'] = False, guess
        local_infile_dict['run_type'] = 'energy'
        local_infile_dict['name'] = name
        local_infile_dict['levelshifta'], local_infile_dict['levelshiftb'] = 0.25, 0.25
        local_infile_dict['method'] = func

        write_input(local_infile_dict)
        os.chdir(home)
        jobscripts.append(os.path.join(PATH, name + '_jobscript'))
    return jobscripts


def prep_thermo(path):
    # Given a path to the outfile of a finished run, this preps the files for a thermo calculation
    # Uses the wavefunction from the previous calculation as an initial guess
    # Returns a list of the PATH(s) to the jobscript(s) to start the solvent sp calculations(s)
    home = os.getcwd()
    path = convert_to_absolute_path(path)

    results = read_outfile(path)
    infile_dict = read_infile(path)

    base = os.path.split(path)[0]

    optimxyz = os.path.join(base, 'scr', 'optim.xyz')
    extract_optimized_geo(optimxyz)

    # Now, start generating the new directory
    name = results['name'] + '_thermo'
    PATH = os.path.join(base, name)
    if os.path.isdir(PATH):
        return ['Thermo Calculation Directory already exists']

    os.mkdir(PATH)
    os.chdir(PATH)

    shutil.copyfile(os.path.join(base, 'scr', 'optimized.xyz'), os.path.join(PATH, name + '.xyz'))
    if infile_dict['spinmult'] == 1:
        shutil.copyfile(os.path.join(base, 'scr', 'c0'), os.path.join(PATH, 'c0'))
        write_jobscript(name, custom_line='# -fin c0')
    if infile_dict['spinmult'] != 1:
        shutil.copyfile(os.path.join(base, 'scr', 'ca0'), os.path.join(PATH, 'ca0'))
        shutil.copyfile(os.path.join(base, 'scr', 'cb0'), os.path.join(PATH, 'cb0'))
        write_jobscript(name, custom_line=['# -fin ca0\n', '# -fin cb0\n'])

    local_infile_dict = copy.copy(infile_dict)
    local_infile_dict['guess'] = True
    local_infile_dict['run_type'] = 'frequencies'
    local_infile_dict['name'] = name

    write_input(local_infile_dict)

    os.chdir(home)

    return [os.path.join(PATH, name + '_jobscript')]


def prep_ultratight(path):
    # Given a path to the outfile of a finished run, this preps a run with tighter convergence criteria
    # Uses the wavefunction and geometry from the previous calculation as an initial guess
    # Returns a list of the PATH(s) to the jobscript(s) to start the solvent sp calculations(s)
    home = os.getcwd()
    path = convert_to_absolute_path(path)

    results = read_outfile(path)
    if not results['finished']:
        raise Exception('This calculation does not appear to be complete! Aborting...')

    infile_dict = read_infile(path)

    base = os.path.split(path)[0]

    optimxyz = os.path.join(base, 'scr', 'optim.xyz')
    extract_optimized_geo(optimxyz)

    # Now, start generating the new directory
    name = results['name'] + '_ultratight'
    PATH = os.path.join(base, name)

    if not os.path.isdir(PATH):  # First time that ultratight has been run, create necessary files
        os.mkdir(PATH)
        os.chdir(PATH)

        if os.path.exists(name + '.in') or os.path.exists(name + '.out') or os.path.exists(name + '_jobscript'):
            raise Exception('This tightened convergence run appears to already exist. Aborting...')

        shutil.copyfile(os.path.join(base, 'scr', 'optimized.xyz'), os.path.join(PATH, name + '.xyz'))
        if infile_dict['spinmult'] == 1:
            shutil.copyfile(os.path.join(base, 'scr', 'c0'), os.path.join(PATH, 'c0'))
            write_jobscript(name, custom_line='# -fin c0')
        elif infile_dict['spinmult'] != 1:
            shutil.copyfile(os.path.join(base, 'scr', 'ca0'), os.path.join(PATH, 'ca0'))
            shutil.copyfile(os.path.join(base, 'scr', 'cb0'), os.path.join(PATH, 'cb0'))
            write_jobscript(name, custom_line=['# -fin ca0\n', '# -fin cb0\n'])

        criteria = ['2.25e-04', '1.5e-04', '0.9e-03', '0.6e-03', '0.5e-06', '1.5e-05']

        local_infile_dict = copy.copy(infile_dict)
        local_infile_dict['guess'] = True
        local_infile_dict['convergence_thresholds'] = criteria
        local_infile_dict['name'] = name

        write_input(local_infile_dict)

        # Make an empty .out file to prevent the resubmission module from mistakenly submitting this job twice
        f = open(name + '.out', 'w')
        f.close()

        os.chdir(home)

        return [os.path.join(PATH, name + '_jobscript')]

    else:  # This has been run before, further tighten the convergence criteria
        os.chdir(PATH)
        infile_dict = read_infile(os.path.join(PATH, name + '.out'))
        criteria = [str(float(i) / 2.) for i in infile_dict['convergence_thresholds']]

        local_infile_dict = copy.copy(infile_dict)
        local_infile_dict['guess'] = True
        local_infile_dict['convergence_thresholds'] = criteria
        local_infile_dict['name'] = name
        write_input(local_infile_dict)

        extract_optimized_geo(os.path.join(PATH, 'scr', 'optim.xyz'))
        shutil.copy(os.path.join(PATH, 'scr', 'optimized.xyz'), os.path.join(PATH, name + '.xyz'))

        os.chdir(home)

        return [os.path.join(PATH, name + '_jobscript')]


def prep_hfx_resample(path, hfx_values=[0, 5, 10, 15, 20, 25, 30]):
    # Given a path to the outfile of a finished run, this preps the files for hfx resampling
    # Uses the wavefunction from the gas phase calculation as an initial guess
    # Returns a list of the PATH(s) to the jobscript(s) to start the resampling calculations(s)
    home = os.getcwd()
    path = convert_to_absolute_path(path)
    base = os.path.split(path)[0]

    results = read_outfile(path)
    if not results['finished']:
        raise Exception('This calculation does not appear to be complete! Aborting...')

    # Check the state of the calculation and ensure than hfx resampling is valid
    infile_dict = read_infile(path)
    if infile_dict['method'] != 'b3lyp':
        raise Exception('HFX resampling may not behave well for methods other than b3lyp!')
    if not infile_dict['hfx']:
        infile_dict['hfx'] = 20
    if infile_dict['hfx'] not in hfx_values:
        raise Exception('HFX resampling list does not contain the original hfx value!')

    # Now, start generating the base directory to hold all the hfx resampling values
    name = results['name'] + '_HFXresampling'
    hfx_path = os.path.join(base, name)
    if not os.path.isdir(hfx_path):
        os.mkdir(hfx_path)
    os.chdir(hfx_path)

    # Make the directory for the original calculation
    subname = name + '_' + str(infile_dict['hfx'])
    PATH = os.path.join(hfx_path, subname)
    if not os.path.isdir(PATH):
        os.mkdir(PATH)
    os.chdir(PATH)

    if not os.path.exists(os.path.join(PATH, subname + '.out')):
        shutil.copyfile(path, subname + '.out')
        shutil.copyfile(path.rsplit('.', 1)[0] + '_jobscript', subname + '_jobscript')
        shutil.copyfile(path.rsplit('.', 1)[0] + '.in', subname + '.in')
        shutil.copytree(os.path.join(os.path.split(path)[0], 'scr'), 'scr')
        if os.path.exists(path.rsplit('.', 1)[0] + '.xyz'):
            shutil.copyfile(path.rsplit('.', 1)[0] + '.xyz', subname + '.xyz')

    # Find the hfx resampling values that we're ready to generate
    hfx_values_to_generate = []
    existing_resampled_values = glob.glob(os.path.join(hfx_path, name + '_*'))
    for existing in existing_resampled_values:
        hfx = int(existing.rsplit('_', 1)[1])
        subname = name + '_' + str(hfx)
        outfile_path = os.path.join(existing, subname + '.out')
        if os.path.exists(outfile_path):
            if read_outfile(outfile_path)['finished']:
                hfx_values_to_generate.append(hfx - 5)
                hfx_values_to_generate.append(hfx + 5)

    hfx_values_to_generate = list(set(hfx_values_to_generate))
    hfx_values_to_generate = [i for i in hfx_values_to_generate if i in hfx_values]

    # Now generate the additional hfx resampling values
    jobscripts = []
    for hfx in hfx_values_to_generate:
        subname = name + '_' + str(hfx)
        if os.path.exists(os.path.join(hfx_path, subname)):  # skip over values that we've already done
            continue

        os.mkdir(os.path.join(hfx_path, subname))
        os.chdir(os.path.join(hfx_path, subname))

        higher_hfx = subname.rsplit('_', 1)[0] + '_' + str(int(subname.rsplit('_', 1)[1]) + 5)
        lower_hfx = subname.rsplit('_', 1)[0] + '_' + str(int(subname.rsplit('_', 1)[1]) - 5)
        if os.path.exists(os.path.join(hfx_path, higher_hfx)):
            source_dir = os.path.join(hfx_path, higher_hfx)
        else:
            source_dir = os.path.join(hfx_path, lower_hfx)

        optimxyz = os.path.join(source_dir, 'scr', 'optim.xyz')
        extract_optimized_geo(optimxyz)

        shutil.copy(os.path.join(source_dir, 'scr', 'optimized.xyz'), subname + '.xyz')
        if infile_dict['spinmult'] == 1:
            shutil.copy(os.path.join(source_dir, 'scr', 'c0'), 'c0')
            write_jobscript(subname, custom_line='# -fin c0')
        elif infile_dict['spinmult'] != 1:
            shutil.copyfile(os.path.join(source_dir, 'scr', 'ca0'), os.path.join('ca0'))
            shutil.copyfile(os.path.join(source_dir, 'scr', 'cb0'), os.path.join('cb0'))
            write_jobscript(subname, custom_line=['# -fin ca0\n', '# -fin cb0\n'])

        local_infile_dict = copy.copy(infile_dict)
        local_infile_dict['guess'] = True
        local_infile_dict['hfx'] = hfx / 100.
        local_infile_dict['name'] = subname
        write_input(local_infile_dict)
        jobscripts.append(os.path.join(os.getcwd(), subname + '_jobscript'))

    os.chdir(home)

    return jobscripts
