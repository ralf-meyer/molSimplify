import os
import glob
import copy
import numpy as np
import subprocess
import pandas as pd
import shutil
import time
from molSimplify.job_manager.classes import resub_history,textfile

def not_nohup(path):
    #The nohup.out file gets caught in the find statement
    #use this function so that we only get TeraChem.outs
    endpath = os.path.split(path)[-1]
    if 'nohup.out' in endpath:
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
    for key in dictionary.keys():
        if type(dictionary[key])==list:
            for entry in dictionary[key]:
                if entry in new_dict.keys():
                    raise Exception('Dictionary inversion failed, values do not serve as unique keys')
                new_dict[entry] = key
        else:
            if dictionary[key] in new_dict.keys():
                raise Exception('Dictionary inversion failed, values do not serve as unique keys')
            new_dict[dictionary[key]]=key
    return new_dict

def call_bash(string, error = False,version = 1):
    if version == 1:
        p = subprocess.Popen(string.split(),stdout = subprocess.PIPE,stderr = subprocess.PIPE)
    elif version == 2:
        p = subprocess.Popen(string,stdout = subprocess.PIPE,stderr = subprocess.PIPE, shell = True)
    out,err = p.communicate() 
 
    out = out.split('\n')
    if out[-1] == '':
        out = out[:-1]
    
    if error:
        return out,err
    else:
        return out

def convert_to_absolute_path(path): 
    if path[0] != '/':
        path = os.path.join(os.getcwd(),path)
    
    return path


def list_active_jobs():
    #  @return A list of active jobs for the current user. By job name
    
    job_report = textfile() 
    job_report.lines = call_bash("qstat -r")
    
    name = job_report.wordgrab('jobname:',2)
        
    return name

def check_completeness(directory = 'in place', max_resub = 5):
    ## Takes a directory, returns lists of finished, failed, and in-progress jobs
    outfiles = find('*.out',directory)
    outfiles = filter(not_nohup,outfiles)
    
    results_tmp = map(read_outfile,outfiles)
    results_tmp = zip(outfiles,results_tmp)
    results_dict = dict()
    for outfile,tmp in results_tmp:
        results_dict[outfile]=tmp
        
    try:
        active_jobs = list_active_jobs()
    except:
        active_jobs = [] #Try/accept necesssary to use locally
        print('WARNING, active jobs list not queried!')
    
    
    def check_finished(path,results_dict=results_dict):
        #Return True if the outfile corresponds to a complete job, False otherwise
        if results_dict[path]['finished']:
            return True
        else:
            return False
           
            
    def check_active(path,active_jobs=active_jobs):
        #Given a path, checks if it's in the queue currently:
        name = os.path.split(path)[-1]
        name = name.rsplit('.',1)[0]
        if name in active_jobs:
            return True
        else:
            return False
                
            
    def check_needs_resub(path):
        if os.path.isfile(path.rsplit('.',1)[0]+'.pickle'):
            history = resub_history()
            history.read(path)
            return history.needs_resub
        else:
            return False
            
    def check_waiting(path):
        if os.path.isfile(path.rsplit('.',1)[0]+'.pickle'):
            history = resub_history()
            history.read(path)
            if history.waiting:
                return True
        return False
    
    def grab_waiting(path):
        if os.path.isfile(path.rsplit('.',1)[0]+'.pickle'):
            history = resub_history()
            history.read(path)
            return history.waiting
        raise Exception('Attempting to grab a "waiting" criteria that does not exist')
            
    def check_chronic_failure(path):
        if os.path.isfile(path.rsplit('.',1)[0]+'.pickle'):
            history = resub_history()
            history.read(path)
            if history.resub_number > max_resub-1:
                return True
        else:
            return False
            
    def check_spin_contaminated(path,results_dict=results_dict):
        results = results_dict[path]
        if results['finished']:
            if type(results['s_squared_ideal']) == float:
                if abs(results['s_squared'] - results['s_squared_ideal']) > 1:
                    return True
        return False
        
    def check_thermo_grad_error(path,results_dict=results_dict):
        results = results_dict[path]
        if results['thermo_grad_error']:
            return True
        else:
            return False
        
    active_jobs = filter(check_active,outfiles)
    finished = filter(check_finished,outfiles)
    needs_resub = filter(check_needs_resub,outfiles)
    waiting = filter(check_waiting,outfiles)
    spin_contaminated = filter(check_spin_contaminated,outfiles)
    thermo_grad_errors = filter(check_thermo_grad_error,outfiles)
    chronic_errors = filter(check_chronic_failure,outfiles)
    errors = list(set(outfiles) - set(active_jobs) - set(finished))
    
    #Sort out conflicts in order of reverse priority
    #A job only gets labelled as finished if it's in no other category
    #A job always gets labelled as active if it fits that criteria, even if it's in every other category too
    finished = list(set(finished)- set(spin_contaminated) - set(needs_resub) - set(errors) - set(thermo_grad_errors) - set(waiting) - set(chronic_errors) - set(active_jobs))
    spin_contaminated = list(set(spin_contaminated) - set(needs_resub) - set(errors) - set(thermo_grad_errors) - set(waiting) - set(chronic_errors) - set(active_jobs))
    needs_resub = list(set(needs_resub) - set(errors) - set(thermo_grad_errors) - set(waiting) - set(chronic_errors) - set(active_jobs))
    errors = list(set(errors) - set(thermo_grad_errors) - set(waiting) - set(chronic_errors) - set(active_jobs))
    thermo_grad_errors = list(set(thermo_grad_errors) - set(waiting) - set(chronic_errors) - set(active_jobs))
    waiting = list(set(waiting) - set(chronic_errors) - set(active_jobs))
    chronic_errors = list(set(chronic_errors) - set(active_jobs))
    active_jobs = list(set(active_jobs))
    
    results = {'Finished':finished,'Active':active_jobs,'Error':errors,'Resub':needs_resub,
            'Spin_contaminated':spin_contaminated, 'Chronic_error':chronic_errors, 
            'Thermo_grad_error':thermo_grad_errors, 'Waiting':waiting}
    
    inverted_results = invert_dictionary(results)
    waiting = [{i:grab_waiting(i)} for i in waiting]
    results['Waiting'] = waiting
    
    return results,inverted_results
    

def find(key,directory = 'in place'):
    ## Looks for all files with a matching key in their name within directory
    #  @return A list of paths
    if directory == 'in place':
        directory = os.getcwd()
    
    bash = 'find '+directory+' -name '+key

    return call_bash(bash)
    
def qsub(jobscript_list):
    ## Takes a list of paths to jobscripts (like output by find) and submits them with qsub
    
    #Handles cases where a single jobscript is submitted instead of a list
    if type(jobscript_list) != list:
        jobscript_list = [jobscript_list]
    
    for i in jobscript_list:
        home = os.getcwd()
        os.chdir(os.path.split(i)[0])
        jobscript = os.path.split(i)[1]
        stdout = call_bash('qsub '+jobscript)
        print stdout
        os.chdir(home)
        time.sleep(1)
        
    return None

def read_outfile(outfile_path):
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
                print 'Warning! Nohup file caught in outfile processing'
                print outfile_path
                counter = 0
            else:
                raise ValueError('.out file type not recognized!')
    output_type = ['TeraChem','ORCA'][counter]
    if output_type == 'ORCA':
        raise Exception(outfile_path)
    
    name = None
    finished = None
    charge = None
    finalenergy = None
    min_energy = None
    s_squared = None
    s_squared_ideal = None
    scf_error = False
    time = None
    thermo_grad_error = False

    if output_type == 'TeraChem':
        
        name,charge = output.wordgrab(['Startfile','charge:'],[4,2],first_line=True)
        name = name.rsplit('.',1)[0]
        finalenergy,s_squared,s_squared_ideal,time,thermo_grad_error = output.wordgrab(['FINAL','S-SQUARED:','S-SQUARED:','processing',
                                                                                        'Maximum component of gradient is too large'],[2,2,4,3,0],last_line=True)
        if thermo_grad_error:
            thermo_grad_error = True
        else:
            thermo_grad_error = False
        if s_squared_ideal:
            s_squared_ideal = float(s_squared_ideal.strip(')'))
            
        min_energy = output.wordgrab('FINAL',2,min_value = True)[0]
        
        is_finished = output.wordgrab(['finished:'],'whole_line',last_line=True)[0]
        if is_finished:
            if is_finished[0] == 'Job' and is_finished[1] == 'finished:':
                finished = is_finished[2:]
        
        is_scf_error = output.wordgrab('DISS','whole line')
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
    return return_dict

def read_infile(outfile_path):
    root = outfile_path.rsplit('.',1)[0]
    inp = textfile(root+'.in')
    charge,spinmult,solvent,run_type,levelshifta,levelshiftb,method,hfx,basis = inp.wordgrab(['charge ', 'spinmult ','pcm ',
                                                                                        'run ','levelshiftvala ',
                                                                                        'levelshiftvalb ','method ',
                                                                                        'HFX ', 'basis '],[1,1,0,1,1,1,1,1,1],last_line=True)
    charge,spinmult = int(charge),int(spinmult)
    if solvent:
        solvent = True
    else:
        solvent = False
    if method[1] == 'u':
        method = method[1:]
        
    convergence_thresholds = inp.wordgrab(['min_converge_gmax ','min_converge_grms ','min_converge_dmax ','min_converge_drms ','min_converge_e ','convthre '],
                                          [1]*6,last_line=True)
    return charge,spinmult,solvent,run_type,levelshifta,levelshiftb,method,hfx,basis,convergence_thresholds

def read_charges(PATH):
    #Takes the path to either the outfile or the charge_mull.xls and returns the charges
    PATH = convert_to_absolute_path(PATH)
    if len(PATH.rsplit('.',1)) > 1:
        if PATH.rsplit('.',1)[1] == 'out':
            PATH = os.path.join(os.path.split(PATH)[0],'scr','charge_mull.xls')
    try:
        charge_mull = textfile(PATH)
        split_lines = [i.split() for i in charge_mull.lines]
        charges = [i[1] + ' '+ i[2] for i in split_lines]
        return charges
    except:
        return []

def read_mullpop(PATH):
    #Takes the path to either the outfile or the mullpop and returns the mullikan populations
    PATH = convert_to_absolute_path(PATH)
    if len(PATH.rsplit('.',1)) > 1:
        if PATH.rsplit('.',1)[1] == 'out':
            PATH = os.path.join(os.path.split(PATH)[0],'scr','mullpop')
    try:
        mullpop = textfile(PATH)
        split_lines = [i.split() for i in mullpop.lines]
        if len(split_lines[2]) == 6:
            pops = [i[1] + ' ' + i[5] for i in split_lines[1:-2]]
        else:
            pops = [i[1] + ' ' + i[5] + ' ' + i[9] for i in split_lines[2:-2]]
            
        return pops
    except:
        return []
    
def create_summary(directory='in place'):
    #Returns a pandas dataframe which summarizes all outfiles in the directory, defaults to cwd
    outfiles = find('*.out',directory)
    outfiles = filter(not_nohup,outfiles)
    results = map(read_outfile,outfiles)
    summary = pd.DataFrame(results)
    
    return summary


def prep_vertical_ip(path, solvent = False):
    #Given a path to the outfile of a finished run, this preps the files for a corresponding vertical IP run
    #Returns a list of the PATH(s) to the jobscript(s) to start the vertical IP calculations(s)
    home = os.getcwd()
    path = convert_to_absolute_path(path)
    
    results = read_outfile(path)
    if not results['finished']:
        raise Exception('This calculation does not appear to be complete! Aborting...')
        
    charge,spin,solvent,run_type,levelshifta,levelshiftb,method,hfx,basis,convergence_thresholds = read_infile(path)
    
    if spin == 1:
        new_spin = [2]
    else:
        new_spin = [spin-1,spin+1]
    
    base = os.path.split(path)[0]
    
    optimxyz = os.path.join(base,'scr','optim.xyz')
    extract_optimized_geo(optimxyz)
    
    jobscripts = []
    for calc in new_spin:
        name = results['name']+'_vertIP_'+str(calc)
        PATH = os.path.join(base,name)
        if os.path.isdir(PATH):
            jobscripts.append('File for vert EA spin '+str(calc)+'already exists')
        else:
            os.mkdir(PATH)
            os.chdir(PATH)
            
            shutil.copyfile(os.path.join(base,'scr','optimized.xyz'),os.path.join(PATH,name+'.xyz'))
            
            write_input(name,charge+1,calc,solvent = solvent, guess = False, 
                run_type = 'energy', method = method, levela = levelshifta, 
                levelb = levelshiftb, thresholds = convergence_thresholds, hfx = hfx, basis = basis)
            write_jobscript(name)
            
            jobscripts.append(os.path.join(PATH,name+'_jobscript'))
    
    os.chdir(home)
    
    return jobscripts


def prep_vertical_ea(path, solvent = False):
    #Given a path to the outfile of a finished run, this preps the files for a corresponding vertical IP run
    #Returns a list of the PATH(s) to the jobscript(s) to start the vertical IP calculations(s)
    home = os.getcwd()
    path = convert_to_absolute_path(path)
    
    results = read_outfile(path)
    if not results['finished']:
        raise Exception('This calculation does not appear to be complete! Aborting...')
    
    charge,spin,solvent,run_type,levelshifta,levelshiftb,method,hfx,basis,convergence_thresholds = read_infile(path)
    
    if spin == 1:
        new_spin = [2]
    else:
        new_spin = [spin-1,spin+1]
    
    base = os.path.split(path)[0]
    
    optimxyz = os.path.join(base,'scr','optim.xyz')
    extract_optimized_geo(optimxyz)
    
    jobscripts = []
    for calc in new_spin:
        name = results['name']+'_vertEA_'+str(calc)
        PATH = os.path.join(base,name)
        if os.path.isdir(PATH):
            jobscripts.append('File for vert EA spin '+str(calc)+'already exists')
        else:
            os.mkdir(PATH)
            os.chdir(PATH)
            
            shutil.copyfile(os.path.join(base,'scr','optimized.xyz'),os.path.join(PATH,name+'.xyz'))
            
            write_input(name,charge-1,calc,solvent = solvent, guess = False, 
                run_type = 'energy', method = method, levela = levelshifta, 
                levelb = levelshiftb, thresholds = convergence_thresholds, hfx = hfx, basis = basis)
            write_jobscript(name)
            
            jobscripts.append(os.path.join(PATH,name+'_jobscript'))
    
    os.chdir(home)
    
    return jobscripts

def prep_solvent_sp(path):
    #Given a path to the outfile of a finished run, this preps the files for a single point solvent run
    #Uses the wavefunction from the gas phase calculation as an initial guess
    #Returns a list of the PATH(s) to the jobscript(s) to start the solvent sp calculations(s)
    home = os.getcwd()
    path = convert_to_absolute_path(path)
    
    results = read_outfile(path)
    if not results['finished']:
        raise Exception('This calculation does not appear to be complete! Aborting...')
    
    charge,spin,solvent,run_type,levelshifta,levelshiftb,method,hfx,basis,convergence_thresholds = read_infile(path)
    
    base = os.path.split(path)[0]
    
    optimxyz = os.path.join(base,'scr','optim.xyz')
    extract_optimized_geo(optimxyz)
    
    #Now, start generating the new directory
    name = results['name']+'_solventSP'
    PATH = os.path.join(base,name)
    if os.path.isdir(PATH):
        return ['Directory for solvent single point already exists']
        
    os.mkdir(PATH)
    os.chdir(PATH)    
    
    shutil.copyfile(os.path.join(base,'scr','optimized.xyz'),os.path.join(PATH,name+'.xyz'))
    if spin == 1:
        shutil.copyfile(os.path.join(base,'scr','c0'),os.path.join(PATH,'c0'))
        write_jobscript(name,custom_line = '# -fin c0')
    if spin != 1:
        shutil.copyfile(os.path.join(base,'scr','ca0'),os.path.join(PATH,'ca0'))
        shutil.copyfile(os.path.join(base,'scr','cb0'),os.path.join(PATH,'cb0'))
        write_jobscript(name,custom_line = ['# -fin ca0\n','# -fin cb0\n'])
    
    write_input(name,charge,spin,solvent = True, guess = True, 
                run_type = 'energy', method = method, levela = levelshifta, 
                levelb = levelshiftb, hfx = hfx,thresholds = convergence_thresholds, basis = basis)
    
    os.chdir(home)
    
    return [os.path.join(PATH,name+'_jobscript')]
    
def prep_thermo(path):
    #Given a path to the outfile of a finished run, this preps the files for a thermo calculation
    #Uses the wavefunction from the previous calculation as an initial guess
    #Returns a list of the PATH(s) to the jobscript(s) to start the solvent sp calculations(s)
    home = os.getcwd()
    path = convert_to_absolute_path(path)
    
    results = read_outfile(path)
    if not results['finished']:
        raise Exception('This calculation does not appear to be complete! Aborting...')
    
    charge,spin,solvent,run_type,levelshifta,levelshiftb,method,hfx,basis,convergence_thresholds = read_infile(path)
    
    base = os.path.split(path)[0]
    
    optimxyz = os.path.join(base,'scr','optim.xyz')
    extract_optimized_geo(optimxyz)
    
    #Now, start generating the new directory
    name = results['name']+'_thermo'
    PATH = os.path.join(base,name)
    if os.path.isdir(PATH):
        return ['Thermo Calculation Directory already exists']
        
    os.mkdir(PATH)
    os.chdir(PATH)    
    
    shutil.copyfile(os.path.join(base,'scr','optimized.xyz'),os.path.join(PATH,name+'.xyz'))
    if spin == 1:
        shutil.copyfile(os.path.join(base,'scr','c0'),os.path.join(PATH,'c0'))
        write_jobscript(name,custom_line = '# -fin c0')
    if spin != 1:
        shutil.copyfile(os.path.join(base,'scr','ca0'),os.path.join(PATH,'ca0'))
        shutil.copyfile(os.path.join(base,'scr','cb0'),os.path.join(PATH,'cb0'))
        write_jobscript(name,custom_line = ['# -fin ca0\n','# -fin cb0\n'])
    
    write_input(name,charge,spin,solvent = solvent, guess = True, 
                run_type = 'frequencies', method = method, levela = levelshifta, 
                levelb = levelshiftb, hfx = hfx, basis = basis)
    
    os.chdir(home)
    
    return [os.path.join(PATH,name+'_jobscript')]
    
def prep_ultratight(path):
    #Given a path to the outfile of a finished run, this preps a run with tighter convergence criteria
    #Uses the wavefunction and geometry from the previous calculation as an initial guess
    #Returns a list of the PATH(s) to the jobscript(s) to start the solvent sp calculations(s)
    home = os.getcwd()
    path = convert_to_absolute_path(path)
    
    results = read_outfile(path)
    if not results['finished']:
        raise Exception('This calculation does not appear to be complete! Aborting...')
    
    charge,spin,solvent,run_type,levelshifta,levelshiftb,method,hfx,basis,convergence_thresholds = read_infile(path)
    
    base = os.path.split(path)[0]
    
    optimxyz = os.path.join(base,'scr','optim.xyz')
    extract_optimized_geo(optimxyz)
    
    #Now, start generating the new directory
    name = results['name']+'_ultratight'
    PATH = os.path.join(base,name)
    
    if not os.path.isdir(PATH): #First time that ultratight has been run, create necessary files
        os.mkdir(PATH)
        os.chdir(PATH)    
        
        if os.path.exists(name+'.in') or os.path.exists(name+'.out') or os.path.exists(name+'_jobscript'):
            raise Exception('This tightened convergence run appears to already exist. Aborting...')
            
        shutil.copyfile(os.path.join(base,'scr','optimized.xyz'),os.path.join(PATH,name+'.xyz'))
        if spin == 1:
            shutil.copyfile(os.path.join(base,'scr','c0'),os.path.join(PATH,'c0'))
            write_jobscript(name,custom_line = '# -fin c0')
        if spin != 1:
            shutil.copyfile(os.path.join(base,'scr','ca0'),os.path.join(PATH,'ca0'))
            shutil.copyfile(os.path.join(base,'scr','cb0'),os.path.join(PATH,'cb0'))
            write_jobscript(name,custom_line = ['# -fin ca0\n','# -fin cb0\n'])
        
        criteria = ['2.25e-04','1.5e-04','0.9e-03','0.6e-03','0.5e-06','1.5e-05']
        
        write_input(name,charge,spin,solvent = solvent, guess = True, 
                run_type = run_type, method = method, levela = levelshifta, 
                levelb = levelshiftb, thresholds = criteria, hfx = hfx, basis = basis)
        
        #Make an empty .out file to prevent the resubmission module from mistakenly submitting this job twice
        f = open(name+'.out','w')
        f.close()

        os.chdir(home)
        
        return [os.path.join(PATH,name+'_jobscript')]
    
    else: #This has been run before, further tighten the convergence criteria
        os.chdir(PATH)
        charge,spin,solvent,run_type,levelshifta,levelshiftb,method,hfx,basis,criteria = read_infile(os.path.join(PATH,name+'.out'))
        criteria = [str(i/2.) for i in criteria]
        
        write_input(name,charge,spin,solvent = solvent, guess = True, 
                run_type = run_type, method = method, levela = levelshifta, 
                levelb = levelshiftb, thresholds = criteria, hfx = hfx, basis = basis)
        
        
        os.chdir(home)
        
        return [os.path.join(PATH,name+'_jobscript')]
        
def prep_hfx_resample(path,hfx_values = [0,0.05,0.10,0.15,0.20,0.25,0.30]):
    #Given a path to the outfile of a finished run, this preps the files for hfx resampling
    #Uses the wavefunction from the gas phase calculation as an initial guess
    #Returns a list of the PATH(s) to the jobscript(s) to start the resampling calculations(s)
    home = os.getcwd()
    path = convert_to_absolute_path(path)
    base = os.path.split(path)[0]
    
    results = read_outfile(path)
    if not results['finished']:
        raise Exception('This calculation does not appear to be complete! Aborting...')
    
    #Check the state of the calculation and ensure than hfx resampling is valid
    charge,spin,solvent,run_type,levelshifta,levelshiftb,method,hfx,basis,convergence_thresholds = read_infile(path)
    if method != 'b3lyp':
        raise Exception('HFX resampling may not behave well for methods other than b3lyp!')
    if not hfx:
        hfx = 0.20
    if hfx not in hfx_values:
        raise Exception('HFX resampling list does not contain the original hfx value!')
    
    #Now, start generating the base directory to hold all the hfx resampling values
    name = results['name']+'_HFXresampling'
    hfx_path = os.path.join(base,name)
    if not os.path.isdir(hfx_path):
        os.mkdir(hfx_path)
    os.chdir(hfx_path)    
    
    #Make the directory for the original calculation
    subname = name+'_'+str(int(100*hfx))
    PATH = os.path.join(hfx_path,subname)
    if not os.path.isdir(PATH):
        os.mkdir(PATH)
    os.chdir(PATH)
    
    if not os.path.exists(os.path.join(PATH,subname+'.out')):
        shutil.copyfile(path,subname+'.out')
        shutil.copyfile(path.rsplit('.',1)[0]+'_jobscript',subname+'_jobscript')
        shutil.copyfile(path.rsplit('.',1)[0]+'.in',subname+'.in')
        shutil.copyfile(path.rsplit('.',1)[0]+'.xyz',subname+'.xyz')
        shutil.copytree(os.path.join(os.path.split(path)[0],'scr'),'scr')
    
    #Find the hfx resampling values that we're ready to generate
    hfx_values_to_generate = []
    existing_resampled_values = glob.glob(os.path.join(hfx_path,name+'_*'))
    for existing in existing_resampled_values:
        hfx = float(existing.rsplit('_',1)[1])/100.
        subname = name+'_'+str(int(100*hfx))
        outfile_path = os.path.join(existing,subname+'.out')
        if os.path.exists(outfile_path):
            if read_outfile(outfile_path)['finished']:
                hfx_values_to_generate.append(hfx-0.05)
                hfx_values_to_generate.append(hfx+0.05)
    
    hfx_values_to_generate = [int(100.*i)/100. for i in hfx_values_to_generate] #We're about to check equality to a float, this makes sure that it behaves well
    hfx_values_to_generate = list(set(hfx_values_to_generate))
    hfx_values_to_generate = [i for i in hfx_values_to_generate if i in hfx_values]
    
    #Now generate the additional hfx resampling values
    jobscripts = []
    for hfx in hfx_values_to_generate:
        subname = name+'_'+str(int(100*hfx))
        if os.path.exists(os.path.join(hfx_path,subname)): #skip over values that we've already done
            continue
            
        os.mkdir(os.path.join(hfx_path,subname))
        os.chdir(os.path.join(hfx_path,subname))
        
        higher_hfx = subname.rsplit('_',1)[0]+'_'+str(int(subname.rsplit('_',1)[1])+5)
        lower_hfx = subname.rsplit('_',1)[0]+'_'+str(int(subname.rsplit('_',1)[1])-5)
        if os.path.exists(os.path.join(hfx_path,higher_hfx)):
            source_dir = os.path.join(hfx_path,higher_hfx)
        else:
            source_dir = os.path.join(hfx_path,lower_hfx)
            
        optimxyz = os.path.join(source_dir,'scr','optim.xyz')
        extract_optimized_geo(optimxyz)
        
        shutil.copy(os.path.join(source_dir,'scr','optimized.xyz'),subname+'.xyz')
        if spin == 1:
            shutil.copy(os.path.join(source_dir,'scr','c0'),'c0')
            write_jobscript(name,custom_line = '# -fin c0')
        if spin != 1:
            shutil.copyfile(os.path.join(source_dir,'scr','ca0'),os.path.join('ca0'))
            shutil.copyfile(os.path.join(source_dir,'scr','cb0'),os.path.join('cb0'))
            write_jobscript(name,custom_line = ['# -fin ca0\n','# -fin cb0\n'])
    
        write_input(name,charge,spin,solvent = solvent, guess = True, 
            run_type = run_type, method = method, levela = levelshifta, 
            levelb = levelshiftb, hfx = hfx,thresholds = convergence_thresholds, basis = basis)
        jobscripts.append(os.path.join(os.getcwd(),subname+'_jobscript'))
    
    os.chdir(home)
    
    return jobscripts

def extract_optimized_geo(PATH, custom_name = False):
    #Given the path to an optim.xyz file, this will extract optimized.xyz, which contains only the last frame
    #The file is written to the same directory as contained optim.xyz
    
    optim = open(PATH,'r')
    lines = optim.readlines()
    optim.close()
    lines.reverse()
    if len(lines) == 0:
        lines = []
        print('optim.xyz is emptry for: ' +PATH)
    else:
        for i in range(len(lines)):
            if 'frame' in lines[i].split():
                break
        lines = lines[:i+2]
        lines.reverse()
    homedir = os.path.split(PATH)[0]

    if custom_name:
        basedir = os.path.split(homedir)[0]
        xyz = glob.glob(os.path.join(basedir,'*.xyz'))[0]
        xyz = os.path.split(xyz)[-1]
        name = xyz[:-4]+'_optimized.xyz'
    else:
        name = 'optimized.xyz'
    
    optimized = open(os.path.join(homedir,name),'w')
    for i in lines:
        optimized.write(i)
    optimized.close()
    
    return lines
    
def pull_optimized_geos(PATHs = []):
    #Given a PATH or list of PATHs to optim.xyz files, these will grab the optimized geometries and move them into a local folder called optimized_geos
    if type(PATHs) != list:
        PATHs = [PATHs]
    if len(PATHs) == 0:
        PATHs = find('optim.xyz') #Default behaviour will pull all optims in the current directory
        
    if os.path.isdir('opimized_geos'):
        raise Exception('optimized_geos already exists!')
    
    os.mkdir('optimized_geos')
    home = os.getcwd()
    for Path in PATHs:
        extract_optimized_geo(Path)
        scr_path = os.path.split(Path)[0]
        homedir = os.path.split(scr_path)[0]
        initial_xyz = glob.glob(os.path.join(homedir,'*.xyz'))
        if len(initial_xyz) != 1:
            print 'Name could not be identified for: '+Path
            print 'Naming with a random 6 digit number'
            name = str(np.random.randint(999999))+'_optimized.xyz'
        else:
            initial_xyz = initial_xyz[0]
            name = os.path.split(initial_xyz)[-1]
            name = name.rsplit('.',1)[0]
            name += '_optimized.xyz'
        
        shutil.move(os.path.join(os.path.split(Path)[0],'optimized.xyz'),os.path.join(home,'optimized_geos',name))
    
def write_input(name,charge,spinmult,run_type = 'energy', method = 'b3lyp', solvent = False, 
                guess = False, custom_line = None, levela = 0.25, levelb = 0.25,
                thresholds = None, basis = 'lacvps_ecp', hfx = None,alternate_coordinates = False):
    #Writes a generic terachem input file
    #solvent indicates whether to set solvent calculations True or False
    
    if spinmult != 1:
        method = 'u'+method
        
    if alternate_coordinates:
        coordinates = alternate_coordinates
    else:
        coordinates = name+'.xyz'
        
    input_file = open(name+'.in','w')
    text = ['levelshiftvalb '+str(levelb)+'\n',
            'levelshiftvala '+str(levela)+'\n',
            'nbo yes\n',
            'run '+run_type+'\n',
            'scf diis+a\n',
            'coordinates '+coordinates+'\n',
            'levelshift yes\n',
            'gpus 1\n',
            'spinmult '+str(spinmult) +'\n',
            'scrdir ./scr\n',
            'basis '+basis+'\n',
            'timings yes\n',
            'charge '+str(charge)+'\n',
            'method '+method+'\n',
            'new_minimizer yes\n',
            'end']
    
    if custom_line:
        text = text[:15]+[custom_line+'\n']+text[15:]
    if guess:
        if type(guess) == bool:
            if spinmult == 1:
                guess = 'c0'
            else:
                guess = 'ca0 cb0'
        text = text[:-1] + ['guess ' + guess + '\n',
                            'end']
    if run_type != 'ts':
        text = text[:-1] + ['maxit 500\n',
                            'end']
                            
    if type(thresholds) == list:
        if thresholds[0]:
            tight_thresholds ="min_converge_gmax "+thresholds[0]+"\nmin_converge_grms "+thresholds[1]+"\nmin_converge_dmax "+thresholds[2]+"\nmin_converge_drms "+thresholds[3]+"\nmin_converge_e "+thresholds[4]+"\nconvthre "+thresholds[5]
            text = text[:-1]+['\n',tight_thresholds+'\n','end']
    
    if hfx:
        text = text[:-1] + ['\n',
                            'HFX '+str(hfx)+'\n',
                            '\n']
                            
    if solvent: #Adds solvent correction, if requested
        text = text[:-1] + ['\n',
                            'pcm cosmo\n',
                            'epsilon 80\n',
                            'pcm_radii read\n',
                            'pcm_radii_file /home/harperd/pcm_radii\n',
                            'end']
    for lines in text:
        input_file.write(lines)
    input_file.close()
    
def write_jobscript(name,custom_line = None,alternate_infile = False,alternate_coordinates = False,time_limit='96:00:00',sleep= False):
    #Writes a generic terachem jobscript
    #custom line allows the addition of extra lines just before the export statement
    if not alternate_infile:
        alternate_infile = name
        
    if alternate_coordinates:
        coordinates = alternate_coordinates
    else:
        coordinates = name+'.xyz'
        
    jobscript = open(name+'_jobscript','w')
    text = ['#$ -S /bin/bash\n',
            '#$ -N '+name+'\n',
            '#$ -cwd\n',
            '#$ -R y\n',
            '#$ -l h_rt='+time_limit+'\n',
            '#$ -l h_rss=8G\n',
            '#$ -q gpus|gpusnew|gpusnewer\n',
            '#$ -l gpus=1\n',
            '#$ -pe smp 1\n',
            '# -fin '+alternate_infile+'.in\n',
            '# -fin ' + coordinates + '\n',
            '# -fout scr/\n',
            'export OMP_NUM_THREADS=1\n',
            'terachem '+alternate_infile+'.in '+'> $SGE_O_WORKDIR/' + name + '.out\n']
    if custom_line:
        if type(custom_line) == list:
            text = text[:12]+custom_line+text[12:]
        else:
            text = text[:12]+[custom_line+'\n']+text[12:]
    if sleep:
        text = text+['sleep 300\n'] #For very short jobs, sleep for 5 minutes after completion to keep the queue happy
    for i in text:
        jobscript.write(i)
    jobscript.close()
