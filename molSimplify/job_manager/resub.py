#!/usr/bin/env python

import os
import glob
import numpy as np
import subprocess
import shutil
import time
import molSimplify.job_manager.tools as tools
import molSimplify.job_manager.moltools as moltools
import time
import sys
from molSimplify.job_manager.classes import resub_history

## Save the outfile within the resub_history pickel object
def save_run(outfile_path):
    history = resub_history()
    history.read(outfile_path)
    
    f = open(outfile_path,'r')
    out_lines = f.readlines()
    f.close()
    history.outfiles.append(out_lines)
    
    infile_path = outfile_path.rsplit('.',1)[0]+'.in'
    f = open(infile_path,'r')
    in_lines = f.readlines()
    f.close()
    history.infiles.append(in_lines)
    
    jobscript_path = outfile_path.rsplit('.',1)[0]+'_jobscript'
    f = open(jobscript_path,'r')
    job_lines = f.readlines()
    f.close()
    history.jobscripts.append(job_lines)
    
    xyz_path = outfile_path.rsplit('.',1)[0]+'.xyz'
    f = open(xyz_path,'r')
    xyz_lines = f.readlines()
    f.close()
    history.xyzs.append(xyz_lines)
    
    history.save()

## Archive the scr file so it isn't overwritten in future resubs
#  @param rewrite_inscr Determines whether to copy this runs wfn and optimized geometry to the inscr directory
def save_scr(outfile_path, rewrite_inscr = True):
    root = os.path.split(outfile_path)[0]
    scr_path = os.path.join(root,'scr')
    
    if os.path.isdir(scr_path):
        #extract the optimized geometry, if it exists
        optim = glob.glob(os.path.join(scr_path,'optim.xyz'))
        if len(optim) > 0:
            tools.extract_optimized_geo(optim[0])
        
        if rewrite_inscr:
            #save the files necessary to resub the job to a folder called inscr
            save_paths = []
            save_paths.extend(glob.glob(os.path.join(scr_path,'c0')))
            save_paths.extend(glob.glob(os.path.join(scr_path,'ca0')))
            save_paths.extend(glob.glob(os.path.join(scr_path,'cb0')))
            save_paths.extend(glob.glob(os.path.join(scr_path,'optimized.xyz')))
            if os.path.isdir(os.path.join(root,'inscr')):
                shutil.rmtree(os.path.join(root,'inscr'))
            os.mkdir(os.path.join(root,'inscr'))
            for path in save_paths:
                shutil.copy(path,os.path.join(root,'inscr',os.path.split(path)[-1]))
        
        #archive the scr under a new name so that we can write a new one
        old_scrs = glob.glob(scr_path+'_*')
        old_scrs = [int(i[-1]) for i in old_scrs]
        if len(old_scrs) > 0:
            new_scr = str(max(old_scrs)+1)
        else:
            new_scr = '0'
        shutil.move(scr_path,scr_path+'_'+new_scr)
        
        return scr_path+'_'+new_scr

#Read the global and local configure files to determine the derivative jobs requested and the settings for job recovery
#The global configure file should be in the same directory where resub() is called
#The local configure file should be in the same directory as the .out file
def read_configure(home_directory,outfile_path):
    
    def load_configure_file(directory):
        def strip_new_line(string):
            if string[-1] == '\n':
                return string[:-1]
            else:
                return string
                
        if directory == 'in place':
            directory = os.getcwd()
            
        configure = os.path.join(directory,'configure')
        if os.path.isfile(configure):
            f = open(configure,'r')
            configure = f.readlines()
            f.close()
            configure = map(strip_new_line,configure)
            return configure
        else:
            return []
    
    home_configure = load_configure_file(home_directory)
    if outfile_path:
        local_configure = load_configure_file(os.path.split(outfile_path)[0])
    else:
        local_configure = []
    
    #Determine which derivative jobs are requested
    solvent,vertEA,vertIP,thermo,dissociation,hfx_resample = False,False,False,False,False,False
    if 'solvent' in home_configure or 'Solvent' in home_configure or 'solvent' in local_configure or 'Solvent' in local_configure:
        solvent = True
    if 'vertEA' in home_configure or 'VertEA' in home_configure or 'vertEA' in local_configure or 'VertEA' in local_configure:
        vertEA = True
    if 'vertIP' in home_configure or 'VertIP' in home_configure or 'vertIP' in local_configure or 'VertIP' in local_configure:
        vertIP = True
    if 'thermo' in home_configure or 'Thermo' in home_configure or 'thermo' in local_configure or 'Thermo' in local_configure:
        thermo = True
    if 'dissociation' in home_configure or 'Dissociation' in home_configure or 'dissociation' in local_configure or 'Dissociation' in local_configure:
        dissociation = True
    if 'hfx_resample' in home_configure or 'HFX_resample' in home_configure or 'hfx_resample' in local_configure or 'HFX_resample' in local_configure:
        hfx_resample = True
    
    #Determine global settings for this run
    max_jobs,max_resub,levela,levelb,method,hfx,octahedral,sleep = False,False,False,False,False,False,True,False
    for configure in [home_configure,local_configure]:
        for line in home_configure:
            if 'max_jobs' in line.split(':'):
                max_jobs = int(line.split(':')[-1]) - 1
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
            if 'octahedral' in line.split(':'):
                octahedral = line.split(':')[-1]
            if 'sleep' in line.split(':'):
                sleep = int(line.split(':')[-1])
    #If global settings not specified, choose defaults:
        if not max_jobs:
            max_jobs = 50 - 1 
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
        #Octahedral defaults to True in original variable initiation
                
    return {'solvent':solvent,'vertEA':vertEA,'vertIP':vertIP,'thermo':thermo,'dissociation':dissociation,
            'hfx_resample':hfx_resample,'max_jobs':max_jobs,'max_resub':max_resub,'levela':levela,
            'levelb':levelb,'method':method,'hfx':hfx,'octahedral':octahedral,'sleep':sleep}

def resub(directory = 'in place',max_jobs = 50,max_resub = 5):
    
    #Takes a directory, resubmits errors, scf failures, and spin contaminated cases
    completeness,inverted_completeness = tools.check_completeness(directory,max_resub)
    errors = completeness['Error'] #These are calculations which failed to complete
    need_resub = completeness['Resub'] #These are calculations with level shifts changed or hfx exchange changed
    spin_contaminated = completeness['Spin_contaminated'] #These are finished jobs with spin contaminated solutions
    active = completeness['Active'] #These are jobs which are currently running
    thermo_grad_error = completeness['Thermo_grad_error'] #These are thermo jobs encountering the thermo grad error
    waiting = completeness['Waiting'] #These are jobs which are or were waiting for another job to finish before continuing.
    finished = completeness['Finished']
    
    #Prep derivative jobs such as thermo single points, vertical IP, and ligand dissociation energies
    prep_derivative_jobs(directory,finished)
    
    #Resub scf convergence errors and unidentified errors
    resubmitted = [] #Resubmitted list gets True if the job is submitted or False if not. Contains booleans, not job identifiers.
    for error in errors:
        if len(active)+np.sum(resubmitted) > max_jobs:
            continue
        results = tools.read_outfile(error)
        if results['scf_error']:
            resub_tmp = resub_scf(error)
            if resub_tmp:
                print('SCF error identified in job: '+os.path.split(error)[-1]+' -Resubmitting with adjusted levelshifts')
                print('')
            resubmitted.append(resub_tmp)
        else:
            resub_tmp = simple_resub(error)
            if resub_tmp:
                print('Unidentified error in job: '+os.path.split(error)[-1]+' -Resubmitting')
                print('')
            resubmitted.append(resub_tmp)
            
    #Resub spin contaminated cases
    for error in spin_contaminated:
        if len(active)+np.sum(resubmitted) > max_jobs:
            continue
        resub_tmp = resub_spin(error)
        if resub_tmp:
            print('Spin contamination identified in job: '+os.path.split(error)[-1]+' -Resubmitting with adjusted HFX')
            print('')
        resubmitted.append(resub_tmp)
        
    #Resub jobs with atypical parameters used to aid convergence
    for error in need_resub:
        if len(active)+np.sum(resubmitted) > max_jobs:
            continue
        resub_tmp = clean_resub(error)
        if resub_tmp:
            print('Job '+os.path.split(error)[-1]+' needs to be rerun with typical paramters. -Resubmitting')
            print('')
        resubmitted.append(resub_tmp)
        
    #Create a job with a tighter convergence threshold for failed thermo jobs
    for error in thermo_grad_error:
        if len(active)+np.sum(resubmitted) > max_jobs:
            continue
        resub_tmp = resub_tighter(error)
        if resub_tmp:
            print('Job '+os.path.split(error)[-1]+' needs a better initial geo. Creating a geometry run with tighter convergence criteria')
            print('')
        resubmitted.append(resub_tmp)
        
    #Look at jobs in "waiting," resume them if the job they were waiting for is finished
    #Currently, this should only ever be thermo jobs waiting for an ultratight job
    for waiting_dict in waiting:
        if len(waiting_dict.keys()) > 1:
            raise Exception('Waiting job list improperly constructed')
        job = waiting_dict.keys()[0]
        waiting_for = waiting_dict[job]
        if waiting_for in finished:
	    history = load_history(job)
            history.waiting = None
            history.save()
            results_for_this_job = tools.read_outfile(job)
            if results_for_this_job['thermo_grad_error']:
                resubmitted.append(resub_thermo(job))
            else:
                raise Exception('A method for resuming job: '+job+' is not defined')
        else:
            resubmitted.append(False)
        
    #Submit jobs which haven't yet been submitted
    to_submit = []
    jobscripts = tools.find('*_jobscript')
    for job in jobscripts:
        if not os.path.isfile(job.rsplit('_',1)[0]+'.out'):
            to_submit.append(job)
    submitted = []
    for job in to_submit:
        if len(submitted)+len(active)+np.sum(resubmitted) > max_jobs:
            continue
        print('Initial sumbission for job: '+os.path.split(job)[-1])
        tools.qsub(job)
        submitted.append(True)
    
    number_resubmitted = np.sum(np.array(resubmitted+submitted))
    # ~ print str(number_resubmitted)+' Jobs submitted'
    return int(number_resubmitted)
    
def simple_resub(outfile_path):
    #Resubmits a job without changing parameters. Particularly useful for CUDA errors.
    save_scr(outfile_path, rewrite_inscr = False)
    save_run(outfile_path)
    history = resub_history()
    history.read(outfile_path)
    history.resub_number += 1
    history.notes.append('Resubbed for unknown error')
    history.save()
    
    root = outfile_path.rsplit('.',1)[0]

    tools.qsub(root+'_jobscript')
    return True
    
def clean_resub(outfile_path):
    #Resubmits a job with default parameters, useful for undoing level shift or hfx alterations
    save_scr(outfile_path)
    save_run(outfile_path)
    history = resub_history()
    history.read(outfile_path)
    history.resub_number += 1
    history.status = 'Normal'
    history.notes.append('Needs clean resub')
    history.needs_resub = False
    history.save()
    
    root = outfile_path.rsplit('.',1)[0]
    name = os.path.split(root)[-1]
    directory = os.path.split(outfile_path)[0]
    charge,spinmult,solvent,run_type,levelshifta,levelshiftb,method,hfx,basis,criteria = tools.read_infile(outfile_path)
    
    home = os.getcwd()
    if len(directory) > 0: #if the string is blank, then we're already in the correct directory
        os.chdir(directory)
        
    if os.path.isfile('inscr/optimized.xyz'):
        coordinates = 'inscr/optimized.xyz' #Should trigger for optimization runs
    elif os.path.isfile(name+'.xyz'):
        coordinates = name+'.xyz' #Should trigger for single point runs
    else:
        raise ValueError('No coordinates idenfied for clean in resubmission in directory '+os.getcwd())
    
    configure_dict = read_configure('in_place',outfile_path)
    
    if spinmult == 1:
        tools.write_input(name=name,charge=charge,spinmult=spinmult,solvent = solvent,run_type = run_type, 
                          guess = 'inscr/c0', alternate_coordinates = coordinates,
                          thresholds = criteria, basis = basis, method = configure_dict['method'],
                          levela = configure_dict['levela'], levelb = configure_dict['levelb'])
                          
    else:
        tools.write_input(name=name,charge=charge,spinmult=spinmult,solvent = solvent,run_type = run_type, 
                          guess = 'inscr/ca0 inscr/cb0', alternate_coordinates =coordinates,
                          thresholds = criteria, basis = basis, method = configure_dict['method'],
                          levela = configure_dict['levela'], levelb = configure_dict['levelb'])
    tools.write_jobscript(name,custom_line = '# -fin inscr/')
    os.chdir(home)
    tools.qsub(root+'_jobscript')
    return True
    
def resub_spin(outfile_path):
    #resubmits a spin contaminated job with blyp to help convergence to a non-spin contaminated solution
    history = resub_history()
    history.read(outfile_path)
    resubbed_before = False
    if 'Spin contaminated, lowering HFX to aid convergence' in history.notes:
        resubbed_before = True
        
    if not resubbed_before:
        save_scr(outfile_path)
        save_run(outfile_path)
        history.resub_number += 1
        history.status = 'HFX altered to assist convergence'
        history.needs_resub = True
        history.notes.append('Spin contaminated, lowering HFX to aid convergence')
        history.save()
        
        
        root = outfile_path.rsplit('.',1)[0]
        name = os.path.split(root)[-1]
        directory = os.path.split(outfile_path)[0]
        charge,spinmult,solvent,run_type,levelshifta,levelshiftb,method,hfx,basis,criteria = tools.read_infile(outfile_path)
        
        home = os.getcwd()
        if len(directory) > 0: #if the string is blank, then we're already in the correct directory
            os.chdir(directory)
        tools.write_input(name=name,charge=charge,spinmult=spinmult,solvent = solvent,run_type = run_type, method = 'blyp',
                          thresholds = criteria, basis = basis, levela = levelshifta, levelb = levelshiftb)
        
        tools.write_jobscript(name)
        os.chdir(home)
        tools.qsub(root+'_jobscript')
        return True
        
    else:
        history.status = os.path.split(outfile_path)[-1]+' has been submitted with lower HFX and still converges to a spin contaminated solution'
        history.save()
        return False
        
def resub_scf(outfile_path):
    #Resubmits a job that's having trouble converging the scf with different level shifts (1.0 and 0.1)
    history = resub_history()
    history.read(outfile_path)
    resubbed_before = False
    if 'SCF convergence error, level shifts adjusted to aid convergence' in history.notes:
        resubbed_before = True
        
    if not resubbed_before:
        save_scr(outfile_path)
        save_run(outfile_path)
        history.resub_number += 1
        history.status = 'Level shifts adjusted to assist convergence'
        history.needs_resub = True
        history.notes.append('SCF convergence error, level shifts adjusted to aid convergence')
        history.save()
        
        root = outfile_path.rsplit('.',1)[0]
        name = os.path.split(root)[-1]
        directory = os.path.split(outfile_path)[0]
        charge,spin,solvent,run_type,levelshifta,levelshiftb,method,hfx,basis,criteria = tools.read_infile(outfile_path)
        
        home = os.getcwd()
        if len(directory) > 0: #if the string is blank, then we're already in the correct directory
            os.chdir(directory)
        tools.write_input(name=name,charge=charge,spinmult=spinmult,solvent = solvent,run_type = run_type,
                          levela = 1.0, levelb = 0.1, method = method, thresholds = criteria, hfx = hfx, basis = basis)
                          
        tools.write_jobscript(name)
        os.chdir(home)
        tools.qsub(root+'_jobscript')
        return True
        
    else:
        history.status = os.path.split(outfile_path)[-1]+' has been submitted with levels shifted and is still encountering an scf error'
        history.save()
        return False

def resub_tighter(outfile_path):
    #Takes the path to the outfile of a thermo job with the gradient error problem
    #Finds the parent job and resubmits it with a tighter scf convergence criteria
    
    name = os.path.split(outfile_path)[-1].rsplit('.',1)[0]
    parent_name = name.rsplit('_',1)[0]
    parent_directory = os.path.split(os.path.split(outfile_path)[0])[0]
    parent_path = os.path.join(parent_directory,parent_name+'.out')
    ultratight_path = os.path.join(parent_directory,parent_name+'_ultratight',parent_name+'_ultratight.out')
    
    if os.path.exists(ultratight_path): #This ultratight resubmission has happend before, need to archive the results
        save_scr(ultratight_path, rewrite_inscr = False)
        save_run(ultratight_path)
        history = resub_history()
        history.read(ultratight_path)
        history.resub_number += 1
        history.status = 'Running with tightened convergence thresholds'
        history.needs_resub = False
        history.notes.append('Further tightening convergence thresholds')
        history.save()
        
    jobscript = tools.prep_ultratight(parent_path) #Prep tighter convergence run
    tools.qsub(jobscript) #Submit tighter convergence run
    
    #Set the original thermo run to wait for the ultratight run to finish
    history = resub_history()
    history.read(outfile_path)
    history.waiting = ultratight_path
    history.save()
    
    return True
    
def resub_thermo(outfile_path):
    #Similar to simple resub, but specific for addressing thermo gradient errors
    #Checks for the existance of an ultratight version of this run. If it exists, uses the most up to date version for the new thermo run
    
    save_scr(outfile_path, rewrite_inscr = False)
    save_run(outfile_path)
    history = resub_history()
    history.read(outfile_path)
    history.resub_number += 1
    history.status = 'Normal'
    history.notes.append('Resubmitting thermo, possibly with a better initial geo')
    history.needs_resub = False
    history.save()
    
    
    name = os.path.split(outfile_path)[-1]
    name = name.rsplit('.',1)[0]
    directory = os.path.split(outfile_path)[0]
    parent_name = name.rsplit('_',1)[0]
    parent_directory = os.path.split(os.path.split(outfile_path)[0])[0]
    ultratight_dir = os.path.join(parent_directory,parent_name+'_ultratight')
    
    charge,spinmult,solvent,run_type,levelshifta,levelshiftb,method,hfx,basis,criteria = tools.read_infile(outfile_path)
    
    if os.path.exists(ultratight_dir):
        if os.path.exists(os.path.join(ultratight_dir,'scr','optim.xyz')):
            tools.extract_optimized_geo(os.path.join(ultratight_dir,'scr','optim.xyz'))
            shutil.copy(os.path.join(ultratight_dir,'scr','optimized.xyz'),outfile_path.rsplit('.',1)[0]+'.xyz')
        else:
            raise Exception('Unable to identify the ultratight geometry for run: '+outfile_path)
            
        if spinmult == 1 and os.path.exists(os.path.join(ultratight_dir,'scr','c0')):
            shutil.copy(os.path.join(ultratight_dir,'scr','c0'),os.path.join(directory,'c0'))
        elif spinmult != 1 and os.path.exists(os.path.join(ultratight_dir,'scr','ca0')) and os.path.exists(os.path.join(ultratight_dir,'scr','cb0')):
            shutil.copy(os.path.join(ultratight_dir,'scr','ca0'),os.path.join(directory,'ca0'))
            shutil.copy(os.path.join(ultratight_dir,'scr','cb0'),os.path.join(directory,'cb0'))
        else:
            raise Exception('Unable to find wavefunction files for ultratight geometry for run: '+outfile_path)
    else:
        raise Exception('An ultratight run does not exist for this thermo file. Consider calling simple_resub() or resub_tighter() instead of resub_thermo()')
        
    jobscript = outfile_path.rsplit('.',1)[0]+'_jobscript'
    tools.qsub(jobscript)
    return True
    

def prep_derivative_jobs(directory,list_of_outfiles):
    
    def check_original(job):
        #Returns true if this job is an original job
        #Returns false if this job is already a derivative job
        
        name = os.path.split(job)[-1]
        name = name.rsplit('.',1)[0]
        name = name.split('_')
        
        if 'solventSP' in name or 'vertEA' in name or 'vertIP' in name or 'thermo' in name or 'kp' in name or 'rm' in name or 'ultratight' in name or 'HFXresampling' in name:
            return False
        else:
            return True
        
    
    jobs = filter(check_original,list_of_outfiles)
    
    for job in jobs:
        results = moltools.read_run(job)
        if not results['Is_Oct']:
            print job+' Does not appear to be octahedral! Not generating derivative jobs...'
            continue
        
        configure_dict = read_configure(directory,job)
        
        if configure_dict['solvent']:
            tools.prep_solvent_sp(job)
        if configure_dict['vertEA']:
            tools.prep_vertical_ea(job)
        if configure_dict['vertIP']:
            tools.prep_vertical_ip(job)
        if configure_dict['thermo']:
            tools.prep_thermo(job)
        if configure_dict['hfx_resample']:
            tools.prep_hfx_resample(job)
        if configure_dict['dissociation']:
            moltools.prep_ligand_breakown(job)
        
def reset(outfile_path):
    #Returns the run to the state it was after a given run
    #reference runs as integers, starting from zero
    
    pickle_path = outfile_path.rsplit('.',1)[0]+'.pickle'
    if os.path.isfile(pickle_path):
        
        print 'Resetting run: '+os.path.split(outfile_path)[-1].rsplit('.',1)[0]
        old_path = os.path.join(os.path.split(outfile_path)[0],'pre_reset')
        if not os.path.isdir(old_path):
            os.mkdir(old_path)
        
        #remove all files from states after the specified state
        move = []
        identifier = 1
        while True:
            move.extend(glob.glob(os.path.join(os.path.split(outfile_path)[0],'scr_'+str(identifier))))
            identifier += 1
            if len(glob.glob(os.path.join(os.path.split(outfile_path)[0],'scr_'+str(identifier)))) == 0:
                break #break when all scr_? files are found.
                
        shutil.move(outfile_path,outfile_path[:-4]+'.old') #rename old out so it isn't found in .out searches
        move.append(outfile_path[:-4]+'.old')
        scr_path = os.path.join(os.path.split(outfile_path)[0],'scr')
        move.append(scr_path)
        for path in move:
            #move the paths to their new location, Random numbers prevent clashes
            shutil.move(path,os.path.join(old_path,str(np.random.randint(999999999))+'_'+os.path.split(path)[-1]))
            
        
        history = resub_history()
        history.read(pickle_path)
        outfile = history.outfiles[0]
        writer = open(outfile_path,'w')
        for i in outfile:
            writer.write(i)
        writer.close()
        shutil.move(scr_path+'_0',scr_path)
        shutil.move(pickle_path,os.path.join(old_path,str(np.random.randint(999999999))+'_resub_history'))

def load_history(PATH):
    #takes the path to either an outfile or the resub_history pickle
    #returns the resub_history class object
    
    history = resub_history()
    history.read(PATH)
    return history
    
def main():
    while True:
        print('**********************************')
        print("****** Assessing Job Status ******")
        print('**********************************')
        
        configure_dict = read_configure('in place',None)
                
        number_resubmitted = resub(max_jobs = configure_dict['max_jobs'],max_resub = configure_dict['max_resub']) 
        
        print('**********************************')
        print("******** "+str(number_resubmitted)+" Jobs Submitted ********")
        print('**********************************')
        print 'sleeping for: '+str(configure_dict['sleep'])
	sys.stdout.flush()
       
        time.sleep(configure_dict['sleep']) #sleep for time specified in configure. If not specified, default to 7200 seconds (2 hours)
        
        #Terminate the script if it is no longer submitting jobs
        completeness_status,_ = tools.check_completeness()
        active_jobs = completeness_status['Active']
        if len(active_jobs) == 0 and number_resubmitted ==0:
            break
        
    print('**********************************')
    print("****** Normal Terminatation ******")
    print('**********************************')
    
if __name__ == '__main__':
    main()
        

