#!/usr/bin/env python

import os
import glob
import numpy as np
import shutil
import time
import sys
import molSimplify.job_manager.tools as tools
import molSimplify.job_manager.moltools as moltools
import molSimplify.job_manager.recovery as recovery
from molSimplify.job_manager.classes import resub_history

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

## Save the outfile within the resub_history pickel object
def save_run(outfile_path, rewrite_inscr = True):
    
    save_scr(outfile_path, rewrite_inscr = rewrite_inscr)
    
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

def kill_jobs(kill_names,message1='Killing job: ',message2=' early'):
    # This function takes a list of job names and kills the jobs associated with them, if the jobs are active
    if type(kill_names) != list:
        kill_names = [kill_names]

    active_jobs,active_ids = tools.list_active_jobs(ids=True)
    active_jobs = zip(active_jobs,active_ids)

    jobs_to_kill = [[name,id_] for name,id_ in active_jobs if name in kill_names]

    for name,id_ in jobs_to_kill:
        print message1+name+message2
        tools.call_bash('qdel '+str(id_))

def resub(directory = 'in place'):
    #Takes a directory, resubmits errors, scf failures, and spin contaminated cases
    
    configure_dict = tools.read_configure(directory,None)
    print 'Global Configure File Found:'
    print configure_dict

    max_resub = configure_dict['max_resub']
    max_jobs = configure_dict['max_jobs']
    
    #Get the state of all jobs being managed by this instance of the job manager
    completeness = moltools.check_completeness(directory,max_resub)
    errors = completeness['Error'] #These are calculations which failed to complete
    scf_errors = completeness['SCF_Error'] #These are calculations which failed to complete, appear to have an scf error, and hit wall time
    need_resub = completeness['Resub'] #These are calculations with level shifts changed or hfx exchange changed
    spin_contaminated = completeness['Spin_contaminated'] #These are finished jobs with spin contaminated solutions
    active = completeness['Active'] #These are jobs which are currently running
    thermo_grad_error = completeness['Thermo_grad_error'] #These are thermo jobs encountering the thermo grad error
    waiting = completeness['Waiting'] #These are jobs which are or were waiting for another job to finish before continuing.
    bad_geos = completeness['Bad_geos'] #These are jobs which finished, but converged to a bad geometry.
    finished = completeness['Finished']

    #Kill SCF errors in progress, which are wasting computational resources
    all_scf_errors = completeness['SCF_Errors_Including_Active'] #These are all jobs which appear to have scf error, including active ones
    scf_errors_to_kill = [scf_err for scf_err in all_scf_errors if scf_err not in scf_errors]
    names_to_kill = [os.path.split(scf_err)[-1].rsplit('.',1)[0] for scf_err in scf_errors_to_kill]
    kill_jobs(names_to_kill,message1='Job: ',message2=' appears to have an scf error. Killing this job early')
    
    #Prep derivative jobs such as thermo single points, vertical IP, and ligand dissociation energies
    needs_derivative_jobs = filter(tools.check_original,finished)
    prep_derivative_jobs(directory,needs_derivative_jobs)
    
    resubmitted = [] #Resubmitted list gets True if the job is submitted or False if not. Contains booleans, not job identifiers.

    #Resub unidentified errors
    for error in errors:
        if len(active)+np.sum(resubmitted) >= max_jobs:
            continue
        resub_tmp = recovery.simple_resub(error)
        if resub_tmp:
            print('Unidentified error in job: '+os.path.split(error)[-1]+' -Resubmitting')
            print('')
        resubmitted.append(resub_tmp)

    #Resub scf convergence errors
    for error in scf_errors:
        if len(active)+np.sum(resubmitted) >= max_jobs:
            continue
        local_configure = tools.read_configure(directory,None)
        if 'scf' in local_configure['job_recovery']:
            resub_tmp = recovery.resub_scf(error)
            if resub_tmp:
                print('SCF error identified in job: '+os.path.split(error)[-1]+' -Resubmitting with adjusted levelshifts')
                print('')
            resubmitted.append(resub_tmp)
    
    #Resub jobs which converged to bad geometries with additional constraints
    for error in bad_geos:
        if len(active)+np.sum(resubmitted) >= max_jobs:
            continue
        local_configure = tools.read_configure(directory,None)
        if 'bad_geo' in local_configure['job_recovery']:
            resub_tmp = recovery.resub_bad_geo(error,directory)
            if resub_tmp:
                print('Bad final geometry in job: '+os.path.split(error)[-1]+' -Resubmitting from initial structure with additional constraints')
                print('')
            resubmitted.append(resub_tmp)
        
    #Resub spin contaminated cases
    for error in spin_contaminated:
        if len(active)+np.sum(resubmitted) >= max_jobs:
            continue
        local_configure = tools.read_configure(directory,None)
        if 'spin_contaminated' in local_configure['job_recovery']:
            resub_tmp = recovery.resub_spin(error)
            if resub_tmp:
                print('Spin contamination identified in job: '+os.path.split(error)[-1]+' -Resubmitting with adjusted HFX')
                print('')
            resubmitted.append(resub_tmp)
        
    #Resub jobs with atypical parameters used to aid convergence
    for error in need_resub:
        if len(active)+np.sum(resubmitted) >= max_jobs:
            continue
        resub_tmp = recovery.clean_resub(error)
        if resub_tmp:
            print('Job '+os.path.split(error)[-1]+' needs to be rerun with typical paramters. -Resubmitting')
            print('')
        resubmitted.append(resub_tmp)
        
    #Create a job with a tighter convergence threshold for failed thermo jobs
    for error in thermo_grad_error:
        if len(active)+np.sum(resubmitted) >= max_jobs:
            continue
        local_configure = tools.read_configure(directory,None)
        if 'thermo_grad_error' in local_configure['job_recovery']:
            resub_tmp = recovery.resub_tighter(error)
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
                resubmitted.append(recovery.resub_thermo(job))
            else:
                raise Exception('A method for resuming job: '+job+' is not defined')
        else:
            resubmitted.append(False)
        
    #Submit jobs which haven't yet been submitted
    to_submit = []
    jobscripts = tools.find('*_jobscript')
    active_jobs = tools.list_active_jobs()
    for job in jobscripts:
        if not os.path.isfile(job.rsplit('_',1)[0]+'.out') and not os.path.split(job.rsplit('_',1)[0])[-1] in active_jobs:
            to_submit.append(job)
    submitted = []
    for job in to_submit:
        if len(submitted)+len(active)+np.sum(resubmitted) >= max_jobs:
            continue
        print('Initial sumbission for job: '+os.path.split(job)[-1])
        tools.qsub(job)
        submitted.append(True)
    
    number_resubmitted = np.sum(np.array(resubmitted+submitted))
    # ~ print str(number_resubmitted)+' Jobs submitted'
    return int(number_resubmitted)
    

def prep_derivative_jobs(directory,list_of_outfiles):
    
    for job in list_of_outfiles:
        configure_dict = tools.read_configure(directory,job)

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
    
    pickle_path = outfile_path.rsplit('.',1)[0]+'.pickle'
    if os.path.isfile(pickle_path):
        
        print 'Resetting run: '+os.path.split(outfile_path)[-1].rsplit('.',1)[0]
        old_path = os.path.join(os.path.split(outfile_path)[0],'pre_reset')
        if not os.path.isdir(old_path):
            os.mkdir(old_path)
        
        #Find all the stdout and stderr files related to previous runs.
        queue_output = glob.glob(outfile_path.rsplit('.',1)[0]+'.e*')
        queue_output.extend(glob.glob(outfile_path.rsplit('.',1)[0]+'.pe*'))
        queue_output.extend(glob.glob(outfile_path.rsplit('.',1)[0]+'.po*'))
        queue_output.extend(glob.glob(outfile_path.rsplit('.',1)[0]+'.o*'))
        queue_output = [i for i in queue_output if i[-1] in ['1','2','3','4','5','6','7','8','9','0']]
        
        #remove all files from states after the specified state
        move = []
        identifier = 1
        while True:
            move.extend(glob.glob(os.path.join(os.path.split(outfile_path)[0],'scr_'+str(identifier))))
            identifier += 1
            if len(glob.glob(os.path.join(os.path.split(outfile_path)[0],'scr_'+str(identifier)))) == 0:
                break #break when all scr_? files are found.
                
        shutil.move(outfile_path,outfile_path[:-4]+'.old') #rename old out so it isn't found in .out searches
        shutil.move(outfile_path[:-4]+'_jobscript',outfile_path[:-4]+'_oldjob') #rename old jobscript so it isn't thought to be  job that hasn't started yet
        move.append(outfile_path[:-4]+'.old')
        move.append(outfile_path[:-4]+'.xyz')
        move.append(outfile_path[:-4]+'.in')
        move.append(outfile_path[:-4]+'_oldjob')
        if os.path.isdir(os.path.join(os.path.split(outfile_path)[0],'inscr')):
            move.append(os.path.join(os.path.split(outfile_path)[0],'inscr'))
        move.extend(queue_output)
        scr_path = os.path.join(os.path.split(outfile_path)[0],'scr')
        move.append(scr_path)
        for path in move:
            #move the paths to their new location, Random numbers prevent clashes
            try:
                shutil.move(path,os.path.join(old_path,str(np.random.randint(999999999))+'_'+os.path.split(path)[-1]))
            except:
                print 'No file found for: '+path
            
        
        #Rewrite the .xyz, .in, jobscript, and .out file to be the same as they were after the first run
        history = resub_history()
        history.read(pickle_path)
        outfile = history.outfiles[0]
        infile = history.infiles[0]
        jobscript = history.jobscripts[0]
        xyz = history.xyzs[0]
        writer = open(outfile_path,'w')
        for i in outfile:
            writer.write(i)
        writer.close()
        writer = open(outfile_path.rsplit('.',1)[0]+'.in','w')
        for i in infile:
            writer.write(i)
        writer.close()
        writer = open(outfile_path.rsplit('.',1)[0]+'.xyz','w')
        for i in xyz:
            writer.write(i)
        writer.close()
        writer = open(outfile_path.rsplit('.',1)[0]+'_jobscript','w')
        for i in jobscript:
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
        time1 = time.time()
                
        number_resubmitted = resub() 
        
        print('**********************************')
        print("******** "+str(number_resubmitted)+" Jobs Submitted ********")
        print('**********************************')
        
        print('job cycle took: '+str(time.time()-time1))
        configure_dict = tools.read_configure('in place',None)
        print('sleeping for: '+str(configure_dict['sleep']))
        sys.stdout.flush()
        time.sleep(configure_dict['sleep']) #sleep for time specified in configure. If not specified, default to 7200 seconds (2 hours)
        
        #Terminate the script if it is no longer submitting jobs
        completeness_status = tools.check_completeness()
        active_jobs = completeness_status['Active']
        if len(active_jobs) == 0 and number_resubmitted ==0:
            break
        
    print('**********************************')
    print("****** Normal Terminatation ******")
    print('**********************************')
    
if __name__ == '__main__':
    main()
        

