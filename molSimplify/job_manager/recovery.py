#!/usr/bin/env python
import os
import shutil
import molSimplify.job_manager.tools as tools
import molSimplify.job_manager.moltools as moltools
from molSimplify.job_manager.classes import resub_history


def simple_resub(outfile_path):
    #Resubmits a job without changing parameters. Particularly useful for CUDA errors.
    save_run(outfile_path, rewrite_inscr = False)
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
    infile_dict = tools.read_infile(outfile_path)
    
    home = os.getcwd()
    if len(directory) > 0: #if the string is blank, then we're already in the correct directory
        os.chdir(directory)
        
    if os.path.isfile('inscr/optimized.xyz'):
        coordinates = 'inscr/optimized.xyz' #Should trigger for optimization runs
    elif os.path.isfile(name+'.xyz'):
        coordinates = name+'.xyz' #Should trigger for single point runs
    else:
        raise ValueError('No coordinates idenfied for clean in resubmission in directory '+os.getcwd())
    
    configure_dict = tools.read_configure('in_place',outfile_path)
    
    infile_dict['coordinates']=coordinates
    infile_dict['method'] = configure_dict['method']
    infile_dict['levelshifta'],infile_dict['levelshiftb'] = configure_dict['levela'],configure_dict['levelb']
    infile_dict['dispersion'] = configure_dict['dispersion']
    infile_dict['convergence_threshold'] = False

    if spinmult == 1:
        infile_dict['guess'] = 'inscr/c0'
        tools.write_input(infile_dict)  
    else:
        infile_dict['guess'] = 'inscr/ca0 inscr/cb0'
        tools.write_input(infile_dict)

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
        history.status = os.path.split(outfile_path)[-1]+' has been submitted with lower HFX and still converges to a spin contaminated solution'
        history.save()
    if 'Needs clean resub' in history.notes:
        resubbed_before = True
        history.status = os.path.split(outfile_path)[-1]+' job recovery has failed - requesting resub_spin() after clean resubmission round'
        history.save()
    if 'HFXresampling' in outfile_path:
        resubbed_before = True
        history.status = os.path.split(outfile_path)[-1]+' is spin contaminated, but submitting with lower HFX does not make sense for HFX resampling jobs'
        history.save()    

    if not resubbed_before:
        save_run(outfile_path, rewrite_inscr = False)
        history = resub_history()
        history.read(outfile_path)
        history.resub_number += 1
        history.status = 'HFX altered to assist convergence'
        history.needs_resub = True
        history.notes.append('Spin contaminated, lowering HFX to aid convergence')
        history.save()
        
        
        root = outfile_path.rsplit('.',1)[0]
        name = os.path.split(root)[-1]
        directory = os.path.split(outfile_path)[0]
        infile_dict = tools.read_infile(outfile_path)
        
        home = os.getcwd()
        if len(directory) > 0: #if the string is blank, then we're already in the correct directory
            os.chdir(directory)

        infile_dict['method'] = 'blyp'
        tools.write_input(infile_dict)
        
        tools.write_jobscript(name)
        os.chdir(home)
        tools.qsub(root+'_jobscript')
        return True
        
    else:
        return False
        
def resub_scf(outfile_path):
    #Resubmits a job that's having trouble converging the scf with different level shifts (1.0 and 0.1)
    history = resub_history()
    history.read(outfile_path)
    resubbed_before = False
    if 'SCF convergence error, level shifts adjusted to aid convergence' in history.notes:
        resubbed_before = True
        history.status = os.path.split(outfile_path)[-1]+' has been submitted with levels shifted and is still encountering an scf error'
        history.save()
    if 'Needs clean resub' in history.notes:
        resubbed_before = True
        history.status = os.path.split(outfile_path)[-1]+' job recovery has failed - requesting resub_scf() after clean resubmission round'
        history.save()
        
    if not resubbed_before:
        save_run(outfile_path, rewrite_inscr = False)
        history = resub_history()
        history.read(outfile_path)
        history.resub_number += 1
        history.status = 'Level shifts adjusted to assist convergence'
        history.needs_resub = True
        history.notes.append('SCF convergence error, level shifts adjusted to aid convergence')
        history.save()
        
        root = outfile_path.rsplit('.',1)[0]
        name = os.path.split(root)[-1]
        directory = os.path.split(outfile_path)[0]
        infile_dict = tools.read_infile(outfile_path)
        
        home = os.getcwd()
        if len(directory) > 0: #if the string is blank, then we're already in the correct directory
            os.chdir(directory)
        infile_dict['levelshifta'],infile_dict['levelshiftb'] = 1.0,0,1
        tools.write_input(infile_dict)
                          
        tools.write_jobscript(name)
        os.chdir(home)
        tools.qsub(root+'_jobscript')
        return True
        
    else:
        return False
        
def resub_bad_geo(outfile_path,home_directory):
    #Resubmits a job that's converged to a bad geometry with additional contraints
    history = resub_history()
    history.read(outfile_path)
    resubbed_before = False
    if 'Bad geometry detected, adding constraints and trying again' in history.notes:
        resubbed_before = True
        history.status = os.path.split(outfile_path)[-1]+" has been submitted with additional constraints and still isn't a good geometry"
        history.save()
    if 'Needs clean resub' in history.notes:
        resubbed_before = True
        history.status = os.path.split(outfile_path)[-1]+' job recovery has failed - requesting resub_bad_geo after clean resubmission round'
        history.save()
        
    if not resubbed_before:
        save_run(outfile_path, rewrite_inscr = True)
        history = resub_history()
        history.read(outfile_path)
        history.resub_number += 1
        history.status = 'Constraints added to help convergence'
        history.needs_resub = True
        history.notes.append('Bad geometry detected, adding constraints and trying again')
        history.save()
        
        root = outfile_path.rsplit('.',1)[0]
        name = os.path.split(root)[-1]
        directory = os.path.split(outfile_path)[0]
        infile_dict = tools.read_infile(outfile_path)
        
        if infile_dict['constraints']:
            raise Exception('resub.py does not currently support the use of external atom constraints. These will be overwritten by clean_resub() during job recovery')
        
        goal_geo = tools.read_configure(home_directory,outfile_path)['geo_check']
        if not goal_geo:
            raise Exception('Goal geometry not specified, job '+outfile_path+' should not have been labelled bad geo!')
        else:
            metal_index,bonded_atom_indices = moltools.get_metal_and_bonded_atoms(outfile_path,goal_geo)
            #convert indexes from zero-indexed to one-indexed
            metal_index += 1
            bonded_atom_indices = [index + 1 for index in bonded_atom_indices]
            #Convert to TeraChem input syntax
            constraints = ['bond '+str(metal_index)+'_'+str(index)+'\n' for index in bonded_atom_indices]
        
        home = os.getcwd()
        if len(directory) > 0: #if the string is blank, then we're already in the correct directory
            os.chdir(directory)

        infile_dict['constraints'] = constraints
        tools.write_input(infile_dict)
                          
        tools.write_jobscript(name)
        os.chdir(home)
        tools.qsub(root+'_jobscript')
        return True
        
    else:
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
        save_run(ultratight_path, rewrite_inscr = False)
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
    
    save_run(outfile_path, rewrite_inscr = False)
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
    
    infile_dict = tools.read_infile(outfile_path)
    
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