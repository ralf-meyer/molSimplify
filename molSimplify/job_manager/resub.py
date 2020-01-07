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


def kill_jobs(kill_names, message1='Killing job: ', message2=' early'):
    # This function takes a list of job names and kills the jobs associated with them, if the jobs are active
    if type(kill_names) != list:
        kill_names = [kill_names]

    active_jobs, active_ids = tools.list_active_jobs(ids=True)
    active_jobs = zip(active_jobs, active_ids)

    jobs_to_kill = [[name, id_] for name, id_ in active_jobs if name in kill_names]

    for name, id_ in jobs_to_kill:
        print message1 + name + message2
        tools.call_bash('qdel ' + str(id_))


def prep_derivative_jobs(directory, list_of_outfiles):
    for job in list_of_outfiles:
        configure_dict = tools.read_configure(directory, job)

        if configure_dict['solvent']:
            tools.prep_solvent_sp(job, configure_dict['solvent'])
        if configure_dict['functionalsSP']:
            tools.prep_functionals_sp(job, configure_dict['functionalsSP'])
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


def resub(directory='in place'):
    # Takes a directory, resubmits errors, scf failures, and spin contaminated cases

    configure_dict = tools.read_configure(directory, None)
    print 'Global Configure File Found:'
    print configure_dict

    max_resub = configure_dict['max_resub']
    max_jobs = configure_dict['max_jobs']

    # Get the state of all jobs being managed by this instance of the job manager
    completeness = moltools.check_completeness(directory, max_resub)
    errors = completeness['Error']  # These are calculations which failed to complete
    scf_errors = completeness[
        'SCF_Error']  # These are calculations which failed to complete, appear to have an scf error, and hit wall time
    need_resub = completeness['Resub']  # These are calculations with level shifts changed or hfx exchange changed
    spin_contaminated = completeness['Spin_contaminated']  # These are finished jobs with spin contaminated solutions
    active = completeness['Active']  # These are jobs which are currently running
    thermo_grad_error = completeness['Thermo_grad_error']  # These are thermo jobs encountering the thermo grad error
    waiting = completeness[
        'Waiting']  # These are jobs which are or were waiting for another job to finish before continuing.
    bad_geos = completeness['Bad_geos']  # These are jobs which finished, but converged to a bad geometry.
    finished = completeness['Finished']
    nactive = tools.get_number_active() #number of active jobs, counting bundled jobs as a single job

    # Kill SCF errors in progress, which are wasting computational resources
    all_scf_errors = completeness['SCF_Errors_Including_Active']  # These are all jobs which appear to have scf error, including active ones
    scf_errors_to_kill = [scf_err for scf_err in all_scf_errors if scf_err not in scf_errors]
    names_to_kill = [os.path.split(scf_err)[-1].rsplit('.', 1)[0] for scf_err in scf_errors_to_kill]
    kill_jobs(names_to_kill, message1='Job: ', message2=' appears to have an scf error. Killing this job early')

    # Prep derivative jobs such as thermo single points, vertical IP, and ligand dissociation energies
    needs_derivative_jobs = filter(tools.check_original, finished)
    prep_derivative_jobs(directory, needs_derivative_jobs)

    resubmitted = []  # Resubmitted list gets True if the job is submitted or False if not. Contains booleans, not job identifiers.

    # Resub unidentified errors
    for error in errors:
        if nactive + np.sum(resubmitted) >= max_jobs:
            continue
        resub_tmp = recovery.simple_resub(error)
        if resub_tmp:
            print('Unidentified error in job: ' + os.path.split(error)[-1] + ' -Resubmitting')
            print('')
        resubmitted.append(resub_tmp)

    # Resub scf convergence errors
    for error in scf_errors:
        if nactive + np.sum(resubmitted) >= max_jobs:
            continue
        local_configure = tools.read_configure(directory, None)
        if 'scf' in local_configure['job_recovery']:
            resub_tmp = recovery.resub_scf(error)
            if resub_tmp:
                print('SCF error identified in job: ' + os.path.split(error)[
                    -1] + ' -Resubmitting with adjusted levelshifts')
                print('')
            resubmitted.append(resub_tmp)

    # Resub jobs which converged to bad geometries with additional constraints
    for error in bad_geos:
        if nactive + np.sum(resubmitted) >= max_jobs:
            continue
        local_configure = tools.read_configure(directory, None)
        if 'bad_geo' in local_configure['job_recovery']:
            resub_tmp = recovery.resub_bad_geo(error, directory)
            if resub_tmp:
                print('Bad final geometry in job: ' + os.path.split(error)[
                    -1] + ' -Resubmitting from initial structure with additional constraints')
                print('')
            resubmitted.append(resub_tmp)

    # Resub spin contaminated cases
    for error in spin_contaminated:
        if nactive + np.sum(resubmitted) >= max_jobs:
            continue
        local_configure = tools.read_configure(directory, None)
        if 'spin_contaminated' in local_configure['job_recovery']:
            resub_tmp = recovery.resub_spin(error)
            if resub_tmp:
                print('Spin contamination identified in job: ' + os.path.split(error)[
                    -1] + ' -Resubmitting with adjusted HFX')
                print('')
            resubmitted.append(resub_tmp)

    # Resub jobs with atypical parameters used to aid convergence
    for error in need_resub:
        if nactive + np.sum(resubmitted) >= max_jobs:
            continue
        resub_tmp = recovery.clean_resub(error)
        if resub_tmp:
            print('Job ' + os.path.split(error)[-1] + ' needs to be rerun with typical paramters. -Resubmitting')
            print('')
        resubmitted.append(resub_tmp)

    # Create a job with a tighter convergence threshold for failed thermo jobs
    for error in thermo_grad_error:
        if nactive + np.sum(resubmitted) >= max_jobs:
            continue
        local_configure = tools.read_configure(directory, None)
        if 'thermo_grad_error' in local_configure['job_recovery']:
            resub_tmp = recovery.resub_tighter(error)
            if resub_tmp:
                print('Job ' + os.path.split(error)[
                    -1] + ' needs a better initial geo. Creating a geometry run with tighter convergence criteria')
                print('')
            resubmitted.append(resub_tmp)

    # Look at jobs in "waiting," resume them if the job they were waiting for is finished
    # Currently, this should only ever be thermo jobs waiting for an ultratight job
    for waiting_dict in waiting:
        if nactive + np.sum(resubmitted) >= max_jobs:
            continue
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
                raise Exception('A method for resuming job: ' + job + ' is not defined')
        else:
            resubmitted.append(False)

    # Submit jobs which haven't yet been submitted
    to_submit = []
    jobscripts = tools.find('*_jobscript')
    active_jobs = tools.list_active_jobs(home_directory=directory,parse_bundles=True)
    for job in jobscripts:
        if not os.path.isfile(job.rsplit('_', 1)[0] + '.out') and not os.path.split(job.rsplit('_', 1)[0])[
                                                                          -1] in active_jobs:
            to_submit.append(job)

    short_jobs_to_submit = [i for i in to_submit if tools.check_short_single_point(i)]
    long_jobs_to_submit = [i for i in to_submit if i not in short_jobs_to_submit]
    if len(short_jobs_to_submit) > 0:
        bundled_jobscripts = tools.bundle_jobscripts(os.getcwd(),short_jobs_to_submit)
    else:
        bundled_jobscripts = []
    to_submit = long_jobs_to_submit + bundled_jobscripts

    submitted = []
    for job in to_submit:
        if len(submitted) + nactive + np.sum(resubmitted) >= max_jobs:
            continue
        print('Initial sumbission for job: ' + os.path.split(job)[-1])
        tools.qsub(job)
        submitted.append(True)

    number_resubmitted = np.sum(np.array(resubmitted + submitted))
    # ~ print str(number_resubmitted)+' Jobs submitted'
    return int(number_resubmitted), int(len(completeness['Active']))


def main():
    counter = 0
    while True:
        print('**********************************')
        print("****** Assessing Job Status ******")
        print('**********************************')
        time1 = time.time()
        fil = open('complete', 'w')
        fil.write('Active')
        fil.close()

        number_resubmitted, number_active = resub()

        print('**********************************')
        print("******** " + str(number_resubmitted) + " Jobs Submitted ********")
        print('**********************************')

        print('job cycle took: ' + str(time.time() - time1))
        configure_dict = tools.read_configure('in place', None)
        print('sleeping for: ' + str(configure_dict['sleep']))
        sys.stdout.flush()
        time.sleep(configure_dict[
                       'sleep'])  # sleep for time specified in configure. If not specified, default to 7200 seconds (2 hours)

        # Terminate the script if it is no longer submitting jobs
        if number_resubmitted == 0 and number_active == 0:
            counter += 1
        else:
            counter = 0
        if counter >= 3:
            break

    print('**********************************')
    print("****** Normal Terminatation ******")
    print('**********************************')
    fil = open('complete', 'w')
    fil.write('True')
    fil.close()


if __name__ == '__main__':
    main()
