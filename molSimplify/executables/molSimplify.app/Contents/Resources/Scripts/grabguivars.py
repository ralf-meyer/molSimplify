# Written by Tim Ioannidis for HJK Group
# Dpt of Chemical Engineering, MIT

##############################################################
########## This script processes the input file  #############
##############################################################

# import std modules
import glob, os, re, argparse, sys, time
from io import *

###########################################
########## check true or false  ###########
###########################################
def checkTrue(arg):
    if 'n' in arg.lower() or '0' in arg.lower() or 'f' in arg.lower():
        return False
    else:
        return True
        
#########################################
########## collect ligands  #############
#########################################
def getligands(gui):
    ligs,ligoccs,lcats,kHs,MLb,lang,lname = '','','','','','',''
    if len(gui.lig0.text().replace(' ','')) > 0:
        ligs += gui.lig0.text().replace(' ','')
        ligoccs += str(gui.lig0occ.value())
        kHs += gui.lig0H.currentText()
        lcats += gui.lig0conn.text().replace(' ','')
        MLb += gui.lig0ML.text().replace(' ','')
        lang += gui.lig0an.text().replace(' ','')
        lname += gui.lig0nam.text().replace(' ','')
    if len(gui.lig1.text().replace(' ','')) > 0:
        ligs += ','+gui.lig1.text().replace(' ','')
        ligoccs += ','+str(gui.lig1occ.value())
        kHs += ','+gui.lig1H.currentText()
        lcats += '/'+gui.lig1conn.text().replace(' ','')
        MLb += ','+gui.lig1ML.text().replace(' ','')
        lang += ','+gui.lig1an.text().replace(' ','')
        lname += ','+gui.lig1nam.text().replace(' ','')
    if len(gui.lig2.text().replace(' ','')) > 0:
        ligs += ','+gui.lig2.text().replace(' ','')
        ligoccs += ','+str(gui.lig2occ.value())
        kHs += ','+gui.lig2H.currentText()
        lcats += '/'+gui.lig2conn.text().replace(' ','')
        MLb += ','+gui.lig2ML.text().replace(' ','')
        lang += ','+gui.lig2an.text().replace(' ','')
        lname += ','+gui.lig2nam.text().replace(' ','')
    if len(gui.lig3.text().replace(' ','')) > 0:
        ligs += ','+gui.lig3.text().replace(' ','')
        ligoccs += ','+str(gui.lig3occ.value())
        kHs += ','+gui.lig3H.currentText()
        lcats += '/'+gui.lig3conn.text().replace(' ','')
        MLb += ','+gui.lig3ML.text().replace(' ','')
        lang += ','+gui.lig3an.text().replace(' ','')
        lname += ','+gui.lig3nam.text().replace(' ','')
    if len(gui.lig4.text().replace(' ','')) > 0:
        ligs += ','+gui.lig4.text().replace(' ','')
        ligoccs += ','+str(gui.lig4occ.value())
        kHs += ','+gui.lig4H.currentText()
        lcats += ','+gui.lig4conn.text().replace(' ','')
        MLb += '/'+gui.lig4ML.text().replace(' ','')
        lang += ','+gui.lig4an.text().replace(' ','')
        lname += ','+gui.lig4nam.text().replace(' ','')
    if len(gui.lig5.text().replace(' ','')) > 0:
        ligs += ','+gui.lig5.text().replace(' ','')
        ligoccs += ','+str(gui.lig5occ.value())
        kHs += ','+gui.lig5H.currentText()
        lcats += '/'+gui.lig5conn.text().replace(' ','')
        MLb += ','+gui.lig5ML.text().replace(' ','')
        lang += ','+gui.lig5an.text().replace(' ','')
        lname += ','+gui.lig5nam.text().replace(' ','')
    if len(gui.lig6.text().replace(' ','')) > 0:
        ligs += ','+gui.lig6.text().replace(' ','')
        ligoccs += ','+str(gui.lig6occ.value())
        kHs += ','+gui.lig6H.currentText()
        lcats += '/'+gui.lig6conn.text().replace(' ','')
        MLb += ','+gui.lig6ML.text().replace(' ','')
        lang += ','+gui.lig6an.text().replace(' ','')
        lname += ','+gui.lig6nam.text().replace(' ','')
    if len(gui.lig7.text().replace(' ','')) > 0:
        ligs += ','+gui.lig7.text().replace(' ','')
        ligoccs += ','+str(gui.lig7occ.value())
        kHs += ','+gui.lig7H.currentText()
        lcats += '/'+gui.lig7conn.text().replace(' ','')
        MLb += ','+gui.lig7ML.text().replace(' ','')
        lang += ','+gui.lig7an.text().replace(' ','')
        lname += ','+gui.lig7nam.text().replace(' ','')
    return ligs,ligoccs,lcats,kHs,MLb,lang,lname


#########################################
########## set ligands  #############
#########################################
def setligands(gui,ligs,ligoccs,lcats,kHs,MLb,lang,lname):
    ligs = re.split(',|\t|&',ligs)
    ligoccs = re.split(',|\t|&',ligoccs)
    lcats = re.split('/|\t|&',lcats)
    kHs = re.split(',|\t|&',kHs)
    MLb = re.split(',|\t|&',MLb)
    lang = re.split(',|\t|&',lang)
    lname = re.split(',|\t|&',lname)
    if len(ligs) > 0:
        gui.lig0.setText(ligs[0])
        if len(ligoccs) > 0 and ligoccs[0]!='':
            gui.lig0occ.setValue(int(ligoccs[0]))
        if len(lcats) > 0:
            gui.lig0conn.setText(lcats[0])
        if len(kHs) > 0:
            if checkTrue(kHs[0]):
                gui.lig0H.setCurrentText('yes')
            else:
                gui.lig0H.setCurrentText('no')
        if len(MLb) > 0:
            gui.lig0ML.setText(MLb[0])
        if len(lang) > 0:
            gui.lig0an.setText(lang[0])
        if len(lname) > 0:
            gui.lig0nam.setText(lname[0])
    if len(ligs) > 1:
        gui.lig1.setText(ligs[1])
        if len(ligoccs) > 1 and ligoccs[1]!='':
            gui.lig1occ.setValue(int(ligoccs[1]))
        if len(lcats) > 1:
            gui.lig1conn.setText(lcats[1])
        if len(kHs) > 1:
            if checkTrue(kHs[1]):
                gui.lig1H.setCurrentText('yes')
            else:
                gui.lig1H.setCurrentText('no')
        if len(MLb) > 1:
            gui.lig1ML.setText(MLb[1])
        if len(lang) > 1:
            gui.lig1an.setText(lang[1])
        if len(lname) > 1:
            gui.lig1nam.setText(lname[1])
        gui.lig0add.hide()
        gui.lig1.setDisabled(False)
        gui.lig1occ.setDisabled(False)
        gui.lig1conn.setDisabled(False)
        gui.lig1H.setDisabled(False)
        gui.lig1ML.setDisabled(False)
        gui.lig1an.setDisabled(False)
        gui.lig1nam.setDisabled(False)
        gui.lig1add.show()
        gui.lig1add.setDisabled(False)
    if len(ligs) > 2:
        gui.lig2.setText(ligs[2])
        if len(ligoccs) > 2 and ligoccs[2]!='':
            gui.lig2occ.setValue(int(ligoccs[2]))
        if len(lcats) > 2:
            gui.lig2conn.setText(lcats[2])
        if len(kHs) > 2:
            if checkTrue(kHs[2]):
                gui.lig2H.setCurrentText('yes')
            else:
                gui.lig2H.setCurrentText('no')
        if len(MLb) > 2:
            gui.lig2ML.setText(MLb[2])
        if len(lang) > 2:
            gui.lig2an.setText(lang[2])
        if len(lname) > 2:
            gui.lig2nam.setText(lname[2])
        gui.lig1add.hide()
        gui.lig2.setDisabled(False)
        gui.lig2occ.setDisabled(False)
        gui.lig2conn.setDisabled(False)
        gui.lig2H.setDisabled(False)
        gui.lig2ML.setDisabled(False)
        gui.lig2an.setDisabled(False)
        gui.lig2nam.setDisabled(False)
        gui.lig2add.show()
        gui.lig2add.setDisabled(False)
    if len(ligs) > 3:
        gui.lig3.setText(ligs[3])
        if len(ligoccs) > 3 and ligoccs[3]!='':
            gui.lig3occ.setValue(int(ligoccs[3]))
        if len(lcats) > 3:
            gui.lig3conn.setText(lcats[3])
        if len(kHs) > 3:
            if checkTrue(kHs[3]):
                gui.lig3H.setCurrentText('yes')
            else:
                gui.lig3H.setCurrentText('no')
        if len(MLb) > 3:
            gui.lig3ML.setText(MLb[3])
        if len(lang) > 3:
            gui.lig3an.setText(lang[3])
        if len(lname) > 3:
            gui.lig3nam.setText(lname[3])
        gui.lig2add.hide()
        gui.lig3.setDisabled(False)
        gui.lig3occ.setDisabled(False)
        gui.lig3conn.setDisabled(False)
        gui.lig3H.setDisabled(False)
        gui.lig3ML.setDisabled(False)
        gui.lig3an.setDisabled(False)
        gui.lig3nam.setDisabled(False)
        gui.lig3add.show()
        gui.lig3add.setDisabled(False)
    if len(ligs) > 4:
        gui.lig4.setText(ligs[4])
        if len(ligoccs) > 4 and ligoccs[4]!='':
            gui.lig4occ.setValue(int(ligoccs[4]))
        if len(lcats) > 4:
            gui.lig4conn.setText(lcats[4])
        if len(kHs) > 4:
            if checkTrue(kHs[4]):
                gui.lig4H.setCurrentText('yes')
            else:
                gui.lig4H.setCurrentText('no')
        if len(MLb) > 4:
            gui.lig4ML.setText(MLb[4])
        if len(lang) > 4:
            gui.lig4an.setText(lang[4])
        if len(lname) > 2:
            gui.lig4nam.setText(lname[4])
        gui.lig3add.hide()
        gui.lig4.setDisabled(False)
        gui.lig4occ.setDisabled(False)
        gui.lig4conn.setDisabled(False)
        gui.lig4H.setDisabled(False)
        gui.lig4ML.setDisabled(False)
        gui.lig4an.setDisabled(False)
        gui.lig4nam.setDisabled(False)
        gui.lig4add.show()
        gui.lig4add.setDisabled(False)
    if len(ligs) > 5:
        gui.lig5.setText(ligs[5])
        if len(ligoccs) > 5 and ligoccs[5]!='':
            gui.lig5occ.setValue(int(ligoccs[5]))
        if len(lcats) > 5:
            gui.lig5conn.setText(lcats[5])
        if len(kHs) > 5:
            if checkTrue(kHs[5]):
                gui.lig5H.setCurrentText('yes')
            else:
                gui.lig5H.setCurrentText('no')
        if len(MLb) > 5:
            gui.lig5ML.setText(MLb[5])
        if len(lang) > 5:
            gui.lig5an.setText(lang[5])
        if len(lname) > 5:
            gui.lig5nam.setText(lname[5])
        gui.lig4add.hide()
        gui.lig5.setDisabled(False)
        gui.lig5occ.setDisabled(False)
        gui.lig5conn.setDisabled(False)
        gui.lig5H.setDisabled(False)
        gui.lig5ML.setDisabled(False)
        gui.lig5an.setDisabled(False)
        gui.lig5nam.setDisabled(False)
        gui.lig5add.show()
        gui.lig5add.setDisabled(False)
    if len(ligs) > 6:
        gui.lig6.setText(ligs[6])
        if len(ligoccs) > 6 and ligoccs[6]!='':
            gui.lig6occ.setValue(int(ligoccs[6]))
        if len(lcats) > 6:
            gui.lig6conn.setText(lcats[6])
        if len(kHs) > 6:
            if checkTrue(kHs[6]):
                gui.lig6H.setCurrentText('yes')
            else:
                gui.lig6H.setCurrentText('no')
        if len(MLb) > 6:
            gui.lig6ML.setText(MLb[6])
        if len(lang) > 6:
            gui.lig6an.setText(lang[6])
        if len(lname) > 6:
            gui.lig6nam.setText(lname[6])
        gui.lig5add.hide()
        gui.lig6.setDisabled(False)
        gui.lig6occ.setDisabled(False)
        gui.lig6conn.setDisabled(False)
        gui.lig6H.setDisabled(False)
        gui.lig6ML.setDisabled(False)
        gui.lig6an.setDisabled(False)
        gui.lig6nam.setDisabled(False)
        gui.lig6add.setDisabled(False)
        gui.lig6add.show()
        gui.lig6add.setDisabled(False)
    if len(ligs) > 7:
        gui.lig7.setText(ligs[7])
        if len(ligoccs) > 7 and ligoccs[7]!='':
            gui.lig7occ.setValue(int(ligoccs[7]))
        if len(lcats) > 7:
            gui.lig7conn.setText(lcats[7])
        if len(kHs) > 7:
            if checkTrue(kHs[7]):
                gui.lig7H.setCurrentText('yes')
            else:
                gui.lig7H.setCurrentText('no')
        if len(MLb) > 7:
            gui.lig7ML.setText(MLb[7])
        if len(lang) > 7:
            gui.lig7an.setText(lang[7])
        if len(lname) > 7:
            gui.lig7nam.setText(lname[7])
        gui.lig6add.hide()
        gui.lig7.setDisabled(False)
        gui.lig7occ.setDisabled(False)
        gui.lig7conn.setDisabled(False)
        gui.lig7H.setDisabled(False)
        gui.lig7ML.setDisabled(False)
        gui.lig7an.setDisabled(False)
        gui.lig7nam.setDisabled(False)
        gui.lig7add.setDisabled(False)

#####################################################
########## write options to input file  #############
#####################################################
def writeinputc(args,fname):
    f = open(fname,'w')
    f.write("# Input file generated from GUI options\n")
    for key, val in args.iteritems():
        if len(val) > 0:
            vals = val.splitlines()
            for v in vals:
                f.write(key+' '+v+'\n')
    f.close()

###########################################################
########## write options to postp input file  #############
###########################################################
def writeinputp(args,fname):
    f = open(fname,'w')
    f.write("# Input file generated from GUI options\n")
    for key, val in args.iteritems():
            f.write(key+' '+val+'\n')
    f.close()
    
#####################################################
########## write options to input file  #############
#####################################################
def writeinputf(args):
    f = open(args['-rundir']+'/geninput.inp','w')
    f.write("# Input file generated from GUI options\n")
    ff = open(args['-rundir']+'/molSimp.log','a')
    dd = "# Input file generated from GUI options at "+ time.strftime('%m/%d/%Y %H:%M')+'\n'
    ff.write(dd)
    for key, val in args.iteritems():
        if len(val) > 0:
            vals = val.splitlines()
            for v in vals:
                f.write(key+' '+v+'\n')
                ff.write(key+' '+v+'\n')
    f.close()
    ff.close()
    
#########################################################
########## grabs GUI options to input file  #############
#########################################################
def grabguivars(gui):
    # list with arguments
    args = dict()
    ### general structure generation options ###
    args['-core'] = gui.etcore.text()
    args['-ccatoms'] = gui.etccat.text()
    ligs,ligoccs,lcats,kHs,MLb,lang,lname=getligands(gui)
    noligs = len(ligs.split(','))-1
    args['-lig'] = ligs
    if len(ligoccs) > noligs:
        args['-ligocc'] = ligoccs
    if len(MLb) > noligs:
        args['-MLbonds'] = MLb
    if len(lang) > noligs:
        args['-pangles'] = lang
    if len(kHs) > noligs:
        args['-keepHs'] = kHs
    if len(lcats) > noligs:
        args['-smicat'] = lcats
    if len(lname) > noligs:
        args['-sminame'] = lname
    if gui.replig.getState():
        args['-replig'] = '1'
    if gui.ligfloc.getState():
        args['-ligloc'] = '1'
    if gui.ligfalign.getState():
        args['-ligalign'] = '1'
    if gui.chkgenall.getState():
        args['-genall'] = '1'
    args['-coord'] = gui.dcoord.currentText()
    args['-geometry'] = gui.dcoordg.currentText()
    if gui.randomchk.getState():
        args['-lignum'] = gui.etlignum.currentText()
        args['-liggrp'] = gui.etliggrp.currentText()
        args['-ligctg'] = gui.etligctg.currentText()
        args['-rgen'] = gui.etrgen.text()
        if gui.randkHs.getState():
            args['-rkHs'] = 'yes'
    args['-distort'] = str(gui.sdist.value())
    args['-oxstate'] = gui.doxs.currentText()
    args['-spin'] = gui.dspin.currentText()
    args['-rundir'] = gui.etrdir.text()
    rdir = args['-rundir']
    if rdir[1] == '/':
         args['-rundir'] = rdir[:-1]
    args['-suff'] = gui.etsuff.text()
    ### binding molecule options ###
    if gui.chkM.getState():
        args['-bind'] = gui.etbind.text()
        args['-bcharge'] = gui.etchbind.text()
        args['-nbind'] = gui.etnbind.text()
        if gui.chsep.getState():
            args['-bsep'] = 'yes'
        args['-nambsmi'] = gui.etbsmi.text()
        args['-maxd'] = gui.etplacemax.text()
        args['-mind'] = gui.etplacemin.text()
        args['-place'] = gui.dmolp.currentText()
        args['-bphi'] = gui.etplacephi.text()
        args['-btheta'] = gui.etplacetheta.text()
        args['-bref'] = gui.etmaskbind.text()
    ### force field optimization ###
    if gui.chkFF.getState():
        args['-ff'] = gui.dff.currentText()
        args['-ffoption'] = gui.dffba.currentText()
    ### Quantum Chemistry options ###
    if gui.chkI.getState():
        args['-qccode'] = gui.qcode.currentText()
        if args['-qccode'].lower() in 'terachem':
            args['-charge'] = gui.etqctch.text()
            if gui.chch.getState():
                args['-calccharge'] = 'yes'
            args['-spin'] = gui.etqctspin.text()
            print gui.qctcalc.currentText()
            args['-runtyp'] = gui.qctcalc.currentText()
            args['-method'] = gui.etqctmethod.text()
            args['-basis'] = gui.etqctbasis.text()
            args['-dispersion'] = gui.qctsel.currentText()
            args['-qoption'] = gui.qceditor.toPlainText()
        elif args['-qccode'].lower() in 'gamess':
            args['-charge'] = gui.etqcgch.text()
            if gui.chch.getState():
                args['-calccharge'] = 'yes'
            args['-spin'] = gui.etqcgspin.text()
            args['-runtyp'] = gui.qcgcalc.currentText()
            args['-method'] = gui.etqcgmethod.text()
            args['-gbasis'] = gui.etqcgbasis.text()
            args['-ngauss'] = gui.etqcngauss.text()
            args['-npfunc'] = gui.etqcnpfunc.text()
            args['-ndfunc'] = gui.etqcndfunc.text()
            args['-sysoption'] = gui.qcgedsys.toPlainText()
            args['-ctrloption'] = gui.qcgedctrl.toPlainText()
            args['-scfoption'] = gui.qcgedscf.toPlainText()
            args['-statoption'] = gui.qcgedstat.toPlainText()
        elif args['-qccode'].lower() in 'qchem':
            args['-charge'] = gui.etqcQch.text()
            if gui.chch.getState():
                args['-calccharge'] = 'yes'
            args['-spin'] = gui.etqcQspin.text()
            args['-runtyp'] = gui.qcQcalc.currentText()
            args['-basis'] = gui.etqcQbasis.text()
            args['-remoption'] = gui.qcQeditor.toPlainText()
            args['-exchange'] = gui.etqcQex.text()
            args['-correlation'] = gui.etqcQcor.text()
            if gui.chQun:
                args['-unrestricted'] = '1'
    ### jobscript options ###
    if gui.chkJ.getState():
        args['-jsched'] = gui.scheduler.currentText()
        args['-jname'] = gui.etjname.text()
        args['-memory'] = gui.etjmem.text()
        args['-wtime'] = gui.etjwallt.text()
        args['-queue'] = gui.etjqueue.text()
        args['-gpus'] = gui.etjgpus.text()
        args['-cpus'] = gui.etjcpus.text()
        args['-modules'] = gui.etjmod.text()
        args['-joption'] = gui.etjopt.toPlainText()
        args['-jcommand'] = gui.jcomm.toPlainText()
    writeinputf(args)
    return args

#########################################################
##########  grabs GUI options for terachem  #############
#########################################################
def grabguivarstc(gui):
    globs = globalvars()
    # list with arguments
    args = dict()
    args['-charge'] = gui.etqctch.text()
    args['-spin'] = gui.etqctspin.text()
    args['-runtyp'] = gui.qctcalc.currentText()
    args['-method'] = gui.etqctmethod.text()
    args['-basis'] = gui.etqctbasis.text()
    args['-dispersion'] = gui.qctsel.currentText()
    args['-qoption'] = gui.qceditor.toPlainText()
    ### write input file ###
    writeinputc(args,globs.installdir+'/Data/.tcdefinput.inp')
    return args

######################################################
##########  grabs GUI options for GAMESS #############
######################################################
def grabguivarsgam(gui):
    globs = globalvars()
    # list with arguments
    args = dict()
    args['-charge'] = gui.etqcgch.text()
    args['-spin'] = gui.etqcgspin.text()
    args['-runtyp'] = gui.qcgcalc.currentText()
    args['-method'] = gui.etqcgmethod.text()
    args['-gbasis'] = gui.etqcgbasis.text()
    args['-ngauss'] = gui.etqcngauss.text()
    args['-npfunc'] = gui.etqcnpfunc.text()
    args['-ndfunc'] = gui.etqcndfunc.text()
    args['-sysoption'] = gui.qcgedsys.toPlainText()
    args['-ctrloption'] = gui.qcgedctrl.toPlainText()
    args['-scfoption'] = gui.qcgedscf.toPlainText()
    args['-statoption'] = gui.qcgedstat.toPlainText()
    ### write input file ###
    writeinputc(args,globs.installdir+'/Data/.gamdefinput.inp')
    return args

#####################################################
##########  grabs GUI options for QChem #############
#####################################################
def grabguivarsqch(gui):
    globs = globalvars()
    # list with arguments
    args = dict()
    args['-charge'] = gui.etqcQch.text()
    args['-spin'] = gui.etqcQspin.text()
    args['-runtyp'] = gui.qcQcalc.currentText()
    args['-basis'] = gui.etqcQbasis.text()
    args['-remoption'] = gui.qcQeditor.toPlainText()
    args['-exchange'] = gui.etqcQex.text()
    args['-correlation'] = gui.etqcQcor.text()
    if gui.chQun:
        args['-unrestricted'] = '1'
    ### write input file ###
    writeinputc(args,globs.installdir+'/Data/.qchdefinput.inp')
    return args

#########################################################
##########  grabs GUI options for jobscript #############
#########################################################
def grabguivarsjob(gui):
    globs = globalvars()
    # list with arguments
    args = dict()
    args['-jname'] = gui.etjname.text()
    args['-memory'] = gui.etjmem.text()
    args['-wtime'] = gui.etjwallt.text()
    args['-queue'] = gui.etjqueue.text()
    args['-gpus'] = gui.etjgpus.text()
    args['-cpus'] = gui.etjcpus.text()
    args['-modules'] = gui.etjmod.text()
    args['-joption'] = gui.etjopt.toPlainText()
    args['-jcommand'] = gui.jcomm.toPlainText()
    ### write input file ###
    writeinputc(args,globs.installdir+'/Data/.jobdefinput.inp')
    return args

############################################################
########## grabs GUI db options to input file  #############
############################################################
def grabdbguivars(gui):
    args = dict()
    ### database search options ###
    args['-dbsim'] = gui.etcDBsmi.text()
    args['-dbcatoms'] = gui.etcDBcatoms.text()
    args['-dbresults'] = gui.etcDBnres.text()
    args['-dboutputf'] = gui.etcDBoutf.text()+gui.cDBdent.currentText()
    args['-dbbase'] = gui.cDBsel.currentText()
    args['-dbsmarts'] = gui.etcDBsmarts.text()
    args['-dbfinger'] = gui.cDBsf.currentText()
    args['-dbatoms'] = gui.etcDBsatoms0.text()+'<'+gui.etcDBsatoms1.text()
    args['-dbbonds'] = gui.etcDBsbonds0.text()+'<'+gui.etcDBsbonds1.text()
    args['-dbarbonds'] = gui.etcDBsabonds0.text()+'<'+gui.etcDBsabonds1.text()
    args['-dbsbonds'] = gui.etcDBsbondss0.text()+'<'+gui.etcDBsbondss1.text()
    args['-dbmw'] = gui.etcDBmw0.text()+'<'+gui.etcDBmw1.text()
    if len(args['-dboutputf'].replace(' ',''))==4:
        args['-dboutputf'] = ''
    ### check input
    if args['-dbatoms']=='<':
        args['-dbatoms'] = ''
    if args['-dbbonds']=='<':
        args['-dbbonds'] = ''
    if args['-dbarbonds']=='<':
        args['-dbarbonds'] = ''
    if args['-dbsbonds']=='<':
        args['-dbsbonds'] = ''
    if args['-dbmw']=='<':
        args['-dbmw'] = ''
    rdir = gui.etrdir.text()
    if rdir[-1]=='/':
        rdir = rdir[:-1]
    args['-rundir'] = rdir
    ### write input file ###
    writeinputc(args,rdir+'/dbinput.inp')
    return args
    
    
###################################################################
########## grabs GUI options to input file for postp  #############
###################################################################
def grabguivarsP(gui):
    args = dict()
    ### post-processing options ###
    args['-postp'] = ''
    args['-postdir'] = gui.etpdir.text()
    args['-postqc'] = gui.pqcode.currentText()
    rdir = args['-postdir']
    if rdir[-1]=='/':
        rdir = rdir[:-1]
    if gui.psum.getState()==1:
        args['-pres'] = ''
    if gui.pwfnav.getState()==1:
        args['-pwfninfo'] = ''
    if gui.pch.getState()==1: 
        args['-pcharge'] = ''
    if gui.pcub.getState()==1:
        args['-pgencubes'] = ''
    if gui.porbs.getState()==1: 
        args['-porbinfo'] = ''
    if gui.pdeloc.getState()==1:
        args['-pdeloc'] = ''
    if gui.pnbo.getState()==1:
        args['-pnbo'] = ''
    ### write input file ###
    writeinputp(args,rdir+'/postproc.inp')
#################################################
########### loads input file to tc  #############
#################################################
def loadfrominputtc(gui,fname):
    f = open(fname,'r')
    s = f.read()
    s = filter(None,s.splitlines())
    f.close()
    gui.qceditor.setText('')
    ### general structure generation options ###
    for ss in s:
        st = ss.split(None,1)
        if '-charge'==st[0]:
            gui.etqctch.setText(st[-1])
        if '-spin'==st[0]:
            gui.etqctspin.setText(st[-1])
        if '-runtyp'==st[0]:
            gui.qctcalc.setCurrentText(st[-1])
        if '-method'==st[0]:
            gui.etqctmethod.setText(st[-1])
        if '-basis'==st[0]:
            gui.etqctbasis.setText(st[-1])
        if '-dispersion'==st[0]:
            gui.qctsel.setCurrentText(st[-1])
        if '-qoption'==st[0]:
            gui.qceditor.setText(gui.qceditor.toPlainText()+st[-1]+'\n')
        if '-charge'==st[0]:
            gui.etqctch.setText(st[-1])
        if '-runtyp'==st[0]:
            gui.qctcalc.setCurrentText(st[-1])
#####################################################
########### loads input file to GAMESS  #############
#####################################################
def loadfrominputgam(gui,fname):
    f = open(fname,'r')
    s = f.read()
    s = filter(None,s.splitlines())
    f.close()
    gui.qcgedsys.setText('')
    gui.qcgedctrl.setText('')
    gui.qcgedscf.setText('')
    gui.qcgedstat.setText('')
    ### general structure generation options ###
    for ss in s:
        st = ss.split(None,1)
      ### Quantum Chemistry options ###
        if '-charge'==st[0]:
            gui.etqcgch.setText(st[-1])
        if '-spin'==st[0]:
            gui.etqcgspin.setText(st[-1])
        if '-runtyp'==st[0]:
            gui.qcgcalc.setCurrentText(st[-1])
        if '-method'==st[0]:
            gui.etqcgmethod.setText(st[-1])
        if '-gbasis'==st[0]:
            gui.etqcgbasis.setText(st[-1])
        if '-ngauss'==st[0]:
            gui.etqcngauss.setText(st[-1])
        if '-npfunc'==st[0]:
            gui.etqcnpfunc.setText(st[-1])
        if '-ndfunc'==st[0]:
            gui.etqcndfunc.setText(st[-1])
        if '-sysoption'==st[0]:
            gui.qcgedsys.setText(gui.qcgedsys.toPlainText()+st[-1]+'\n')
        if '-ctrloption'==st[0]:
            gui.qcgedctrl.setText(gui.qcgedctrl.toPlainText()+st[-1]+'\n')
        if '-scfoption'==st[0]:
            gui.qcgedscf.setText(gui.qcgedscf.toPlainText()+st[-1]+'\n')
        if '-statoption'==st[0]:
            gui.qcgedstat.setText(gui.qcgedstat.toPlainText()+st[-1]+'\n')
####################################################
########### loads input file to QChem  #############
####################################################
def loadfrominputqch(gui,fname):
    f = open(fname,'r')
    s = f.read()
    s = filter(None,s.splitlines())
    f.close()
    gui.qcQeditor.setText('')
    ### general structure generation options ###
    for ss in s:
        st = ss.split(None,1)
      ### Quantum Chemistry options ###
        if '-basis'==st[0]:
            gui.etqcQbasis.setText(st[-1])
        if '-charge'==st[0]:
            gui.etqcQch.setText(st[-1])
        if '-spin'==st[0]:
            gui.etqcQspin.setText(st[-1])
        if '-runtyp'==st[0]:
            gui.qcQcalc.setCurrentText(st[-1])
        if '-remoption'==st[0]:
            gui.qcQeditor.setText(gui.qcQeditor.toPlainText()+st[-1]+'\n')
        if '-exchange'==st[0]:
            gui.etqcQex.setText(st[-1])
        if '-correlation'==st[0]:
            gui.etqcQcor.setText(st[-1])
        if '-unrestricted'==st[0]:
            gui.chQun.setChecked(True)
#########################################################
########### loads input file to jobscripts  #############
#########################################################
def loadfrominputjob(gui,fname):
    f = open(fname,'r')
    s = f.read()
    s = filter(None,s.splitlines())
    f.close()
    gui.etjopt.setText('')
    gui.jcomm.setText('')
    ### general structure generation options ###
    for ss in s:
        st = ss.split(None,1)
        ### jobscript options ###
        if '-jname'==st[0]:
            gui.etjname.setText(st[-1])
        if '-memory'==st[0]:
            gui.etjmem.setText(st[-1])
        if '-wtime'==st[0]:
            gui.etjwallt.setText(st[-1])
        if '-queue'==st[0]:
            gui.etjqueue.setText(st[-1])
        if '-gpus'==st[0]:
            gui.etjgpus.setText(st[-1])
        if '-cpus'==st[0]:
            gui.etjcpus.setText(st[-1])
        if '-modules'==st[0]:
            gui.etjmod.setText(st[-1])
        if '-joption'==st[0]:
            gui.etjopt.setText(gui.etjopt.toPlainText()+st[-1]+'\n')
        if '-jcommand'==st[0]:
            gui.jcomm.setText(gui.jcomm.toPlainText()+st[-1]+'\n')
#################################################
########## loads input file to GUI  #############
#################################################
def loadfrominputfile(gui,fname):
    f = open(fname,'r')
    s = f.read()
    s = filter(None,s.splitlines())
    f.close()
    db = False
    ligs,ligoccs,lcats,kHs,MLb,lang,lname='','','','','','',''
    # reset edit boxes
    gui.qceditor.setText('')
    gui.qcgedsys.setText('')
    gui.qcgedctrl.setText('')
    gui.qcgedscf.setText('')
    gui.qcgedstat.setText('')
    gui.qcQeditor.setText('')
    gui.etjopt.setText('')
    gui.jcomm.setText('')
    ### general structure generation options ###
    for ss in s:
        st = ss.split(None,1)
        if '-lig' not in st[0]:
            st[-1] = st[-1].split('#')[0]
        if '-core'==st[0]:
            gui.etcore.setText(st[-1])
        if '-lig'==st[0]:
            ligs = st[-1]
        if '-ligocc'==st[0]:
            ligoccs= st[-1]
        if '-replig'==st[0]:
            if checkTrue(st[-1]):
                gui.replig.setChecked(True)
        if '-ligloc'==st[0]:
            if checkTrue(st[-1]):
                gui.ligfloc.setChecked(True)
        if '-ligalign'==st[0]:
            if checkTrue(st[-1]):
                gui.ligfalign.setChecked(True)
        if '-genall'==st[0]:
            if checkTrue(st[-1]):
                gui.chkgenall.setChecked(True)
        if '-MLbonds'==st[0]:
            Mlb = st[-1]
        if '-distort'==st[0]:
            gui.sdist.setValue(int(st[-1]))
        if '-pangles'==st[0]:
            lang = st[-1]
        if '-ccatoms'==st[0]:
            gui.etccat.setText(st[-1])
        if '-coord'==st[0]:
            gui.dcoord.setCurrentText(st[-1])
        if '-geometry'==st[0]:
            gui.dcoordg.setCurrentText(st[-1])
        if '-lignum'==st[0]:
            gui.etlignum.setCurrentText(st[-1])
        if '-liggrp'==st[0]:
            gui.etliggrp.setCurrentText(st[-1])
        if '-ligctg'==st[0]:
            gui.etligctg.setCurrentText(st[-1])
        if '-rkHs'==st[0]:
            if checkTrue(st[-1]):
                gui.randkHs.setChecked(True)
        if '-rgen'==st[0]:
            gui.etrgen.setText(st[-1])
            gui.randomchk.setChecked(True)
        if '-keepHs'==st[0]:
            kHs = st[-1]
        if '-smicat'==st[0]:
            lcats = st[-1]
        if '-sminame'==st[0]:
            lname = st[-1]
        if '-oxstate'==st[0]:
            gui.doxs.setCurrentText(st[-1])
        if '-rundir'==st[0]:
            gui.etrdir.setText(st[-1])
        if '-suff'==st[0]:
            gui.etsuff.setText(st[-1])
        ### binding molecule options ###
        if '-bind'==st[0]:
            gui.chkM.setChecked(True)
            gui.enableemol()
            gui.etbind.setText(st[-1])
        if '-bsep'==st[0]:
            gui.chsep.setChecked(True)
        if '-bcharge'==st[0]:
            gui.etchbind.setText(st[-1])
        if '-nbind'==st[0]:
            gui.etnbind.setText(st[-1])
        if '-nambsmi'==st[0]:
            gui.etbsmi.setText(st[-1])
        if '-maxd'==st[0]:
            gui.etplacemax.setText(st[-1])
        if '-mind'==st[0]:
            gui.etplacemin.setText(st[-1])
        if '-place'==st[0]:
            gui.dmolp.setCurrentText(st[-1])
        if '-bphi'==st[0]:
            gui.etplacephi.setText(st[-1])
        if '-btheta'==st[0]:
            gui.etplacetheta.setText(st[-1])
        if '-bref'==st[0]:
            gui.etmaskbind.setText(st[-1])
        ### force field optimization ###
        ff = False
        if '-ff'==st[0]:
            gui.chkFF.setChecked(True)
            gui.enableffinput()
            gui.dff.setCurrentText(st[-1])
            ff = True
        if '-ffoption'==st[0]:
            sopt = st[-1].split('&')
            b,a,e = False,False,False
            for sopts in sopt:
                ssopt = sopts.replace(' ','')
                if ssopt[0].lower()=='b':
                    b = True
                elif ssopt[0].lower()=='a':
                    a = True
            if b and not a:
                idx = 0
            elif not b and a:
                idx = 1
            elif b and a:
                idx = 2
            else:
                idx = 2
            gui.dffba.setCurrentIndex(idx)
        ### Quantum Chemistry options ###
        if '-qccode'==st[0]:
            gui.chkI.setChecked(True)
            gui.enableqeinput()
            gui.qcode.setCurrentText(st[-1])
        if '-calccharge'==st[0]:
            gui.chch.setChecked(True)
        if '-basis'==st[0]:
            gui.etqctbasis.setText(st[-1])
            gui.etqcQbasis.setText(st[-1])
        if '-dispersion'==st[0]:
            gui.qctsel.setCurrentText(st[-1])
        if '-qoption'==st[0]:
            gui.qceditor.setText(gui.qceditor.toPlainText()+st[-1]+'\n')
        if '-charge'==st[0]:
            gui.etqcgch.setText(st[-1])
            gui.etqctch.setText(st[-1])
            gui.etqcQch.setText(st[-1])
        if '-spin'==st[0]:
            gui.etqcgspin.setText(st[-1])
            gui.etqctspin.setText(st[-1])
            gui.etqcQspin.setText(st[-1])
        if '-runtyp'==st[0]:
            gui.qcgcalc.setCurrentText(st[-1])
            gui.qctcalc.setCurrentText(st[-1])
            gui.qcQcalc.setCurrentText(st[-1])
        if '-method'==st[0]:
            gui.etqctmethod.setText(st[-1])
            gui.etqcgmethod.setText(st[-1])
        if '-gbasis'==st[0]:
            gui.etqcgbasis.setText(st[-1])
        if '-ngauss'==st[0]:
            gui.etqcngauss.setText(st[-1])
        if '-npfunc'==st[0]:
            gui.etqcnpfunc.setText(st[-1])
        if '-ndfunc'==st[0]:
            gui.etqcndfunc.setText(st[-1])
        if '-sysoption'==st[0]:
            gui.qcgedsys.setText(gui.qcgedsys.toPlainText()+st[-1]+'\n')
        if '-ctrloption'==st[0]:
            gui.qcgedctrl.setText(gui.qcgedctrl.toPlainText()+st[-1]+'\n')
        if '-scfoption'==st[0]:
            gui.qcgedscf.setText(gui.qcgedscf.toPlainText()+st[-1]+'\n')
        if '-statoption'==st[0]:
            gui.qcgedstat.setText(gui.qcgedstat.toPlainText()+st[-1]+'\n')
        if '-remoption'==st[0]:
            gui.qcQeditor.setText(gui.qcQeditor.toPlainText()+st[-1]+'\n')
        if '-exchange'==st[0]:
            gui.etqcQex.setText(st[-1])
        if '-correlation'==st[0]:
            gui.etqcQcor.setText(st[-1])
        if '-unrestricted'==st[0]:
            gui.chQun.setChecked(True)
        ### jobscript options ###
        if '-jsched'==st[0]:
            gui.chkJ.setChecked(True)
            gui.enableqeinput()
            gui.scheduler.setCurrentText(st[-1])
        if '-jname'==st[0]:
            gui.etjname.setText(st[-1])
        if '-memory'==st[0]:
            gui.etjmem.setText(st[-1])
        if '-wtime'==st[0]:
            gui.etjwallt.setText(st[-1])
        if '-queue'==st[0]:
            gui.etjqueue.setText(st[-1])
        if '-gpus'==st[0]:
            gui.etjgpus.setText(st[-1])
        if '-cpus'==st[0]:
            gui.etjcpus.setText(st[-1])
        if '-modules'==st[0]:
            gui.etjmod.setText(st[-1])
        if '-joption'==st[0]:
            gui.etjopt.setText(gui.etjopt.toPlainText()+st[-1]+'\n')
        if '-jcommand'==st[0]:
            gui.jcomm.setText(gui.jcomm.toPlainText()+st[-1]+'\n')
        ### database search options ###
        if '-dbsim'==st[0]:
            gui.etcDBsmi.setText(st[-1])
            db = True
        if '-dbcatoms'==st[0]:
            gui.etcDBcatoms.setText(st[-1])
        if '-dbresults'==st[0]:
            gui.etcDBnres.setText(st[-1])
            db = True
        if '-dboutputf'==st[0]:
            ssp = st[-1].split('.')
            gui.etcDBoutf.setText(ssp[0])
            gui.cDBdent.setCurrentText('.'+ssp[-1])
        if '-dbbase'==st[0]:
            gui.cDBsel.setCurrentText(st[-1])
            db = True
        if '-dbsmarts'==st[0]:
            gui.etcDBsmarts.setText(st[-1])
        if '-dbfinger'==st[0]:
            gui.cDBsf.setCurrentText(st[-1])
        if '-dbatoms'==st[0]:
            ssp = st[-1].split('<')
            gui.etcDBsatoms0.setText(ssp[0])
            gui.etcDBsatoms1.setText(ssp[-1])
        if '-dbbonds'==st[0]:
            ssp = st[-1].split('<')
            gui.etcDBsbonds0.setText(ssp[0])
            gui.etcDBsbonds1.setText(ssp[-1])
        if '-dbarbonds'==st[0]:
            ssp = st[-1].split('<')
            gui.etcDBsabonds0.setText(ssp[0])
            gui.etcDBsabonds1.setText(ssp[-1])
        if '-dbsbonds'==st[0]:
            ssp = st[-1].split('<')
            gui.etcDBsbondss0.setText(ssp[0])
            gui.etcDBsbondss1.setText(ssp[-1])
        if '-dbmw'==st[0]:
            ssp = st[-1].split('<')
            gui.etcDBmw0.setText(ssp[0])
            gui.etcDBmw1.setText(ssp[-1])
    setligands(gui,ligs,ligoccs,lcats,kHs,MLb,lang,lname)
    if db:
        gui.searchDBW()
