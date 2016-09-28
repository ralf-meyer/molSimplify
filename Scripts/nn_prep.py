# Written by JP Janet for HJK Group
# Dpt of Chemical Engineering, MIT

##########################################################
####### This script is a collection of helper ############
########  routines for ANN integration in    #############
########           molsimplify               #############
##########################################################

# import custom modules
from geometry import *
from io import *
from Classes.globalvars import *
# import standard modules

def check_ligands(ligs,batlist,dents):
    n_ligs = len(ligs)
    unique_ligans = []
    axial_ligs = []
    equitorial_ligs = []
    ax_dent = []
    eq_dent  =[]
    eq_ligs = []
    valid = True
    if  (set(dents) == [2]):
        print('triple bidentate case\n')
        unique_ligs = []
        if not(n_ligs) == 3:
                ## something unpe
                valid = False 
        for i in range(0,n_ligs):
            this_bat = batlist[i]
            this_lig = ligs[i]
            this_dent = dents[i]
            ## mulitple points
            print('this bat',this_bat)
            if not (this_lig in unique_ligs):
                    print('adding unique ligs',this_lig)
            elif this_lig in unique_ligs:
                   equitorial_ligs.append(this_lig)
                   eq_dent.append(this_dent)
        for uligs in unique_ligs:
            if not (uligs in equitorial_ligs): #only occured once
                axial_ligs.append(this_lig)
                ax_dent.append(this_dent)
        if not (len(unique_ligs) == 3):
            valid = False
    else:
        print('regular case\n')
        for i in range(0,n_ligs):
            this_bat = batlist[i]
            this_lig = ligs[i]
            this_dent = dents[i]
            ## mulitple points
            print('this bat',this_bat)
            if len(this_bat) == 1:
                if (5 in this_bat) or (6 in this_bat):
                    if not (this_lig in axial_ligs):
                        print('adding axial')
                        axial_ligs.append(this_lig)
                        ax_dent.append(this_dent)
                else:
                    if not (this_lig in equitorial_ligs):
                        equitorial_ligs.append(this_lig)
                        eq_dent.append(this_dent)
            else:
                if not (this_lig in equitorial_ligs):
                        equitorial_ligs.append(this_lig)
                        eq_dent.append(this_dent)
    if not (len(axial_ligs) == 1):
        print('axial ligs mismatch: ',axial_ligs,ax_dent)
        valid = False
    if not (len(equitorial_ligs) == 1):
        print('equitorial ligs mismatch: ',equitorial_ligs,eq_dent)
        valid = False
    return valid,axial_ligs,equitorial_ligs,ax_dent,eq_dent

def check_metal(metal,oxidation_state):
    supported_metal_dict = {"Fe":[2,3],"Mn":[2,3],"Cr":[2,3],
                            "Co":[2,3],"Ni":[2]}
    romans={'I':'1','II':'2','III':'3','IV':'4','V':'5','VI':'6'}
    if oxidation_state  in romans.keys():
        oxidation_state= romans[oxidation_state]
    outcome = False
    if metal in supported_metal_dict.keys():
        if oxidation_state in supported_metal_dict[metal]:
            outcome = True
    return outcome

def examine_inputs(args,ligs,occs,dents,batslist,installdir,licores):
    emsg = list()
    input_ok = True 
    metal = args.core
    if not args.geometry == "oct":
        emsg += "\n [ANN] Geometry is not supported at this time"
        input_ok = False
    if not args.oxstate:
        emsg += "\n [ANN] oxidation state must be given"
        inputs_ok = False
    if input_ok:
        oxidation_state = args.oxstate
        check_metal(metal,oxidation_state)
        ## generate key in descriptor space
        this_metal = metal.lower()
        this_ox = oxidation_state
    valid,axial_ligs,equitorial_ligs,ax_dent,eq_dent = check_ligands(ligs,batslist,dents)

    print("\n\n")
    print('Here comes occs')
    print(occs)
    print('Ligands')
    print(ligs)
    print('Here comes dents')
    print(dents)
    print('Here comes bats')
    print(batslist)
    print('lig validity',valid)
    print('ax ligs',axial_ligs)
    print('eq ligs',equitorial_ligs)

    print("\n\n")
    if valid:
            ax_lig3D,emsg = lig_load(installdir,axial_ligs[0],licores) # load ligand
            if emsg:
                    print(esmg)
                    return False
            ax_lig3D.convert2mol3D() ## mol3D representation of ligand

            eq_lig3D,emsg = lig_load(installdir,equitorial_ligs[0],licores) # load ligand
            if emsg:
                    print(esmg)
                    return False
            eq_lig3D.convert2mol3D() ## mol3D representation of ligand

            print(eq_lig3D)
            print(ax_lig3D)

#this_key


