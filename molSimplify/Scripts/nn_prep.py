# Written by JP Janet for HJK Group
# Dpt of Chemical Engineering, MIT

##########################################################
####### This script is a collection of helper ############
########  routines for ANN integration in    #############
########           molsimplify               #############
##########################################################

# import custom modules
from molSimplify.Scripts.geometry import *
from molSimplify.Scripts.io import *
from molSimplify.Classes.globalvars import *
from molSimplify.python_nn.graph_analyze import *
from molSimplify.python_nn.pytry import *
# import standard modules
import openbabel


def get_bond_order(OBMol):
    ## informs the ANN
    ## of the highest bond
    ## order in ligand
    # INPUT:
    #   - OBMol:  OBMol class ligand
    # OUTPUT:
    #   - max_bo: int, max bond order
    OBMol.PerceiveBondOrders()
    max_bo = 0
    for bonds in openbabel.OBMolBondIter(OBMol):
       this_BO = bonds.GetBondOrder()
       if this_BO > max_bo:
           max_bo = this_BO
    return max_bo


def check_ligands(ligs,batlist,dents,tcats):
    ## tests if ligand combination
    ## is compatiable with the ANN
    # INPUT:
    #   - ligs:  list of mol3D class, ligands
    #   - batlist: list of int, occupations 
    #   - dents: list of int, denticity
    #   - tcats: list of int/bool
    # OUTPUT:
    #   - valid: bool
    ## tcats controls
    ## manual overide
    ## of connection atoms

    n_ligs = len(ligs)
    unique_ligans = []
    axial_ligs = []
    equitorial_ligs = []
    ax_dent =0 
    eq_dent  =0
    eq_ligs = []
    eq_tcat = False
    ax_tcat = False
    valid = True
    if  (set(dents) == set([2])):
        print('triple bidentate case\n')
        unique_ligs = []
        ucats = []
        if not(n_ligs) == 3:
                ## something unexpected happened!
                valid = False 
        for i in range(0,n_ligs):
            this_bat = batlist[i]
            this_lig = ligs[i]
            this_dent = dents[i]
            ## mulitple points
            if not (this_lig in unique_ligs):
                    print('adding unique ligs',this_lig)
                    unique_ligs.append(this_lig)
                    ucats.append(this_dent)
            elif (this_lig in unique_ligs) and (not this_lig in equitorial_ligs) :
                   equitorial_ligs.append(this_lig)
                   eq_dent = this_dent
                   eq_tcat = tcats[i]
        if len(unique_ligs) == 1:
            axial_ligs.append(equitorial_ligs[0])
            ax_dent  = 2
            ax_tcat = eq_tcat
        elif len(unique_ligs) == 2:
            for i,uligs in enumerate(unique_ligs):
                if not (uligs in equitorial_ligs): #only occured once
                    axial_ligs.append(this_lig)
                    ax_dent = 2
                    ax_tcat = ucats[i]
        else:
            valid = False
    else:
        for i in range(0,n_ligs):
            this_bat = batlist[i]
            this_lig = ligs[i]
            this_dent = dents[i]
            ## mulitple points
            if len(this_bat) == 1:
                if (5 in this_bat) or (6 in this_bat):
                    if not (this_lig in axial_ligs):
                        axial_ligs.append(this_lig)
                        ax_dent = this_dent
                        ax_tcat = tcats[i]
                else:
                    if not (this_lig in equitorial_ligs):
                        equitorial_ligs.append(this_lig)
                        eq_dent = this_dent
                        eq_tcat = tcats[i]
            else:
                if not (this_lig in equitorial_ligs):
                        equitorial_ligs.append(this_lig)
                        eq_dent = this_dent
                        eq_tcat = tcats[i]
    if not (len(axial_ligs) == 1):
        print('axial ligs do not match, have  ' + str(' '.join(axial_ligs)))
        valid = False
    if not (len(equitorial_ligs) == 1):
        print('equitorial ligs do not match, have ' + str(' '.join(equitorial_ligs)))
        valid = False
    return valid,axial_ligs,equitorial_ligs,ax_dent,eq_dent,ax_tcat,eq_tcat

def check_metal(metal,oxidation_state):
    supported_metal_dict = {"Fe":[2,3],"Mn":[2,3],"Cr":[2,3],
                            "Co":[2,3],"Ni":[2]}
    romans={'I':'1','II':'2','III':'3','IV':'4','V':'5','VI':'6'}
    if oxidation_state  in romans.keys():
        oxidation_state= romans[oxidation_state]
    outcome = False
#    print('incheck',oxidation_state)
    if metal in supported_metal_dict.keys():
#        print('metal in',supported_metal_dict[metal])
        if int(oxidation_state) in supported_metal_dict[metal]:
            outcome = True
    return outcome,oxidation_state


def get_truncated_kier(ligand,connection_atoms):
        ### three hop truncation
        trunc_mol = obtain_truncation(ligand,connection_atoms,3)
        this_kier = kier(trunc_mol) 
        return this_kier
def get_con_at_type(mol,connection_atoms):
    this_type = ""
    been_set = False
    valid = True
    for atoms in connection_atoms:
        this_symbol = mol.getAtom(atoms).symbol()
        if not (this_symbol == this_type):
            if not been_set:
                this_type = this_symbol
            else:
                print('different connects')
                valid = False
    if not this_type in ['C','O','Cl','N','F']:
        valid = False
        print('untrained atom type: ',this_type)
    return valid,this_type


def ANN_preproc(args,ligs,occs,dents,batslist,tcats,licores):
    nn_excitation = []
    r = 0
    emsg = list()
    valid = True 
    metal = args.core
    if not args.geometry == "oct":
#        print('ANN: geom  is: ' + str(args.geometry))
        emsg.append("\n [ANN] Geometry is not supported at this time, MUST give -geometry = oct")
        valid = False
    if not args.oxstate:
        emsg.append("\n [ANN] oxidation state must be given")
        valid = False
    if valid:
        oxidation_state = args.oxstate
        valid, oxidation_state = check_metal(metal,oxidation_state)
        ## generate key in descriptor space
        this_metal = metal.lower()
        ox = int(oxidation_state)
        spin = args.spin
        #print('metal validity',valid)
        if not valid:
            emsg.append("\n [ANN] Oxidation state not available for this metal")

    if valid:
        high_spin = spin_classify(this_metal,spin,ox)
        if not valid:
            emsg.append("\n [ANN] this spin state not available for this metal")

    print('ANN status: '+ str(' '.join(emsg)))
    if valid:
        valid,axial_ligs,equitorial_ligs,ax_dent,eq_dent,ax_tcat,eq_tcat = check_ligands(ligs,batslist,dents,tcats)

       # print("\n")
       # print('Here comes occs')
       # print(occs)
       # print('Ligands')
       # print(ligs)
       #  print('Here comes dents')
       # print(dents)
       # print('Here comes bats')
       # print(batslist)
       # print('lig validity',valid)
       # print('ax ligs',axial_ligs)
       # print('eq ligs',equitorial_ligs)
       # print('spin is',spin)

    if valid:
            ax_lig3D,r_emsg = lig_load(axial_ligs[0],licores) # load ligand
            if r_emsg:
                    emsg += r_emsg
            ax_lig3D.convert2mol3D() ## mol3D representation of ligand
            eq_lig3D,r_emsg = lig_load(equitorial_ligs[0],licores) # load ligand
            if r_emsg:
                    emsg += r_emsg
            eq_lig3D.convert2mol3D() ## mol3D representation of ligand
            if ax_tcat:
                    ax_lig3D.cat = ax_tcat
                    print('custom ax tcat ',ax_tcat)
            if eq_tcat:
                    ax_lig3D.cat = ax_tcat
                    print('custom eq tcat ',eq_tcat)

    if valid:
        valid,ax_type = get_con_at_type(ax_lig3D,ax_lig3D.cat)
    if valid:
        valid,eq_type = get_con_at_type(eq_lig3D,eq_lig3D.cat)

    if valid:
        eq_ki = get_truncated_kier(eq_lig3D,eq_lig3D.cat)
        ax_ki = get_truncated_kier(ax_lig3D,ax_lig3D.cat)
#        print('ax_ki',eq_ki)
#        print('eq_ki',ax_ki)
        eq_EN = get_lig_EN(eq_lig3D,eq_lig3D.cat)
        ax_EN = get_lig_EN(ax_lig3D,ax_lig3D.cat)
 #       print('ax_EN',eq_EN)
 #       print('eq_EN',ax_EN)
        eq_bo = get_bond_order(eq_lig3D.OBmol.OBMol)
        ax_bo = get_bond_order(ax_lig3D.OBmol.OBMol)
        eq_charge = eq_lig3D.OBmol.charge
        ax_charge = ax_lig3D.OBmol.charge

         
        ## preprocess:
        sum_delen  = (2.0)*ax_EN + (4.0)*eq_EN
        if abs(eq_EN) > abs(ax_EN):
                max_delen = eq_EN
        else:
                max_delen = ax_EN
        alpha = 0.2 # default for B3LYP
  #      print('ax_bo',ax_bo)
 #       print('eq_bo',eq_bo)
  #      print('ax_charge',ax_charge)
  #      print('eq_charge',eq_charge)

   #     print('sum_delen',sum_delen)
   #     print('max_delen',max_delen)
   #     print('ax_type',ax_type)
   #     print('eq_type',eq_type)


    if valid:
        sfd = get_sfd()
    #    print(sfd)
        ### scale for ANN by normalizing all values
        alpha = (alpha - sfd['alpha'][0])/sfd['alpha'][1]
        ox = (ox - sfd['ox'][0])/sfd['ox'][1]
        eq_dent = (eq_dent - sfd['eq_dent'][0])/sfd['eq_dent'][1]
        ax_dent = (eq_dent - sfd['ax_dent'][0])/sfd['ax_dent'][1]
        sum_delen = (sum_delen - sfd['sum_delen'][0])/sfd['sum_delen'][1]
        max_delen = (max_delen - sfd['max_delen'][0])/sfd['max_delen'][1]
        eq_bo = (eq_bo - sfd['eq_bo'][0])/sfd['eq_bo'][1]
        ax_bo = (ax_bo - sfd['ax_bo'][0])/sfd['ax_bo'][1]
        eq_ki = (eq_ki - sfd['eq_ki'][0])/sfd['eq_ki'][1]
        ax_ki = (ax_ki - sfd['ax_ki'][0])/sfd['ax_ki'][1]
        nn_excitation = [0,0,0,0,0, # metals co/cr/fe/mn/ni                 #1-4
                   ox,alpha,eq_charge,ax_charge, #ox/alpha/eqlig charge/axlig charge #5-8
                   eq_dent,ax_dent,# eq_dent/ax_dent/ #9-10
                   0,0,0,0, # axlig_connect: Cl,N,O,S #10 -14
                   0,0,0,0, # eqliq_connect: Cl,N,O,S #14-18
                   sum_delen,max_delen, #mdelen, maxdelen #23-24
                   ax_bo,eq_bo, #axlig_bo, eqliq_bo #19-20
                   ax_ki,eq_ki]#axlig_ki, eqliq_kii #21-22
   # print(nn_excitation)
   # print('\n')
    ### discrete variable encodings
    if valid:
        valid,nn_excitation = metal_corrector(nn_excitation,this_metal)
   # print('metal_cor',valid)
    #
    if valid:
        valid,nn_excitation = ax_lig_corrector(nn_excitation,ax_type)
    #print('ax_cor',valid)
    #
    if valid:
        valid,nn_excitation = eq_lig_corrector(nn_excitation,eq_type)
    #print('eq_cor',valid)
    #

    if valid:
        print("*******************************************************************")
        print("************** ANN is engaged and advising on spin ****************")
        print("************** and metal-ligand bond distancess    ****************")
        print("*******************************************************************")
        if high_spin:
            print('You have selected a high-spin state, s = ' + str(spin))
        else:
            print('You have selected a low-spin state, s = ' + str(spin))
        ## engage ANN
        delta = 0 
        delta = get_splitting(nn_excitation)
        if delta[0] < 0 and not high_spin:
            if abs(delta[0]) > 5:
                print('warning, ANN predicts a high spin ground state for this complex')
            else:
                print('warning, ANN predicts a near degenerate ground state for this complex')
        if delta[0] >= 0 and high_spin:
            if abs(delta[0]) > 5:
                print('warning, ANN predicts a low spin ground state for this complex')
            else:
                    print('warning, ANN predicts a near degenerate ground state for this complex')
        print("ANN predicts a spin splitting (HS - LS) of " + str(delta[0]) + ' kcal/mol')
        r = 0
        if not high_spin:
            r = get_ls_dist(nn_excitation)
        else:
            r = get_hs_dist(nn_excitation)

        print('ANN bond length is predicted to be: '+str(r) + ' angstrom')
        print("*******************************************************************")


    return valid,r


def ax_lig_corrector(excitation,con_atom_type):
    ax_lig_index_dictionary = {'Cl':11,'F':11,'N':12,'O':13,'S':14,'C':14}
    ## C is the basic value
    try:
        if not con_atom_type == "C":
            excitation[ax_lig_index_dictionary[con_atom_type]] = 1
        else:
            excitation[ax_lig_index_dictionary[con_atom_type]] = 0
        valid = True
    except:
       valid = False
    return valid,excitation

def eq_lig_corrector(excitation,con_atom_type):
    eq_lig_index_dictionary = {'Cl':15,'F':15,'N':16,'O':17,'S':18,'C':0}
    try:
        if not con_atom_type == "C":
            excitation[eq_lig_index_dictionary[con_atom_type]] = 1
        else:
            excitation[eq_lig_index_dictionary[con_atom_type]] = 0
 
        valid = True
    except:
       valid = False
    return valid,excitation



def metal_corrector(excitation,metal):
    metal_index_dictionary = {'co':0,'cr':1,'fe':2,'mn':3,'ni':4}
    try:
        excitation[metal_index_dictionary[metal]] = 1
        valid = True
    except:
       valid = False
    return valid,excitation
#n = network_builder([25,50,51],"nn_split")

def get_sfd():
    sfd = {"split_energy":[-54.19,142.71],
           "slope":[-174.20,161.58],
           "ls_min":[1.8146,0.6910],
           "hs_min":[1.8882,0.6956],
           "ox":[2,1],
           "alpha":[0,0.3],
           "ax_charge":[-2,2],
           "eq_charge":[-2,2],
           "ax_dent":[1,1],
           "eq_dent":[1,3],
           "sum_delen":[-5.34,12.54],
           "max_delen":[-0.89, 2.09],
           "ax_bo":[0,3],
           "eq_bo":[0.00,3],
           "ax_ki":[0.00, 4.29],
           "eq_ki":[0.00,6.96]}
    return sfd

def spin_classify(metal,spin,ox):
    metal_spin_dictionary = {'co':{2:4,3:5},
                              'cr':{2:5,3:4},
                              'fe':{2:5,3:6},
                              'mn':{2:6,3:5},
                              'ni':{2:3}}
    high_spin = False
        #print(metal_spin_dictionary[metal],ox)
        #print('spin checking: ',str(spin),str(metal_spin_dictionary[metal][ox]))
    if (int(spin) >= int(metal_spin_dictionary[metal][ox])):
        high_spin = True
    return high_spin

def get_splitting(excitation):
    sfd = get_sfd()
    delta = simple_splitting_ann(excitation)
    delta = delta*sfd['split_energy'][1] + sfd['split_energy'][0] 
    return delta
def get_ls_dist(excitation):
    sfd = get_sfd()
    r = simple_ls_ann(excitation)
    r = r*sfd['ls_min'][1] + sfd['ls_min'][0] 
    return r

def get_hs_dist(excitation):
    sfd = get_sfd()
    r = simple_hs_ann(excitation)
    r = r*sfd['hs_min'][1] + sfd['hs_min'][0] 
    return r

