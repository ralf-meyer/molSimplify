def check_ligands(dent_list,occ_list):
    n_ligs = occ_list

    for i in range(0,n_ligs):
        pass
    pass




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

def examine_inputs(args,occs,ligs):
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
    
    print("\n\n")
    print('Here comes occs')
    print(occs)
    print('Ligands')
    print(ligs)
    print("\n\n")

#    this_
#this_key


