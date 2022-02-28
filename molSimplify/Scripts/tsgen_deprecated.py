# flake8: noqa
#### This section comes from the structgen_one section of the code, must 
#### be reformatted to reproduce structgen behavior
def generate_ts(args):
    if (args.tsgen):
        substrate = []
        for i in args.substrate:
            substrate.append(i)
        subcatoms = False
        subcatoms_all = ['0']
        if args.subcatoms:
            if 'hat' in [str(subcatom).lower() for subcatom in args.subcatoms]:
                sub, emsg, subcatoms = substr_load(substrate[0], 0, subcatoms)
                sub.convert2mol3D()
                subcatoms_all = sub.getHs()
                if args.debug:
                    print(subcatoms_all)
            if 'epo' in [str(subcatom).lower() for subcatom in args.subcatoms]:
                sub, emsg, subcatoms = substr_load(substrate[0], 0, subcatoms)
                sub.convert2mol3D()
                subcatoms_all = sub.getC2Cs()
                if args.debug:
                    print(subcatoms_all)
            else:
                subcatoms = []
                for i in args.subcatoms:
                    subcatoms.append(int(i))
        mlig = args.mlig
        mligcatoms = args.mligcatoms
        if not args.subcatoms and len([i for i in list(subcores.keys()) if i.subname == substrate[0]]) > 0:
            nruns = len([i for i in list(subcores.keys())
                         if i.subname == substrate[0]]) > 0
        else:
            nruns = len(subcatoms_all)
        for runs in range(nruns):
            if 'hat' in args.subcatoms or 'epo' in args.subcatoms:
                sub_i = 0
                if type(subcatoms_all[runs]) == list:
                    subcatoms = subcatoms_all[runs]
                else:
                    subcatoms = [subcatoms_all[runs]]
                if args.debug:
                    print(('the subcatoms for run ' +
                           str(runs) + ' is ' + str(subcatoms)))
            else:
                sub_i = runs
            # if args.conformer:
            #     ncycles = 5
            # else:
            #     ncycles = 1
            # for cycles in range(ncycles):
            core3D_i = mol3D()
            core3D_i.copymol3D(core3D)
            core3D_i, complex3D, subcatoms, emsg, this_diag = msubcomplex(
                args, core3D_i, substrate, sub_i, subcatoms, mlig, subcatoms_ext, mligcatoms_ext)
            fname = name_ts_complex(rootdir, name_core, args.geometry, ligands, ligoc,
                                    substrate, subcatoms, mlig, mligcatoms, sernum, args, nconf, sanity)
            if args.debug:
                print(('fname is ' + str(fname)))
            # write xyz file
            core3D_i.writexyz(fname)
            strfiles.append(fname)
            # write report file
            this_diag.set_mol(core3D_i)
            this_diag.write_report(fname+'.report')
            # write input file from command line arguments
            getinputargs(args, fname)
        return strfiles, emsg, this_diag

# Main substrate placement routine for transition state
#  @param args Namespace of arguments
#  @oaram core3D core3D of the mcomplex
#  @param substrate a list of substrates
#  @param sub_i the index of the substrate
#  @param subcatoms the atom index in the substrate where the mcomplex is attached to
#  @param mlig the
#  @param subcatoms_ext the atom index, compounded with the number of atoms in the mcomplex,
#  in the substrate where the mcomplex is attahced to
#  @param mligcatoms_ext the atom index, compounded with the number of atoms in the mcomplex,
#  in the mcoomplex where the substrate is attahced to
#  @return mol3D of built complex, list of all mol3D ligands and core, error messages
def msubcomplex(args, core3D, substrate, sub_i, subcatoms, mlig, subcatoms_ext, mligcatoms_ext):

    globs = globalvars()
    subcores = getsubcores()
    this_diag = run_diag()
    if globs.debug:
        print(('\nGenerating TS complex with substrate and mlig:', substrate, mlig))
    if args.gui:
        args.gui.iWtxt.setText('\nGenerating complex with core:'+args.core+' and ligands: ' + ' '.join(arg.ligs)+'\n' +
                               args.gui.iWtxt.toPlainText())
        args.gui.app.processEvents()
    # import gui options
    if args.gui:
        from Classes.mWidgets import mQDialogWarn
    # initialize variables
    emsg, complex3D = False, []  # function returns
    # currently only one substrate with the occurance of one is supported.
    occs = 1
    catsmi = []     # SMILES substrates connection atoms
    smilessub = 0   # count how many smiles strings
    cats0 = []      # connection atoms for each substrate
    dentl = []      # denticity of substrates
    # tcats = []      # list of connection atoms for all substrates
    connected = []  # indices in core3D of substrate atoms connected to metal
    frozenats, freezeangles = [], False  # ff frozen variables
    MLoptbds = []   # list of bond lengths
    rempi = False   # remove dummy pi orbital center of mass atom
    backbatoms, batslist, bats = [], [], []  # backbond atoms variables
    bondls_m3D, bangles_m3D = [], []  # molecular parameters
    charge = 0
    if args.calccharge:
        charge = core3D.charge
    elif args.charge:
        charge = args.charge[0]
    ox = args.oxstate if args.oxstate else False
    spin = args.spin if args.spin else False
    # load substrate
    sub, subcatoms, emsg = substr_load(substrate[0], sub_i, subcatoms)
    sub.convert2mol3D()
    # # remove dummy atoms in core3D
    for atom_i, atom in enumerate(core3D.atoms):
        asym = atom.sym
        if asym == 'X':
            core3D.deleteatom(atom_i)
    # obtain rxn types
    if len(sub.grps) > 0:
        rxn_types = sub.grps
    else:
        rxn_types = ['inter'] * len(subcatoms)
    # obtaining molecular parameters by rxn types
    bondlss = []
    bangless = []
    bdihedralss = []
    if len(rxn_types) != len(mligcatoms_ext):
        mligcatoms_ext = mligcatoms_ext * len(rxn_types)
    for rxn_type_i, rxn_type in enumerate(rxn_types):
        #
        # load bond data
        #
        # load dicts
        fsr_dict = loaddata_ts('/Data/ML_FSR_for_' + rxn_type + '.dat')
        MLS_dict = loaddata_ts('/Data/MLS_FSR_for_' + rxn_type + '.dat')
        # get indexes
        midxes = core3D.findMetal()  # list of metal indexes
        mligcatom_ext = mligcatoms_ext[rxn_type_i]
        mcoord = core3D.getAtom(midxes[0]).coords()
        d0 = 10
        for i, fidx in enumerate(core3D.getBondedAtoms(mligcatom_ext)):
            fcoord = core3D.getAtom(fidx).coords()
            d = distance(fcoord, fcoord)
            if d < d0:
                d0 = d
                mligfidx_ext = fidx
        subcatom = int(subcatoms[rxn_type_i])
        bnum0 = 4
        d0 = 10
        subfidx = 0
        for i, fidx in enumerate(sub.getBondedAtoms(subcatom)):
            if fidx in subcatoms:
                subfidx = fidx
                break
            else:
                if len(subcatoms) > rxn_type_i + 1:
                    d = distance(sub.getAtomCoords(fidx),
                                 sub.getAtomCoords(subcatoms[rxn_type_i+1]))
                    if d < d0:
                        d0 = d
                        subfidx = fidx
                else:
                    bnum = len(sub.getBondedAtoms(fidx))
                    if bnum < bnum0:
                        bnum0 = bnum
                        subfidx = fidx
        #! TODO this part should be contracted by using the same method
        # get atom symbols
        msym = core3D.getAtom(midxes[0]).sym
        mligcatomsym = core3D.getAtom(mligcatom_ext).sym
        mligfidxsym = core3D.getAtom(mligfidx_ext).sym
        subcatomsym = sub.getAtom(subcatom).sym
        subfidxsym = sub.getAtom(subfidx).sym
        # get valence electron counts
        mligcatomve = int(amassdict[mligcatomsym][3])
        mligfidxve = int(amassdict[mligfidxsym][3])
        subcatomve = int(amassdict[subcatomsym][3])
        subfidxve = int(amassdict[subfidxsym][3])
        # get the number of bonded atoms
        mligcatombnum = len(core3D.getBondedAtoms(mligcatom_ext))
        subcatombnum = len(sub.getBondedAtoms(subcatom))
        mligfidxbnum = len(core3D.getBondedAtoms(mligfidx_ext))
        subfidxbnum = len(sub.getBondedAtoms(subfidx))
        # get remaining valency
        octet = 8 if mligcatomsym != 'H' else 2
        mligcatomval = (octet - mligcatomve - mligcatombnum)
        octet = 8 if subcatomsym != 'H' else 2
        subcatomval = (octet - subcatomve - subcatombnum)
        if mligfidxsym in ['V', 'Cr', 'Mn', 'Fe', 'Co']:
            octet = 18
            multiplier = 2
        elif mligfidxsym in ['Ni', 'Cu']:
            octet = 14
            multiplier = 2
        elif mligfidxsym != 'H':
            octet = 8
            multiplier = 1
        else:
            octet = 2
            multiplier = 1
        mligfidxval = (octet - mligfidxve - mligfidxbnum * multiplier)
        octet = 8 if subfidxsym != 'H' else 2
        subfidxval = (octet - subfidxve - subfidxbnum)
        # preping for bondls
        # get covalent radius
        mr = float(amassdict[msym][2])
        mligcatomr = float(amassdict[mligcatomsym][2])
        mligfidxr = float(amassdict[mligfidxsym][2])
        subcatomr = float(amassdict[subcatomsym][2])
        subfidxr = float(amassdict[subfidxsym][2])
        # get sum of cov radii
        sumr_m_mlig = mr + mligcatomr
        sumr_mlig_mligfidx = mligcatomr + mligfidxr
        sumr_mlig_sub = mligcatomr + subcatomr
        sumr_sub_subfidx = subcatomr + subfidxr
        # keys for bondls, in the order of bondl_core3D, bondl_m3D, and bondl_sub3D
        keys = [[mligcatomsym + str(mligcatomval), mligfidxsym + str(mligfidxval), sumr_mlig_mligfidx],
                [mligcatomsym + str(mligcatomval), subcatomsym +
                 str(subcatomval), sumr_mlig_sub],
                [subcatomsym + str(subcatomval), subfidxsym + str(subfidxval), sumr_sub_subfidx]]
        # obtain bondls and bangles
        # default ratio
        fsr = 1
        bondls = []
        for i, key in enumerate(keys):
            bondl_m3D = False
            sym1 = key[0]
            sym2 = key[1]
            sumr = key[2]
            if i == 1:
                bondl_m3D = True
            fsr, exact_match = get_ts_fsr_database(
                sym1, sym2, fsr_dict, ox, spin, bondl_m3D=bondl_m3D)
            bondl = fsr * sumr
            bondls.append(bondl)
        bondlss.append(bondls)
        print(('keys are ' + str(keys)))
        print(('bondlss are ' + str(bondlss)))
        # prepping for bangles
        # keys for bangles, in  the order of bangle_m3D, bangle_m3Dsub, bangle_core3D, and bangle_sub
        # keys = [[mligcatomsym, mligcatomval], [subcatomsym, subcatomval], [mligfidxsym, mligfidxval],
        #         [subfidxsym, subfidxval]]
        keys = [[mligcatomsym, str(mligcatomval)], [subcatomsym, str(subcatomval)], [mligfidxsym, str(mligfidxval)],
                [subfidxsym, str(subfidxval)]]
        bdihedrals = []
        bangles = []
        for key in keys:
            sym = key[0]
            val = key[1]
            bangle, bdihedral, exact_match = get_ts_MLSangle_database(
                sym, val, MLS_dict, ox, spin)
            bangles.append(bangle)
            bdihedrals.append(bdihedral)
        bangless.append(bangles)
        bdihedralss.append(bdihedrals)
        # bangless = [[120, 180, 180, 109]]
        print(('keys are ' + str(keys)))
        print(('bangless are ' + str(bangless)))
        print(('bdihedralss are ' + str(bdihedralss)))
    # determine if the sub in lig
    sub_in_lig = False
    if args.lig:
        if (args.substrate[0].lower() in [lig_i.lower() for lig_i in args.lig]):
            sub_in_lig = True
    if not sub_in_lig:
        # freeze core
        for i in range(0, core3D.natoms):
            frozenats.append(i)
        # freeze key atoms in substrate
        if args.debug:
            print(('subcatoms after init_substrate is ' + str(subcatoms)))
        if len(subcatoms) > 1:
            for subcatom in subcatoms:
                if isinstance(subcatom, int):
                    frozenats.append(core3D.natoms + subcatom)
        else:
            frozenats.append(core3D.natoms + int(subcatoms[0]))
            for bondedat in sub.getBondedAtoms(int(subcatoms[0])):
                frozenats.append(core3D.natoms + bondedat)
        # compute number of connecting points required
        cpoints_required = len(bondlss)
        # load core and initialize template
        # also adjust the specified distance in the core3D and distance and angle in the connection point
        m3D, core3D, geom, backbatoms, coord, corerefatoms = init_mcomplex_template(args, core3D, cpoints_required,
                                                                                    mligcatoms_ext, bondlss, bangless,
                                                                                    bdihedralss)

        totsub = 0  # total number of substrates added
        subsused = 0
        # initialize ligand
        sub3D, rempi, subcatoms, subpiatoms = init_substrate(
            args, sub, subcatoms, bondlss, bangless, bdihedralss)
        for subcatom in subcatoms:
            frozenats.append(core3D.natoms + subcatom)
            for batidx in sub3D.getBondedAtoms(subcatom):
                frozenats.append(core3D.natoms + batidx)
        if emsg:
            return False, emsg
        for j in range(0, occs):
            denticity = sub.denticity
            # print('denticity is ' + str(denticity))
            if not(substrate == 'x' or substrate == 'X'):
                # add atoms to connected atoms list
                # catoms = sub.cat # connection atoms
                # initatoms = core3D.natoms # initial number of atoms in core3D
                # for at in catoms:
                #     connected.append(initatoms+at)
                # initialize variables
                # metal c oordinates in backbone
                mcoords = core3D.getAtom(midxes[0]).coords()
                atom0, r0, r1, r2, r3 = 0, mcoords, 0, 0, 0  # initialize variables
                coreref = corerefatoms.getAtom(totsub)
                # connecting point in backbone to align ligand to
                # batoms = get_batoms(args,batslist,subsused)
                cpoints = m3D.atoms[m3D.natoms-denticity:]
                # attach ligand depending on the denticity
                # optimize geometry by minimizing steric effects
                if (denticity == 1):
                    sub3D = align_sub(args, cpoints, core3D, coreref, sub3D, subcatoms, mligcatoms_ext[0], bangless,
                                      rempi, subpiatoms)
                elif (denticity == 2):
                    # for mligcatom_ext in mligcatoms_ext:
                    #     if mligcatom_ext not in core3D.findMetal():
                    #         break
                    # print('cpoints are ' + str([cpoint.coords() for cpoint in cpoints]))
                    # print('mligcatom_ext is ' + str(mligcatom_ext))
                    # print('subcatoms are ' + str(subcatoms))
                    sub3D = align_dent2_sub(
                        args, cpoints, m3D, core3D, mligcatom_ext, sub3D, subcatoms, bangless)
                auxm = mol3D()
                auxm.copymol3D(sub3D)
                complex3D.append(auxm)

                if args.cdxml:
                    # combine molecules into a crude msubcomplex
                    core3D_ = mol3D()
                    core3D_.copymol3D(core3D)
                    liglist, ligdents, ligcons = ligand_breakdown(core3D_)
                    ax_ligand_list, eq_ligand_list, ax_natoms_list, eq_natoms_list, ax_con_int_list, eq_con_int_list, \
                        ax_con_list, eq_con_list, built_ligand_list = ligand_assign(
                            core3D_, liglist, ligdents, ligcons)
                    if globs.testTF():
                        # new RACs-ANN
                        from molSimplify.Scripts.tf_nn_prep import invoke_ANNs_from_mol3d
                    # result_dict, {"ls_bl": r_ls, "hs_bl": r_hs, "split": split, "distance": split_dist}
                    result_dict = invoke_ANNs_from_mol3d(
                        core3D_, oxidation_state=3)
                    # for ds_i, ds in enumerate(ds1):
                    #     print('The euclidean distance between the built M-L' + str(ds_i) + ' and the training data is ' + str(min(ds)))
                    eu_dist = result_dict['distance']
                    print(
                        ('The euclidean distance between the built rOH and the training data is ' + str(eu_dist)))
                    midx = core3D_.findMetal()[0]
                    bondls = result_dict['ls_bl'][0]
                    for bondl_i, bondl in enumerate(bondls):
                        if bondl_i < 2:
                            idxes = ax_con_list[bondl_i]
                            for idx in idxes:
                                core3D_.BCM(idx, midx, bondl)
                                # print('Adjusting the distance between %s and %s to %s' %(idx, midx, bondl))
                        else:
                            idxes = eq_con_list[bondl_i-2]
                            for idx in idxes:
                                core3D_.BCM(idx, midx, bondl)
                                # print('Adjusting the distance between %s and %s to %s' % (idx, midx, bondl))
                    # # refined msubcomplex build
                    m3D, core3D, geom, backbatoms, coord, corerefatoms = init_mcomplex_template(args, core3D_,
                                                                                                cpoints_required,
                                                                                                mligcatoms_ext,
                                                                                                bondlss,
                                                                                                bangless)
                    # m3D.writexyz(fpath + '/init_mcomplex.xyz')

                    if (denticity == 1):
                        sub3D = align_sub(args, cpoints, core3D, coreref, sub3D, subcatoms, mligcatoms_ext[0], bangless,
                                          rempi, subpiatoms)

                    elif (denticity == 2):
                        # for mligcatom_ext in mligcatoms_ext:
                        #     if mligcatom_ext not in core3D_.findMetal():
                        #         break
                        # print('mligcatom_ext is ' + str(mligcatom_ext))
                        sub3D = align_dent2_sub(
                            args, cpoints, m3D, core3D, mligcatom_ext, sub3D, subcatoms, bangless)
                auxm = mol3D()
                auxm.copymol3D(sub3D)
                complex3D.append(auxm)
                # if 'a' not in sub.ffopt.lower():
                #     for latdix in range(0,sub3D.natoms):
                #         frozenats.append(latdix+core3D.natoms)
                # combine molecules into a crude msubcomplex
                core3D = core3D.combine(sub3D)
                core3D.convert2mol3D()
                # core3D.writexyz(fpath + '/msubcomplex.xyz')
                # remove dummy cm atom if requested
                # if rempi:
                #     core3D.deleteatom(core3D.natoms-1)
                if args.calccharge:
                    core3D.charge += sub3D.charge
                totsub += denticity
                subsused += 1
        # perform FF optimization if requested
        if 'a' in args.ffoption:
            print('Performing final FF constrained opt for the metal-substrate complex')
            connected = []
            freezeangles = []
            core3D, enc = ffopt(args.ff, core3D, connected, 1, frozenats, freezeangles, MLoptbds, 'Adaptive',
                                args.debug)
            # core3D.writexyz(fpath + '/msubcomplex_ff.xyz')
    else:
        core3D = align_intra_sub(
            args, core3D, subcatoms_ext[0], mligcatoms_ext[0], bondlss, bangless)
        # reserved for future conformer search development
        # for cycle in range(5):
        #     if 'a' in args.ffoption:
        #         print('Performing five conformer search: ' + str(cycle) + ' out of 5.')
        #         connected = []
        #         freezeangles = []
        #         core3D,enc = conformer_search(args.ff,core3D,1,frozenats,freezeangles,MLoptbds,'Adaptive',args.debug)
        #         # core3D,enc = ffopt(args.ff,core3D,connected,1,frozenats,freezeangles,MLoptbds,'Adaptive',args.debug)
    return core3D, complex3D, subcatoms, emsg, this_diag

# Initializes core and template mol3Ds and properties
#  @param args Namespace of arguments
#  @param cpoints_required Number of connecting points required
#  @return mol3D of core, template, geometry, backbone atoms, coordination number, core reference atom index
def init_mcomplex_template(args, core3D, cpoints_required, mligcatoms_ext, bondlss, bangless, bdihedralss):
    # initialize core and template
    m3D = mol3D()
    # container for ordered list of core reference atoms
    corerefatoms = mol3D()
    # geometry load flag
    geom = False
    cpoint_nums = []
    backbatoms = []
    # bring mligcatoms_ext that are corresponding to core to the front
    i = 0
    # print(mligcatoms_ext)
    for idx, mligcatom_ext in enumerate(mligcatoms_ext):
        mligcatomsym = core3D.getAtom(mligcatom_ext).sym
        if mligcatomsym in [core.lower() for core in args.core]:
            mligcatoms_ext[idx], mligcatoms_ext[i] = mligcatoms_ext[i], mligcatoms_ext[idx]
            bondlss[idx], bondlss[i] = bondlss[i], bondlss[idx]
            bangless[idx], bangless[i] = bangless[i], bangless[idx]
            i += 1
    # check mligcatoms
    d_ref = 100
    m3D.copymol3D(core3D)
    for i in range(cpoints_required):
        bondl_core3D, bondl_m3D, bondl_sub = bondlss[i][0], bondlss[i][1], bondlss[i][2]
        bangle_m3D, bangle_core3D, bangle_sub = bangless[i][0], bangless[i][2], bangless[i][3]
        bdihedral_m3D, bdihedral_core3D, bdihedrel_sub = bdihedralss[
            i][0], bdihedralss[i][2], bdihedralss[i][3]
        mligcatom_ext = mligcatoms_ext[i]
        mligcatomcoords = m3D.getAtom(mligcatom_ext).coords()
        mligcatomsym = m3D.getAtom(mligcatom_ext).sym
        midx = m3D.findMetal()[0]
        mcoord = m3D.getAtom(midx).coords()
        # adjust mlig
        # obtain the anchoring atom
        for idx in m3D.getBondedAtoms(mligcatom_ext):
            coord = m3D.getAtom(idx).coords()
            d = distance(mcoord, coord)
            if d < d_ref:
                d_ref = d
                mligcatom_ext_anchor = idx
                mliganchorcoords = coord
        # adjust M-L distance for the incoming substrate
        # if (args.core not in m3D.getAtom(mligcatom_ext).sym):
        core3D.BCM(mligcatom_ext, mligcatom_ext_anchor, bondl_core3D)
        m3D.BCM(mligcatom_ext, mligcatom_ext_anchor, bondl_core3D)
        # check natoms in mlig
        if args.mlig:
            lig, emsg = lig_load(args.mlig[0])
            lig.convert2mol3D()
            natoms_lig = lig.natoms
            # adjust M-{L} angle for the incoming substrate
            if natoms_lig > 1:
                core3D.ACM(mligcatom_ext, mligcatom_ext_anchor,
                           midx, bangle_core3D)
                m3D.ACM(mligcatom_ext, mligcatom_ext_anchor,
                        midx, bangle_core3D)
        # if mligcatomsym in [core.lower() for core in args.core]:
        #
        if i == 0:
            # print('mligcatom_ext_anchor is ' + str(mligcatom_ext_anchor))
            # print('bangle_m3D is ' + str(bangle_m3D))
            # cpoint = getconnectiongivenphi(m3D, mligcatom_ext, mligcatom_ext_anchor, bondl_m3D, bangle_m3D)
            cpoint = getconnectiongivenangles(
                m3D, mligcatom_ext, mligcatom_ext_anchor, bondl_m3D, (bangle_m3D, bdihedral_m3D))
            # cpoint = getconnection(m3D, mligcatom_ext, bondl_m3D)
            # store core reference atom
            conatom3D = atom3D(m3D.getAtom(mligcatom_ext).sym,
                               m3D.getAtom(mligcatom_ext).coords())
            corerefatoms.addAtom(conatom3D)
            if args.debug:
                print((corerefatoms.getAtom(0).symbol()))
            # initiate dummy atom
            dummy_atom = atom3D(Sym='X', xyz=cpoint)
            m3D.addAtom(dummy_atom)
            cpoint_nums.append(m3D.natoms - 1)
            # nums = m3D.findAtomsbySymbol('I')
            backbatoms = getbackbcombsall(cpoint_nums)
            mligcatom_old = mligcatom_ext
        # if the second cpoint is on the same mligcatom as the first one
        elif mligcatom_ext == mligcatom_old:
            cpoint = getconnectiongivenr(
                m3D, mligcatom_ext, mligcatom_ext_anchor, bondl_m3D, bangle_m3D, bondl_sub)
            # print('bondl_sub is ' + str(bondl_sub))
            # store core reference atom
            conatom3D = atom3D(core3D.getAtom(mligcatom_ext).sym,
                               core3D.getAtom(mligcatom_ext).coords())
            corerefatoms.addAtom(conatom3D)
            if args.debug:
                print((corerefatoms.getAtom(0).symbol()))
            # corerefatoms.append(ccatoms[i])N
            # add connecting points to template
            m3D.addAtom(atom3D(Sym='X', xyz=cpoint))
            cpoint_nums.append(m3D.natoms - 1)
            # except IndexError:
            #     pass
            backbatoms = getbackbcombsall(cpoint_nums)
        else:
            cpoint = getconnectiongivenangles(
                m3D, mligcatom_ext, mligcatom_ext_anchor, bondl_m3D, (bangle_m3D, bdihedral_m3D))
            # cpoint = getconnectiongivenphi(m3D, mligcatom_ext, mligcatom_ext_anchor, bondl_m3D, bangle_m3D)
            # store core reference atom
            conatom3D = atom3D(core3D.getAtom(mligcatom_ext).sym,
                               core3D.getAtom(mligcatom_ext).coords())
            corerefatoms.addAtom(conatom3D)
            if args.debug:
                print((corerefatoms.getAtom(0).symbol()))
            # initiate dummy atom
            dummy_atom = atom3D(Sym='X', xyz=cpoint)
            m3D.addAtom(dummy_atom)
            cpoint_nums.append(m3D.natoms - 1)
            # obtain rotation axis
            mliganchorcoords = m3D.getAtomCoords(mligcatom_ext_anchor)
            mligcatoms_ext_anchor2 = [idx for idx in m3D.getBondedAtoms(
                mligcatom_ext_anchor) if idx != mligcatom_ext]
            for idx in mligcatoms_ext_anchor2:
                coord = m3D.getAtom(idx).coords()
                d = distance(mcoord, coord)
                mliganchor2coords = coord
                if d < d_ref:
                    d_ref = d
                    mliganchor2coords = coord
            r_mliganchor_mliganchor2 = vecdiff(
                mliganchorcoords, mliganchor2coords)
            # rotate the [mligcatom_ext, mligcatom_ext_anchor] bond to bring the dummy atoms closer
            d0 = 100
            refidx = sorted(m3D.findAtomsbySymbol('X'))[0]
            refcoords = m3D.getAtomCoords(refidx)
            num_rotation = 0
            m3D_ = mol3D()
            m3D_.copymol3D(m3D)
            theta = 5
            while num_rotation < 72:
                num_rotation += 1
                m3D_.ACM_axis(mligcatom_ext, mligcatom_ext_anchor,
                              r_mliganchor_mliganchor2, theta)
                xidx = sorted(m3D_.findAtomsbySymbol('X'))[-1]
                xcoords = m3D_.getAtomCoords(xidx)
                d = distance(refcoords, xcoords)
                if d < d0:
                    d0 = d
                    m3D = mol3D()
                    m3D.copymol3D(m3D_)
    # set charge from oxidation state if desired
    if args.calccharge:
        if args.oxstate:
            if args.oxstate in list(romans.keys()):
                core3D.charge = int(romans[args.oxstate])
            else:
                core3D.charge = int(args.oxstate)
    # remove X atoms from m3D to generate core3D
    core3D = mol3D()
    core3D.copymol3D(m3D)
    for atidx in range(core3D.natoms)[::-1]:
        asym = core3D.getAtom(atidx).sym
        if asym == 'X':
            core3D.deleteatom(atidx)

    return m3D, core3D, geom, backbatoms, coord, corerefatoms

# Initializes substrate 3D geometry and properties
#  @param args Namespace of arguments
#  @param sub mol3D of ligand
#  @param subcatoms list of connecting atoms in the substrate
#  @return mol3D of ligand, flag for pi-coordination, pi-coordinating atoms
def init_substrate(args, sub, subcatoms, bondlss, bangless, bdihedralss):
    rempi = False
    idx3 = ''
    sub.convert2mol3D()
    sub3D = mol3D()
    sub3D.copymol3D(sub)
    # sub3D.createMolecularGraph(False)
    for rxn_type_i in range(len(bondlss)):
        bondl_m3D = bondlss[rxn_type_i][1]
        bondl_sub = bondlss[rxn_type_i][2]
        bangle_sub = bangless[rxn_type_i][3]
        bdihedral_sub = bdihedralss[rxn_type_i][3]
        # if SMILES string, copy connecting atoms list to mol3D properties
        # if not sub.cat and tcats[i]:
        #     if 'c' in tcats[i]:
        #         lig.cat = [lig.natoms]
        #     else:
        #         lig.cat = tcats[i]
        # change name
        # substrate bond manipulation to prepare for activation
        anchor_atom_idx = int(subcatoms[rxn_type_i])
        ref_to_sub = False
        moved_atom_idx = ''
        for atidx in sub3D.getBondedAtoms(anchor_atom_idx):
            moved_atom_coords = sub3D.getAtomCoords(atidx)
            # check to see if there is another connection atom
            distss = []
            if atidx in subcatoms:
                moved_atom_idx = atidx
                break
            else:
                ref_to_sub = True
                for subcatom in subcatoms:
                    subcoords = sub3D.getAtomCoords(subcatom)
                    d = distance(moved_atom_coords, subcoords)
                    distss.append((d, atidx))
        if ref_to_sub:
            distss.sort()
            moved_atom_idx = distss[0][1]
        elif not moved_atom_idx:
            moved_atom_idx = sub3D.getBondedAtoms(anchor_atom_idx)[0]
        sub3D.BCM(moved_atom_idx, anchor_atom_idx, bondl_sub)
        # angle adjustment
        if sub3D.natoms > 2 and len(sub3D.getBondedAtoms(moved_atom_idx)) > 1 and moved_atom_idx not in subcatoms:
            idx3 = [idx for idx in sub3D.getBondedAtoms(
                moved_atom_idx) if idx != anchor_atom_idx][0]
            sub3D.ACM(anchor_atom_idx, moved_atom_idx, idx3, bangle_sub)
        # # check for pi-coordinating substrate
        subpiatoms = []
        if 'pi' in sub.cat:
            if args.debug:
                print('substrate is identified as a pi-type substrate')
            subcatoms = []
            for k in sub.cat[:-1]:
                if isinstance(k, int):
                    subcatoms.append(k)
                # sub3Dpiatoms.addAtom(sub3D.getAtom(k))
                # sub3Dpiatoms.addAtom(sub3D.getAtom(k))
            # subpiatoms = sub.cat[:-1]
            # sub3D.addAtom(atom3D('C',sub3Dpiatoms.centermass()))
            # if args.debug:
            #     sub3D.printxyz()
            # sub3D.cat = [sub3D.natoms-1]
            rempi = True
        cpoint = getconnectiongivenangles(
            sub3D, anchor_atom_idx, moved_atom_idx, bondl_m3D, (bangle_sub, bdihedral_sub))
        dummy = atom3D(Sym='X', xyz=cpoint)
        sub3D.addAtom(dummy)
        xidxes = sorted(sub3D.findAtomsbySymbol('X'))
        refidx = xidxes[0]
        refcoords = sub3D.getAtomCoords(refidx)
        if len(xidxes) > 1 and idx3:
            anchor_coords = sub3D.getAtomCoords(anchor_atom_idx)
            moved_coords = sub3D.getAtomCoords(moved_atom_idx)
            idx3_coords = sub3D.getAtomCoords(idx3)
            u0 = vecdiff(anchor_coords, moved_coords)
            u1 = vecdiff(moved_coords, idx3_coords)
            # u =
            d0 = 10
            thetas = [0, 0]
            theta1 = 5
            for theta0 in range(0, 360, 5):
                num_rotation = 0
                sub3D_ = mol3D()
                sub3D_.copymol3D(sub3D)
                sub3D_.ACM_axis(anchor_atom_idx, moved_atom_idx, u0, theta0)
                while num_rotation < 72:
                    num_rotation += 1
                    sub3D_.ACM_axis(anchor_atom_idx,
                                    moved_atom_idx, u1, theta1)
                    pidx = xidxes[1]
                    pcoords = sub3D_.getAtomCoords(pidx)
                    d = distance(refcoords, pcoords)
                    # planar = checkplanar(anchor_coords, moved_coords, pcoords, refcoords)
                    if d < d0:
                        d0 = d
                        thetas[0] = theta0
                        thetas[1] = theta1 * num_rotation
            # print('distance is ' + str(d0))
            # print('theta is ' + str(thetas))
            sub3D.ACM_axis(anchor_atom_idx, moved_atom_idx, u0, thetas[0])
            sub3D.ACM_axis(anchor_atom_idx, moved_atom_idx, u1, thetas[1])

    return sub3D, rempi, subcatoms, subpiatoms

# Finds the optimum attachment point for an atom/group to a central atom given the desired bond length and bond angle
#
#  Objective function maximizes the minimum distance between attachment point and other groups bonded to the central atom
#  @param core mol3D of core
#  @param cidx Core connecting atom index
#  @param refidx idx of the atom bonded to the core connecting atom
#  @param BL Optimal core-ligand bond length
#  @return Coordinates of optimum attachment point
def getconnectiongivenangles(core, cidx, refidx, BL, BAs):
    # ncore = core.natoms
    # groups = core.getBondedAtoms(cidx)
    cpoint = []
    ccoords = core.getAtom(cidx).coords()
    refcoords = core.getAtomCoords(refidx)
    anchoridxes = [idx for idx in core.getBondedAtoms(cidx) if idx != refidx]
    # print('cidx is ' + str(cidx))
    # print('refidx is ' + str(refidx))
    # print('anchoridxes are ' + str(anchoridxes))
    r_c_ref = vecdiff(refcoords, ccoords)
    pcoords = getPointu(ccoords, BL, r_c_ref)
    # print('BL is ' + str(BL))
    # theta, u = rotation_params(refcoords, ccoords, pcoords)
    phi = BAs[0] / 180. * pi
    phi_ref = BAs[0]
    theta_ref = BAs[1]
    objopt = 0
    if anchoridxes:
        # print('theta_ref is ' + str(theta_ref))
        anchorcoords = core.getAtomCoords(anchoridxes[0])
        theta, norm = rotation_params(refcoords, ccoords, anchorcoords)
        # brute force search
        objopt = 0
        distss = []
        # print('phi_ref is ' + str(phi_ref))
        # print('theta_ref is ' + str(theta_ref))
        #! TODO this double for loop needs to go
        for itheta in range(0, 360, 10):
            theta = itheta / 180. * pi
            for iphi in range(0, 360, 10):
                phi = iphi / 180. * pi
                pcoords_ = PointTranslateSph(
                    ccoords, r_c_ref, [BL, theta, phi])
                angle = vecangle(norm, r_c_ref)
                angle_phi = vecangle(
                    r_c_ref, [pcoords_[i] - ccoords[i] for i in range(len(ccoords))])
                if angle < 80:
                    angle_theta = theta_ref
                else:
                    angle_theta = vecangle(
                        norm, [pcoords_[i] - ccoords[i] for i in range(len(ccoords))])
                # print('angle is ' + str(angle))
                if abs(angle_theta - theta_ref) < 5 and abs(angle_phi - phi_ref) < 5:
                    dists = []
                    xidxes = core.findAtomsbySymbol('X')
                    # whether the second X should be placed near the first one
                    if xidxes:
                        for xidx in xidxes:
                            dists.append(-1 *
                                         distance(core.getAtomCoords(xidx), pcoords_))
                    else:
                        for ig in range(core.natoms):
                            if ig != cidx:
                                dists.append(
                                    distance(core.getAtomCoords(ig), pcoords_))
                    distss.append((min(dists), pcoords_))
                    # print('distss are ' + str(distss))
                    # print('angle_phi is ' + str(angle_phi) + ', angle_theta is ' + str(angle_theta))
                core3D = mol3D()
                core3D.copymol3D(core)
                atom = atom3D(Sym='X', xyz=pcoords_)
                core3D.addAtom(atom)
        # print('distss are ' + str(sorted(distss)[-1]))
        cpoint = sorted(distss)[-1][1]
    else:
        # print('theta_ref is ' + str(theta_ref))
        distss = []
        for itheta in range(1, 359, 1):
            theta = itheta / 180. * pi
            for iphi in range(0, 360, 10):
                phi = iphi / 180. * pi
                pcoords_ = PointTranslateSph(
                    ccoords, pcoords, [BL, theta, phi])
                angle_phi = vecangle(
                    r_c_ref, [pcoords_[i] - ccoords[i] for i in range(len(ccoords))])
                if abs(angle_phi - phi_ref) < 5:
                    dists = []
                    for ig in range(core.natoms):
                        if ig != cidx:
                            dists.append(
                                distance(core.getAtomCoords(ig), pcoords_))
                    distss.append((min(dists), pcoords_))
        cpoint = sorted(distss)[-1][1]

    return cpoint

# Finds the optimum attachment point for an atom/group to a central atom given the desired bond length and bond angle
#
#  Objective function maximizes the minimum distance between attachment point and other groups bonded to the central atom
#  @param core mol3D of core
#  @param cidx Core connecting atom index
#  @param refidx idx of the atom bonded to the core connecting atom
#  @param BL Optimal core-ligand bond length
#  @return Coordinates of optimum attachment point
def getconnectiongivenr(core, cidx, refidx, BL, BA, r):
    midx = core.findMetal()[0]
    ccoords = core.getAtom(cidx).coords()
    if core.getAtom(cidx).sym in [core.getAtom(midx).sym]:
        refcoords = [0, 0, 0]
        for fidx in core.getBondedAtoms(cidx):
            nfidx = len(core.getBondedAtoms(cidx))
            refcoords[0] += core.getAtom(fidx).coords()[0]/nfidx
            refcoords[1] += core.getAtom(fidx).coords()[1]/nfidx
            refcoords[2] += core.getAtom(fidx).coords()[2]/nfidx
    else:
        refcoords = core.getAtom(refidx).coords()
    phi = float(BA/180.*pi)
    # brute force search
    cpoint = []
    P0coords = core.atoms[-1].coords()
    P = PointTranslateSphgivenr(ccoords, refcoords, [BL, phi], P0coords, r)
    cpoint = P
    return cpoint

# Aligns a substrate's activated bond along the metal-connecting atom axis
#  @param corerefcoords Core reference coordinates
#  @param sub3D mol3D of substrate
#  @param atom0 Substrate connecting atom index
#  @param core3D mol3D of partially built complex
#  @return mol3D of aligned substrate
def align_sub_firstbond(args, corerefcoords, sub3D, atom0, core3D, bangle_m3Dsub):
    # rotate to align center of symmetry
    # globs = globalvars()
    r0 = corerefcoords
    r1 = sub3D.getAtom(atom0).coords()
    sub3Db = mol3D()
    sub3Db.copymol3D(sub3D)
    # auxmol = mol3D()
    # for at in sub3D.getBondedAtomsSmart(atom0):
    #     auxmol.addAtom(sub3D.getAtom(at))
    at = sub3D.getBondedAtomsSmart(atom0)[0]
    r2 = sub3D.getAtomCoords(at)
    theta, u = rotation_params(r0, r1, r2)
    # rotate around axis and get both images
    sub3D = rotate_around_axis(sub3D, r1, u, bangle_m3Dsub-(180-theta))
    if args.debug:
        print(('bangle m3Dsub is ' + str(bangle_m3Dsub)))
        print(('theta is ' + str(theta)))
    sub3D_aligned = mol3D()
    sub3D_aligned.copymol3D(sub3D)

    return sub3D_aligned

# Aligns a monodentate substrate to core connecting atom coordinates
#  @param args Namespace of arguments
#  @param cpoint atom3D containing backbone connecting point
#  @param core3D mol3D of partially built complex
#  @param coreref atom3D of core reference atom
#  @param substrate Name of substrate for dictionary lookup
#  @param sub3D mol3D of substrate
#  @param catoms List of substrate connecting atom indices
#  @param rempi Flag for pi-coordinating substrate
#  @param subpiatoms List of pi-coordinating atom indices in substrate
#  @param MLb Custom M-L bond length (if any)
#  @param ANN_flag Flag for ANN activation
#  @param ANN_bondl ANN-predicted M-L bond length
#  @param this_diag ANN diagnostic object
#  @param MLbonds M-L bond dictionary
#  @param MLoptbds List of final M-L bond lengths
#  @param i Ligand serial number
#  @param EnableAutoLinearBend Flag for enabling automatic bending of linear ligands (e.g. superoxo)
#  @return mol3D of aligned ligand, updated list of M-L bond lengths
def align_sub(args, cpoints, core3D, coreref, sub3D, subcatoms, mligcatoms_ext, bangless, rempi, subpiatoms):
    bangle_m3Dsub = bangless[0][1]
# ,ANN_flag=False,ANN_bondl=[],this_diag=0,MLbonds=dict(),MLoptbds=[],i=0,EnableAutoLinearBend=True):
    corerefcoords = coreref.coords()
    # connection atom in sub3D
    atom0 = int(subcatoms[0])
    subcoords = sub3D.getAtom(atom0)
    if args.debug:
        print(('atom0 is ' + str(atom0)))
    # print substrate coordinates before translation
    if args.debug:
        print(corerefcoords)
        print((cpoints[0].coords()))
        # print(atom0)
        # sub3D.printxyz()
    # translate ligand to overlap with backbone connecting point
    sub3D.alignmol(subcoords, cpoints[0])
    # determine bond length (database/cov rad/ANN)
    # bondl = get_MLdist(args,lig3D,atom0,ligand,coreref,MLb,i,ANN_flag,ANN_bondl,this_diag,MLbonds)
    # MLoptbds = []
    # bondl = 2
    # MLoptbds.append(bondl)
    # align ligand to correct M-L distance
    u = vecdiff(cpoints[0].coords(), corerefcoords)
    # sub3D = aligntoaxis2(sub3D, cpoint.coords(), corerefcoords, u, bondl)
    # print('cpoint is ' + str(cpoints[0].coords()) + '. corerefcoords is ' + str(corerefcoords) + '.')
    sub3D = aligntoaxis(sub3D, cpoints[0].coords(), corerefcoords, u)
    if args.debug:
        print(('length of subpiatoms is ' + str(len(subpiatoms))))
    if rempi and len(subcatoms) > 1:
        # align linear (non-arom.) pi-coordinating ligand
        sub3D = align_linear_pi_sub(
            core3D, mligcatoms_ext, sub3D, atom0, subcatoms, bangle_m3Dsub)
        if args.debug:
            print(
                ('aligning a linear pi ligand at a bangle_m3Dsub of ' + str(bangle_m3Dsub)))
    elif sub3D.natoms > 1:
        # align ligand center of symmetry
        sub3D = align_sub_firstbond(
            args, corerefcoords, sub3D, atom0, core3D, bangle_m3Dsub)
        print(('bangle_m3Dsub is ' + str(bangle_m3Dsub)))
        if sub3D.natoms > 2:
            # check for linear molecule and align
            sub3D = check_rotate_linear_lig(corerefcoords, sub3D, atom0)
            # check for symmetric molecule
            sub3D = check_rotate_symm_lig(corerefcoords, sub3D, atom0, core3D)
        # rotate around M-L axis to minimize steric repulsion
        # sub3D = rotate_MLaxis_minimize_steric_ts(mligcatoms_ext,sub3D,atom0,core3D)
        # rotate around L-Sub axis to minimize steric repulsion
        sub3D = rotate_MLaxis_minimize_steric(
            corerefcoords, sub3D, atom0, core3D)
    # rotate the substrate around the M-mlig axis to reduce repulsion
    mcoords = [core3D.getAtom(idx).coords() for idx in core3D.getBondedAtoms(
        mligcatoms_ext) if idx in core3D.findMetal()][0]
    rmref = vecdiff(mcoords, corerefcoords)
    sub3D_aligned = rotate_MLaxis_minimize_steric_ts(
        sub3D, subcoords, core3D, corerefcoords, rmref)
    # if args.debug:
    #     sub3D_aligned.printxyz()

    return sub3D_aligned

# Aligns an intramolecular monodentate substrate to core connecting atom coordinates
#  @param args Namespace of arguments
#  @param cpoint atom3D containing backbone connecting point
#  @param core3D mol3D of partially built complex
#  @param coreref atom3D of core reference atom
#  @param substrate Name of substrate for dictionary lookup
#  @param sub3D mol3D of substrate
#  @param catoms List of substrate connecting atom indices
#  @param rempi Flag for pi-coordinating substrate
#  @param subpiatoms List of pi-coordinating atom indices in substrate
#  @param MLb Custom M-L bond length (if any)
# @param ANN_flag Flag for ANN activation
# @param ANN_bondl ANN-predicted M-L bond length
#  @param this_diag ANN diagnostic object
#  @param MLbonds M-L bond dictionary
#  @param MLoptbds List of final M-L bond lengths
#  @param i Ligand serial number
#  @param EnableAutoLinearBend Flag for enabling automatic bending of linear ligands (e.g. superoxo)
#  @return mol3D of aligned ligand, updated list of M-L bond lengths
def align_intra_sub(args, core3D, subcatoms_ext, mligcatoms_ext, bondl_core3D, bondl_m3D, bondl_sub, bangless):
    bangle_m3D = bangless[0][0]
    bangle_m3Dsub = bangless[0][1]
    reacting_at_list = []
    reacting_at_list.append(subcatoms_ext)
    reacting_at_list.append(mligcatoms_ext)
    midx_list = core3D.findMetal()
    constr_list = []
    for midx in midx_list:
        constr_list.append(midx)
        for fidx in core3D.getBondedAtoms(midx):
            if fidx not in reacting_at_list:
                constr_list.append(fidx)
    sidx = core3D.getBondedAtoms(subcatoms_ext)[0]
    for i in core3D.getBondedAtoms(subcatoms_ext):
        if i in midx_list:
            sidx = i
    # bondl_sub_i = distance(core3D.getAtom(subcatoms_ext[0]).coords(),core3D.getAtom(sidx).coords())
    # bondl_m3D_i = distance(core3D.getAtom(subcatoms_ext[0]).coords(),core3D.getAtom(mligcatoms_ext[0]).coords())
    # bondl_core3D_i = distance(core3D.getAtom(mligcatoms_ext[0]).coords(),core3D.getAtom(midx_list[0]).coords())
    # bangle_m3D_i,u = rotation_params(core3D.getAtom(midx_list[0]).coords(),core3D.getAtom(mligcatoms_ext[0]).coords(),core3D.getAtom(subcatoms_ext[0]).coords())
    core3D.convert2OBMol()
    OBMol = core3D.OBMol
    # for bond in openbabel.OBMolBondIter(OBMol):
    #     idx1 = bond.GetBeginAtomIdx()
    #     idx2 = bond.GetEndAtomIdx()
    #     # if (idx1-1 in constr_list and idx2-1 in reacting_at_list) or (idx2-1 in constr_list and idx1-1 in reacting_at_list):
    #     #     bond.SetBO(0)
    #     if (idx1-1 == 31) and (idx2-1 == 33):
    #         bond.SetBO(2)

    # core3D.OBMol = OBMol
    # core3D.convert2mol3D()
    # core3D.convert2OBMol()
    # OBMol = core3D.OBMol
    ff = openbabel.OBForceField.FindForceField('mmff94')
    # for a_m3D in numpy.arange(bangle_m3D_i,bangle_m3D,-10).tolist():
    #     for b_sub in numpy.arange(bondl_sub_i,bondl_sub,-0.1).tolist():
    #         for b_m3D in numpy.arange(bondl_m3D_i,bondl_m3D,-0.1).tolist():
    #             for b_core3D in numpy.arange(bondl_core3D_i,bondl_core3D,-0.1).tolist():
    #                 constr = openbabel.OBFFConstraints()
    #                 for idx in constr_list:
    #                     constr.AddAtomConstraint(idx+1)
    #                 for subidx in subcatoms_ext:
    #                     constr.AddDistanceConstraint(subidx+1,sidx+1,b_sub) # bondl_sub
    #                     for mligcatoms in mligcatoms_ext:
    #                         constr.AddDistanceConstraint(subidx+1,mligcatoms+1,b_m3D) # bondl_m3D
    #                     for midx in midx_list:
    #                         constr.AddDistanceConstraint(midx+1,mligcatoms+1,b_core3D) # bondl_core3D
    #                         constr.AddAngleConstraint(midx+1,mligcatoms_ext[0]+1,subidx+1,a_m3D) # bangle_m3D
    constr = openbabel.OBFFConstraints()
    for idx in constr_list:
        constr.AddAtomConstraint(idx+1)
    subidx = subcatoms_ext
    constr.AddDistanceConstraint(subidx+1, sidx+1, bondl_sub)  # bondl_sub
    if args.debug:
        print(('Distance constraint between %s and %s is %s' %
               (subidx, sidx, bondl_sub)))
    mligcatoms = mligcatoms_ext
    constr.AddDistanceConstraint(
        subidx+1, mligcatoms+1, bondl_m3D)  # bondl_m3D
    if args.debug:
        print(('Distance constraint between %s and %s is %s' %
               (subidx, mligcatoms, bondl_m3D)))
    atnos = []
    for midx in midx_list:
        constr.AddDistanceConstraint(
            midx+1, mligcatoms+1, bondl_core3D)  # bondl_core3D
        if args.debug:
            print(('Distance constraint between %s and %s is %s' %
                   (midx, mligcatoms, bondl_core3D)))
        constr.AddAngleConstraint(
            midx+1, mligcatoms_ext+1, subidx+1, bangle_m3D)  # bangle_m3D
        if args.debug:
            print(('Angle constraint among %s, %s, and %s is %s' %
                   (midx, mligcatoms_ext, subidx, bangle_m3D)))
        atno = OBMol.GetAtom(midx+1).GetAtomicNum()
        atnos.append(atno)
        OBMol.GetAtom(midx+1).SetAtomicNum(14)
        # constr.AddDistanceConstraint(1,22,2.202)
    # printing obatom type
    # for obatom in openbabel.OBMolAtomIter(OBMol):
    #     print(obatom.GetType())
    s = ff.Setup(OBMol, constr)
    ff.SteepestDescent(500)
    for i, midx in enumerate(midx_list):
        atno = atnos[i]
        OBMol.GetAtom(midx+1).SetAtomicNum(atno)
    ff.GetCoordinates(OBMol)
    # printing obatom type
    # for obatom in openbabel.OBMolAtomIter(OBMol):
    #     print(obatom.GetType())
    core3D.OBMol = OBMol
    core3D.convert2mol3D()

    return core3D

# Aligns a bidentate substrate to core connecting atom coordinates
#  @param args Namespace of arguments
#  @param cpoint atom3D containing backbone connecting point
#  @param m3D mol3D of backbone template
#  @param core3D mol3D of partially built complex
#  @param mligcatom
#  @param sub3D mol3D of substrate
#  @param subcatoms List of substrate connecting atom indices
#  @param bangle_m3Dsub
#  @return mol3D of aligned ligand, updated lists of frozen atoms and M-L bond lengths
def align_dent2_sub(args, cpoints, m3D, core3D, mligcatom, sub3D, subcatoms, bangless):
    bangle_m3Dsub = bangless[0][1]
    corerefcoords = m3D.getAtom(mligcatom).coords()
    r0 = corerefcoords
    # get cis conformer by rotating rotatable bonds
    #lig3D = find_rotate_rotatable_bond(lig3D,catoms)
    # connection atom
    atom0 = subcatoms[0]
    atom1 = subcatoms[1]
    subcoords = sub3D.getAtom(atom1).coords()
    # translate ligand to match first connecting atom to backbone connecting point
    sub3D.alignmol(sub3D.getAtom(atom0), cpoints[0])
    # rotate the substrate to align with the cpoints
    r1 = sub3D.getAtom(atom0).coords()  # atom0 coord
    r2 = sub3D.getAtom(atom1).coords()  # atom1 coord
    rx2 = cpoints[1].coords()  # second cpoint coordxxx
    theta, u = rotation_params(r2, r1, rx2)
    sub3D = rotate_around_axis(sub3D, r1, u, 180-theta)
    # m3D_ = mol3D()
    # # for atom in m3D.atoms:
    # #     if atom.sym != 'X':
    # #         m3D_.addAtom(atom)
    # m3D_.copymol3D(m3D)
    # m3D_.copymol3D(sub3D)
    # rotate the sub3D along the two coordinating atoms to maximize the overlap between sub3D-Xs and the mcomplex
    r1 = sub3D.getAtom(atom0).coords()  # atom0 coord
    r2 = sub3D.getAtom(atom1).coords()  # atom1 coord
    xidxes_m3D = sorted(m3D.findAtomsbySymbol('X'))
    xidxes_sub3D = sorted(sub3D.findAtomsbySymbol('X'))
    refidx = m3D.getBondedAtoms(xidxes_m3D[0])[0]
    refcoords = m3D.getAtomCoords(refidx)
    xcoords = sub3D.getAtomCoords(xidxes_sub3D[0])
    theta, u = rotation_params(refcoords, r1, xcoords)
    u = vecdiff(r1, r2)
    # print('theta is ' + str(theta))
    if theta < 90:
        theta = 180 - theta
    sub3D = rotate_around_axis(sub3D, r1, u, theta)

    sub3D_aligned = sub3D
    # sub3D_aligned = rotate_MLaxis_minimize_steric_ts(sub3D, subcoords, m3D, corerefcoords, rmref)
    # r1 = sub3D_aligned.getAtom(atom0).coords()  # atom0 coord
    # r2 = sub3D_aligned.getAtom(atom1).coords()  # atom1 coord
    # r12 = vecdiff(r2, r1)
    # subcoords = sub3D_aligned.getAtom(atom1).coords()
    # sub3D_aligned = rotate_MLaxis_minimize_steric_ts(sub3D_aligned, subcoords, m3D, subcoords, r12)

    return sub3D_aligned

# Loads ts M-L-S angle from database and reports if compound is in DB
#  @param sym atomic symbol
#  @param sym valency
#  @param fsr_dict formal shortness ratio dictionary
#  @param ox oxidation state of the metal
#  @param spin spin multiplicity of the metal
#  @return fsr, flag for exact DB match
def get_ts_MLSangle_database(sym, val, MLS_dict, ox=False, spin=False):
    # check for roman letters in oxstate
    if ox:  # if defined put oxstate in keys
        if ox in list(romans.keys()):
            oxs = romans[ox]
        else:
            oxs = ox
    else:
        oxs = '-'
    # check for spin multiplicity
    spin = spin if spin else '-'
    keys = []
    syms = [[sym, val], [val, sym]]
    oss = [[oxs, spin], [oxs, '-'], ['-', spin], ['-', '-']]
    for sym in syms:
        for os in oss:
            key = sym + os
            keys.append(tuple(key))
    found = False
    exact_match = False
    # search for data
    for key in keys:
        if (key in list(MLS_dict.keys())):  # if exact key in dictionary
            bangle = float(MLS_dict[key][0])
            bdihedral = float(MLS_dict[key][1])
            found = True
            if (key == ((sym, val, oxs, spin))):  # exact match
                exact_match = True
            break
    if not found:  # last resort covalent radii
        bangle = 120
        bdihedral = 90
    return bangle, bdihedral, exact_match


# Loads ts M-L fsr from database and reports if compound is in DB
#  @param sym1 atomic symbol 1
#  @param sym2 atomic symbol 2
#  @param fsr_dict formal shortness ratio dictionary
#  @param ox oxidation state of the metal
#  @param spin spin multiplicity of the metal
#  @return fsr, flag for exact DB match
def get_ts_fsr_database(sym1, sym2, fsr_dict, ox=False, spin=False, bondl_m3D=False):
    # check for roman letters in oxstate
    if ox:  # if defined put oxstate in keys
        if ox in list(romans.keys()):
            oxs = romans[ox]
        else:
            oxs = ox
    else:
        oxs = '-'
    # check for spin multiplicity
    spin = spin if spin else '-'
    keys = []
    syms = [[sym1, sym2], [sym2, sym1]]
    oss = [[oxs, spin], [oxs, '-'], ['-', spin], ['-', '-']]
    for sym in syms:
        for os in oss:
            key = sym + os
            keys.append(tuple(key))
    found = False
    exact_match = False
    # search for data
    for key in keys:
        if (key in list(fsr_dict.keys())):  # if exact key in dictionary
            if bondl_m3D:
                fsr = float(fsr_dict[key][1])
            else:
                fsr = float(fsr_dict[key][0])
            found = True
            if (key == ((sym1, sym2, oxs, spin))):  # exact match
                exact_match = True
            break
    if not found:  # last resort covalent radii
        fsr = 1
    return fsr, exact_match

# Aligns a linear pi substrate's connecting point to the metal-substrate axis
#  @param corerefcoords Core reference coordinates
#  @param sub3D mol3D of suband
#  @param atom0 substrate connecting atom index
#  @param subpiatoms List of substrate pi-connecting atom indices
#  @return mol3D of aligned substrate
def align_linear_pi_sub(core3D, mligcatoms_ext, sub3D, atom0, subcatoms, bangle_m3D):
    # align the norm of the pi bond to the L-S vector at atom0 (one of the C in the pi bond)
    ratom0 = sub3D.getAtom(atom0).coords()
    try:
        ratom1 = [sub3D.getAtom(i).coords() for i in sub3D.getBondedAtoms(
            atom0) if i in subcatoms][0]
    except:
        ratom1 = [sub3D.getAtom(i).coords()
                  for i in subcatoms if i is not atom0][0]
    ratom2 = [sub3D.getAtom(i).coords() for i in sub3D.getBondedAtoms(
        atom0) if i not in subcatoms][0]
    theta102, uatom102 = rotation_params(ratom1, ratom0, ratom2)
    rL = core3D.getAtom(mligcatoms_ext).coords()
    rLS = vecdiff(rL, ratom0)
    thetaL0u102, uL0u102 = rotation_params(rL, ratom0, uatom102)
    thetaL0u102 = vecangle(rLS, uatom102)  # aka u1
    # u2 = cross(uatom102,rLS)
    sub3D_aligned = mol3D()
    sub3D_aligned.copymol3D(sub3D)
    # sub3D_aligned = rotate_around_axis(sub3D_aligned, ratom0, u2, 90-theta_uatom102_rLS)
    # rotate sub3D such that the substrate lies perpendicular to the complex
    sub3D_aligned = rotate_around_axis(
        sub3D_aligned, ratom0, uL0u102, -thetaL0u102)
    # rotate sub3D such that the X-C-M angle matches bangle_m3Dsub
    try:
        ratom1 = [sub3D.getAtom(i).coords() for i in sub3D.getBondedAtoms(
            atom0) if i in subcatoms][0]
    except:
        ratom1 = [sub3D.getAtom(i).coords()
                  for i in subcatoms if i is not atom0][0]
    ratom2 = [sub3D.getAtom(i).coords() for i in sub3D.getBondedAtoms(
        atom0) if i not in subcatoms][0]
    theta10L, u10L = rotation_params(ratom1, ratom0, rL)
    sub3D_aligned = rotate_around_axis(
        sub3D_aligned, ratom0, u10L, bangle_m3D+theta10L)
    # agjust the angle among L-C-Ca
    try:
        ratom1 = [sub3D.getAtom(i).coords() for i in sub3D.getBondedAtoms(
            atom0) if i in subcatoms][0]
    except:
        ratom1 = [sub3D.getAtom(i).coords()
                  for i in subcatoms if i is not atom0][0]
    ratom2 = [sub3D.getAtom(i).coords() for i in sub3D.getBondedAtoms(
        atom0) if i not in subcatoms][0]
    # objfuncopt = 90
    # thetaopt = 0
    # for theta in range(0,360,1):
    #     sub3D_aligned = rotate_around_axis(sub3D_aligned,ratom0,rLS, theta)
    #     #objfunc = abs(vecangle(vecdiff(sub3D_aligned.getAtom(atom0).coords(),corerefcoords),vecdiff(sub3D_aligned.getAtom(ligpiatoms[0]).coords(),sub3D_aligned.getAtom(ligpiatoms[1]).coords()))-90)
    #     objfunc = abs(distance(sub3D_aligned.getAtom(ligpiatoms[0]).coords(),corerefcoords) - distance(sub3D_aligned.getAtom(ligpiatoms[1]).coords(),corerefcoords))
    #     if objfunc < objfuncopt:
    #         thetaopt = theta
    #         objfuncopt = objfunc
    #         sub3Dopt = mol3D() # lig3Dopt = sub3D_aligned DOES NOT WORK!!!
    #         sub3Dopt.copymol3D(sub3D_aligned)
    # sub3D_aligned = rotate_around_axis(sub3D_aligned, ratom0, rLS, 30)
    # reduce the steric repulsion between the mcomplex and the substrate by rotating away the ther C in the pi bond
    sub3D_aligned = rotate_MLaxis_minimize_steric(
        rL, sub3D_aligned, atom0, core3D)

    return sub3D_aligned
    
# Rotates aligned ligand about M-L axis to minimize steric clashes with rest of complex
#  @param mol mol3D of the molecule to be rotated
#  @param molatcoords the coordinates of the atom in the rotated molecule
#  @param refmol mol3D of the molecule of reference
#  @param refcoords the coordinates of the atom in the reference molecule
#  @param u the vector of rotation axis
#  @return mol3D of rotated ligand
def rotate_MLaxis_minimize_steric_ts(mol, coords, refmol, refcoords, u):
    dist0 = 0
    dist = 0
    theta0 = 0
    while theta0 < 360:
        dists = []
        mol_ = rotate_around_axis(mol, refcoords, u, theta0)
        for atom in mol_.atoms:
            coords_ = atom.coords()
            if coords_ != coords or mol_.natoms == 1:
                for atom_ref in refmol.atoms:
                    refcoords_ = atom_ref.coords()
                    if refcoords_ != refcoords:
                        dist = distance(coords_, refcoords_)
                        dists.append(dist)
        dist = min(dists)
        if dist > dist0:
            dist0 = dist
            mol_aligned = mol3D()
            mol_aligned.copymol3D(mol_)
        theta0 += 5

    return mol
