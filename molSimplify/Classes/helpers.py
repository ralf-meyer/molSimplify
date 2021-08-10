#  @file helpers.py
#  Contains some useful helper functions for classes.
#
#  Written by HJK Group
#
#  Dpt of Chemical Engineering, MIT

from molSimplify.Classes.AA3D import AA3D
from molSimplify.Classes.mol3D import mol3D
from molSimplify.Classes.atom3D import atom3D

def read_atom(line):
    """ Reads a line of a pdb into an atom dictionary.

    Parameters
    ----------
        line : str
            a line of a pdb file

    Returns
    -------
        atom_dict : dictionary
            dictionary containing attributes of atom when reading the pdb
    """
            
    labels = ['Type', 'SerialNum', 'Name', 'AltLoc', 'ResName', 'ChainID',
              'ResSeq', 'ICode', 'X', 'Y', 'Z', 'Occupancy', 'TempFactor',
              'Element', 'Charge']
    data = [line[0:6].strip(), int(line[6:11]), line[12:16].strip(),
            line[16:17].strip(), line[17:20].strip(), line[21:22].strip(),
            int(line[22:26]), line[26:27].strip(), float(line[30:38]),
            float(line[38:46]), float(line[46:54]), float(line[54:60]),
            float(line[60:66]), line[70:78].strip(), line[78:80].strip()]
    atom_dict = dict(zip(labels, data))
    atom_dict['Element'] = atom_dict['Element'][0] + atom_dict['Element'][1:].lower()
    return atom_dict

def makeMol(a_dict, mols, conf, chains, prev_a_dict, bonds, aa=True):
    """ Creates an AA3D from a_dicts and adds it to the appropriate places.

    Parameters
    ----------
        a_dict : dictionary
            created from the read_atom method
        mols : dictionary
            keyed by chain and location, valued by AA3Ds or mol3Ds
        conf : list
            contains amino acid conformations
        chains : dictionary
            contains amino acid chains
        prev_a_dict : dictionary
            the a_dict from the previous run
        bonds : dictionary
            contains bonding info
        aa : boolean
            True if mols consists of AA3Ds, False if consists of mol3Ds

    Returns
    -------
        mols, conf, chains, prev_a_dict as updated
    """
    loc = a_dict['AltLoc']
    ploc = prev_a_dict["AltLoc"]
    m = 0
    if (a_dict['ChainID'], a_dict['ResSeq']) not in mols.keys():
        if aa:
            m = AA3D(a_dict['ResName'], a_dict['ChainID'], a_dict['ResSeq'],
                     a_dict['Occupancy'], loc)
        else:
            m = mol3D(a_dict['ResName'], loc)
        mols[(a_dict['ChainID'], a_dict['ResSeq'])] = [m]
    if loc != '' and ploc != '':
        l = ord(ploc)
        if loc > chr(l) and len(mols[(a_dict['ChainID'], a_dict['ResSeq'])]) == l-64:
            if aa:
                m = AA3D(a_dict['ResName'], a_dict['ChainID'],
                         a_dict['ResSeq'], a_dict['Occupancy'], loc)
            else:
                m = mol3D(a_dict['ResName'], loc)
            mols[(a_dict['ChainID'], a_dict['ResSeq'])].append(m)
    prev_a_dict = a_dict
    if a_dict['ChainID'] not in chains.keys():
        chains[a_dict['ChainID']] = [] # initialize key of chain dictionary
    if m != 0 and int(float(a_dict['Occupancy'])) != 1 and m not in conf:
        conf.append(m)
    if m != 0 and m not in chains[a_dict['ChainID']] and m not in conf:
        chains[a_dict['ChainID']].append(m)
    if loc == '' or loc == "A" or ploc == '':
        m = mols[(a_dict['ChainID'], a_dict['ResSeq'])][0]
    elif (l-64) < len(mols[(a_dict['ChainID'], a_dict['ResSeq'])]):
        m = mols[(a_dict['ChainID'], a_dict['ResSeq'])][l-64]
    else:
        m = mols[(a_dict['ChainID'], a_dict['ResSeq'])][-1]
    atom = atom3D(Sym=a_dict['Element'], xyz=[a_dict['X'], a_dict['Y'],
                                              a_dict['Z']],
                  Tfactor=a_dict['TempFactor'], occup=a_dict['Occupancy'],
                  greek=a_dict['Name'])
    m.addAtom(atom, a_dict['SerialNum']) # terminal Os may be missing
    if aa:
        m.setBonds()
        bonds.update(m.bonds)
        if m.prev != None:
            bonds[m.n].add(m.prev.c)
        if m.next != None:
            bonds[m.c].add(m.next.n)
    return atom, mols, conf, chains, prev_a_dict, bonds
