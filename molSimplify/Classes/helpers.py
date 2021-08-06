#  @file helpers.py
#  Contains some useful helper functions for classes.
#
#  Written by HJK Group
#
#  Dpt of Chemical Engineering, MIT

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

def make_aa(a_dict, aas, prev_a_dict):
    """ Creates an AA3D from a_dicts and adds it to the appropriate places.

    Parameters
    ----------
        a_dict : dictionary
            created from the read_atom method
        aas : dictionary
            keyed by chain and location, valued by AA3Ds
        conf : list
            contains amino acid conformations
        prev_a_dict : dictionary
            the a_dict from the previous run

    Returns
    -------
        a : AA3D
            AA3D created from the a_dict
    """
    loc = a_dict['AltLoc']
    if loc != '' and prev_a_dict["AltLoc"] != '':
        if loc > chr(l) and len(aas[(a_dict['ChainID'],
                                     a_dict['ResSeq'])]) == l-64:
            a = AA3D(a_dict['ResName'], a_dict['ChainID'],
                 a_dict['ResSeq'], a_dict['Occupancy'], loc)
            aas[(a_dict['ChainID'], a_dict['ResSeq'])].append(a)
    prev_a_dict = a_dict
    if loc != '':
        l = ord(a_dict["AltLoc"])
    if (a_dict['ChainID'], a_dict['ResSeq']) not in aas.keys():
        a = AA3D(a_dict['ResName'], a_dict['ChainID'],
             a_dict['ResSeq'], a_dict['Occupancy'], loc)
        aas[(a_dict['ChainID'], a_dict['ResSeq'])] = [a]
    if a_dict['ChainID'] not in chains.keys():
        chains[a_dict['ChainID']] = [] # initialize key of chain dictionary
    if int(float(a_dict['Occupancy'])) != 1 and a not in conf:
        conf.append(a)
    if a not in chains[a_dict['ChainID']] and a not in conf:
        chains[a_dict['ChainID']].append(a)
    if a_dict["AltLoc"] == '' or a_dict["AltLoc"] == "A":
        a = aas[(a_dict['ChainID'], a_dict['ResSeq'])][0]
    elif (l-65) < len(aas[(a_dict['ChainID'], a_dict['ResSeq'])]):
        a = aas[(a_dict['ChainID'], a_dict['ResSeq'])][l-65]
    else:
        a = aas[(a_dict['ChainID'], a_dict['ResSeq'])][-1]
    return a
