#  @file helpers.py
#  Contains some useful helper functions for classes.
#
#  Written by HJK Group
#
#  Dpt of Chemical Engineering, MIT

def read_atom(line):
    labels = ['Type', 'SerialNum', 'Name', 'AltLoc', 'ResName', 'ChainID', 'ResSeq', 'ICode', 'X', 'Y', 'Z', 'Occupancy', 'TempFactor', 'Element', 'Charge']
    data = [line[0:6].strip(), int(line[6:11]), line[12:16].strip(), line[16:17].strip(), line[17:20].strip(), line[21:22].strip(), int(line[22:26]), line[26:27].strip(), float(line[30:38]), float(line[38:46]), float(line[46:54]), float(line[54:60]), float(line[60:66]), line[70:78].strip(), line[78:80].strip()]
    atom_dict = dict(zip(labels, data))
    atom_dict['Element'] = atom_dict['Element'][0] + atom_dict['Element'][1:].lower()
    return atom_dict

"""
def pStripSpaces(l):
    #When reading a string from a pdb file, fixes some of the spacing issues
    #that result from strange behavior with whitespace.

    #Parameters
    #----------
    #    lst : a list that has just been split by whitespace from a pdb file

    #Returns
    #-------
    #    good_lst : a list with proper spacing, so it can be read properly

    if '+' in l[-1] or '-' in l[-1]: # fix charge of Sym
        l[-1] = l[-1][:(len(l[-1]) - 2)]
    if '0' in l[-1]: # fix number attached
        l[-1] = l[-1][:(len(l[-1]) - 1)]
    if len(l[-1]) == 2:# fix case
        l[-1] = l[-1][0] + l[-1][1].lower() 
    if len(l[-2]) > 6: # occupancy has too many characters
        l2 = l
        l = l2[:-2] + [l2[-2][:4], l2[-2][4:], l2[-1]]
    if len(l[1]) > 3 and len(l) != 11:
        l2 = l
        if len(l[1]) > 4 and (l[1][4] == 'A' or l[1][4] == 'B'):
            l = [l2[0], l2[1][:4], l2[1][4:]] + l2[2:]
        elif l[1][3] == 'A' or l[1][3] == 'B':
            l = [l2[0], l2[1][:3], l2[1][3:]] + l2[2:]
        elif ('1' in l[2] or len(l[2]) == 1) and l[1][3] != "'":
            l = [l2[0], l2[1][:3], l2[1][3:]] + l2[2:]
    if len(l[3]) > 1:
        digits = {'1','2','3','4','5','6','9'} # add more as needed
        l2 = l
        if len(l[2]) != 1 or l[3][1] in digits:
            l = l2[:3] + [l2[3][:1], l2[3][1:]] + l2[4:]
        elif l[-2][0] == '0' and len(l) == 11:
            l = [l2[0], l2[1]+l2[2]] + l2[3:]
        else:
            l = l2[:2] + [l2[2]+l2[3]] + l2[4:]
    if len(l[2]) == 1:
        l2 = l
        if len(l) > 11:
             l = [l2[0], l2[1]+l2[2]] + l2[3:]
        if len(l[3]) == 1 and l[4] == 'b': # 1 lcs exist
            l = l2[:2] + [l2[2]+l2[3]] + l2[4:]
    if len(l) < 11 and len(l[3]) > 1:
        l2 = l
        if l[3][1] == '1' or l[3][1] == '2':
            l = l2[:3] + [l2[3][:1], l2[3][1:]] + l2[4:]
    # fix coordinate spacing
    if '-' in l[5][1:]: 
        y = l[5]
        y = l[5].split('-')
        if len(y) > 2 and y[0] != '': # extra long string case
            l = l[:5] + [y[0], '-'+y[1], '-'+y[2]] + l[6:]
        elif y[0] != '':
            l = l[:5] + [y[0], '-'+y[1]] + l[6:]
        elif len(y) > 3: # extra long string case
            l = l[:5] + ['-'+y[1], '-'+y[2], '-'+y[3]] + l[6:]
        else:
            l = l[:5] + ['-'+y[1], '-'+y[2]] + l[6:]
    if '-' in l[6][1:]:
        y = l[6]
        y = l[6].split('-')
        if y[0] != '':
            l = l[:6] + [y[0], '-'+y[1]] + l[7:]
        else:
            l = l[:6] + ['-'+y[1], '-'+y[2]] + l[7:]
    if '-' in l[7][1:]:
        y = l[7]
        y = l[7].split('-')
        if y[0] != '':
            l = l[:7] + [y[0], '-'+y[1]] + l[8:]
        else:
            l = l[:7] + ['-'+y[1], '-'+y[2]] + l[8:]
    if "" in l:
        l.remove("")
    return l
"""