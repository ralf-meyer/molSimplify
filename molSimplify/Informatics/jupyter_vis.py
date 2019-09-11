"""
Py3Dmol install: (works in both Python2/3 conda environments)
conda install -c rmg py3dmol 
Some Documentation: https://pypi.org/project/py3Dmol/
3DMol.js backend: http://3dmol.csb.pitt.edu/index.html
"""
import py3Dmol
import math as m
from molSimplify.Classes.mol3D import*

def view_xyz(xyz_names,w=400,h=400,columns=2,representation='stick',labelsize=18,labels=True):
    """
    py3Dmol view atoms object(s)
    xyz_names = xyz files that will be rendered in a tiled format in jupyter (list,str)
    w = width of frame (or subframes) in pixels (int)
    h = height of frame (or subframes) in pixels (int)
    cols = number of columns in subframe (int)
    representation = how the molecule will be viewed (str)
    labelsize = size of the data label (in Points) (int)
    labels = turn labels on/off (bool)
    """
    if isinstance(xyz_names,str): # Shift to list if only string passed
        xyz_names = [xyz_names]
    if len(xyz_names) == 1:
        view_ats = py3Dmol.view(width=w,height=h)
        item = xyz_names[0]
        label = item.split('/')[-1].replace('.xyz','') # Get rid of extension for label
        mol = mol3D()
        mol.readfromxyz(item)
        metal_atom_index = mol.findMetal()
        if metal_atom_index:
            label_posits = mol.getAtomCoords(metal_atom_index[0])
        else:
            label_posits = mol.centersym() 
        view_ats.addModel(open(item,'r').read(),'xyz')
        view_ats.setStyle({representation:{'colorscheme':'Jmol'}})
        if labels:
            view_ats.addLabel("{}".format(label), {'position':{'x':'{}'.format(label_posits[0]),
                              'y':'{}'.format(label_posits[1]),'z':'{}'.format(label_posits[2])},
                              'backgroundColor':"'black'",'backgroundOpacity':'0.3',
                              'fontOpacity':'1', 'fontSize':'{}'.format(labelsize),
                              'fontColor':"white",'inFront':'true',})
        view_ats.zoomTo()
        view_ats.show()
    else:
        rows = int(m.ceil(float(len(xyz_names))/columns))
        w = w*columns
        h = h*rows 
        # print(rows,columns)
        view_ats = py3Dmol.view(width=w,height=h,linked=False,viewergrid=(rows,columns))
        x,y = 0,0 # Subframe position
        for i,item in enumerate(xyz_names):
            label = item.split('/')[-1].replace('.xyz','') # Get rid of extension for label
            mol = mol3D()
            mol.readfromxyz(item)
            metal_atom_index = mol.findMetal()
            if metal_atom_index:
                label_posits = mol.getAtomCoords(metal_atom_index[0])
            else:
                label_posits = mol.centersym()
            view_ats.addModel(open(item,'r').read(),'xyz',viewer=(x,y))
            view_ats.setStyle({representation:{'colorscheme':'Jmol'}},viewer=(x,y))
            if labels:
                view_ats.addLabel("{}".format(label), {'position':{'x':'{}'.format(label_posits[0]),
                                  'y':'{}'.format(label_posits[1]),'z':'{}'.format(label_posits[2])},
                                  'backgroundColor':"'black'",'backgroundOpacity':'0.5',
                                  'fontOpacity':'1','fontSize':'{}'.format(labelsize),
                                  'fontColor':"white",'inFront':'true',}, viewer=(x,y))
            view_ats.zoomTo(viewer=(x,y))
            if y+1 < columns:
                y+=1
            elif y+1 == columns:
                x+=1
                y=0
        view_ats.show()
