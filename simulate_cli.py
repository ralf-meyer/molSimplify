import os
import subprocess
from datetime import datetime
start = datetime.now()
#subprocess.call(["python",'main.py','-slab_gen',
#                 '-cif_path', 'Unitcells/gold.cif','-duplication_vector','[4,2,2]'])
# subprocess.call(["python",'ptcotest.py','-slab_gen','-cif_path','data/gold.cif',
#                  '-slab_size','15,15,5','-align_dist','1.89','-duplication_vector','2,2,2'])
#subprocess.call(["python",'ptcotest.py','-slab_gen','-cif_path','data/gold.cif',
#                 '-align_dist','1.89','-duplication_vector','[2,2,2]','-slab_size','[15,15,5]',
#                 '-align_distance_method','custom','miller_index','[1,1,1]','-align_method'])


#subprocess.call(["python",'main.py','-slab_gen','-cif_path','Unitcells/gold.cif',
#                 '-align_dist','2.7','-slab_size','[15,15,5]','-align_distance_method',
#                 'custom','-miller_index','[1,1,1]','-align_method','alignpair',
#                 '-place_on_slab','-target_molecule','/home/jp/Runs/copo.xyz','-surface_atom_type','Au',
#                 '-object_align','Co','-control_angle','40','-angle_control_partner','2',
#                 '-angle_surface_axis','[1,1]'])

#subprocess.call(["python",'main.py','-slab_gen','-cif_path','Unitcells/gold.cif',
#                 '-slab_size','[15,15,5]','-align_distance_method',
#                 'chemisorption','-miller_index','[1,1,1]','-align_method','alignpair',
#                 '-place_on_slab','-target_molecule','/home/jp/Runs/copo.xyz','-surface_atom_type','Au',
#                 '-object_align','Co','-control_angle','40','-angle_control_partner','2',
#                 '-angle_surface_axis','[1,1]'])

#subprocess.call(["python",'main.py','-slab_gen','-cif_path','Unitcells/gold.cif',
#                 '-slab_size','[15,15,5]',
#                 '-miller_index','[1,1,1]'])
#subprocess.call(["python",'main.py','-slab_gen','-cif_path','Unitcells/gold.cif',
#                 '-duplication_vector','[3,3,2]','-miller_index','[1,1,1]'])
#subprocess.call(["python",'main.py','-slab_gen','-cif_path','Unitcells/gold.cif',
#                 '-align_dist','1.89','-slab_size','[15,15,5]','-align_distance_method',
#                 'custom','-align_method','alignpair',
#                 '-place_on_slab','-target_molecule','/home/jp/Runs/copo.xyz','-surface_atom_type','Au',
 #                '-object_align','Co','-control_angle','90','-angle_control_partner','2',
  #               '-angle_surface_axis','[1,1]'])
#subprocess.call(["python",'main.py','-slab_gen','-cif_path','Unitcells/anatase_tio2.cif',
#                 '-slab_size','[12,12,8]','-miller_index','[1,1,0]','-place_on_slab','-target_molecule','/home/jp/Runs/mo5.xyz',
#                 '-align_method','alignpair','-surface_atom_type','Ti','-num_surface_atoms','4',
#                 '-object_align','2,3,4,5','-align_dist','2.5'])
#subprocess.call(["python",'main.py','-slab_gen','-cif_path','Unitcells/anatase_tio2.cif',
#                 '-slab_size','[12,12,8]','-place_on_slab','-miller_index','[1,0,1]','-target_molecule','/home/jp/Runs/mo5.xyz',
#                 '-align_method','alignpair','-surface_atom_type','Ti','-num_surface_atoms','4',
#                 '-object_align','2,3,4,5','-align_dist','2.5'])
#subprocess.call(["python",'main.py','-slab_gen','-cif_path','Unitcells/anatase_tio2.cif',
#                 '-slab_size','[12,12,8]','-miller_index','[1,0,1]','-target_molecule','/home/jp/Runs/mo5.xyz',
#                 '-align_method','alignpair','-surface_atom_type','Ti','-num_surface_atoms','4',
#                 '-object_align','2,3,4,5','-align_dist','2.5'])
#
#subprocess.call(["python",'-m','molSimplify.__main__','-slab_gen','-cif_path','Unitcells/anatase_tio2.cif',
#                 '-slab_size','[9,9,9]','-place_on_slab','-target_molecule','/home/jp/Runs/co.xyz',
#                 '-align_method','alignpair','-surface_atom_type','Ti','-num_surface_atoms','1',
#                 '-object_align','C','-align_dist','2.5'])
#
#subprocess.call(["python",'molSimplify/main.py','-slab_gen','-cif_path','molSimplify/Unitcells/anatase_tio2.cif',
#                 '-slab_size','[9,9,9]','-miller_index','[1,0,0]','-freeze','1'])



#
#mol_name = "sally"
#subprocess.call(["python","molSimplify/main.py","-core","Fe","-coord",'6',"-lig","acac,water","-ligocc",'2,2',
#                 '-geometry','oct','-distort','0','checkdirb','True','-ligalign','False','-calccharge','yes',
#                 '-keepHs','False','-bcharge','0','-rundir','/home/jp/Runs/mpc/\n','-jobdir temp',
#                 '-oxstate','II','-spin','1','-mopac','-name',mol_name+'\n','-qccode TeraChem'])
#subprocess.call(["python","molSimplify/main.py","-core","Fe","-coord",'4',"-lig","diaminomethyl","-ligocc",'2',

#'-geometry','sqp','-distort','0','checkdirb','True','-ligalign','False','-calccharge','yes',
#                 '-keepHs','False','-bcharge','0','-rundir','/home/jp/Runs/\n','-jobdir consti','-debug True',
#                 '-oxstate','II','-spin','5','-name',mol_name+'\n','-qccode TeraChem'])
#mol_name = "frank"
subprocess.call(["python","molSimplify/main.py","-core","Mn","-coord",'6',"-lig","porphyrin water","-ligocc",'1,2',
                '-geometry','oct','-distort','0','checkdirb','True','-ligalign','False','-calccharge','yes',
                 '-keepHs','False','-bcharge','0','-rundir','/home/jp/Runs/\n','-debug True',
                 '-oxstate','II','-spin','1'])
print('**********************************************')
subprocess.call(["python","molSimplify/main.py","-core","Fe","-coord",'6',"-lig","porphyrin water","-ligocc",'1,2',
                '-geometry','oct','-distort','0','checkdirb','True','-ligalign','False','-calccharge','yes',
                 '-keepHs','False','-bcharge','0','-rundir','/home/jp/Runs/\n',
                 '-oxstate','II','-spin','1','-qccode TeraChem'])
#subprocess.call(["python","molSimplify/main.py","-core","Fe","-coord",'4',"-lig","cyclam, water","-ligocc",'1',
#                 '-geometry','sqp','-distort','0','checkdirb','True','-ligalign','False','-calccharge','yes',
#                 '-keepHs','False','-bcharge','0','-rundir','/home/jp/Runs/\n','-jobdir consti','-debug True',
#                 '-oxstate','II','-spin','5','-name',mol_name+'\n','-qccode TeraChem'])


print(datetime.now() - start)
