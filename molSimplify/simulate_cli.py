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
subprocess.call(["python",'-m','molSimplify.__main__','-slab_gen','-cif_path','molSimplify/Unitcells/anatase_tio2.cif',
                 '-slab_size','[9,9,9]','-target_molecule','/home/jp/Runs/co.xyz','-miller_index','[1,0,1]',
                 '-align_method','alignpair','-surface_atom_ind','39','-num_surface_atoms','1','-debug','True',
                 '-object_align','C','-align_dist','0.5','-shave_extra_layers','1'])
#
#subprocess.call(["python",'-m','molSimplify.__main__','-slab_gen','-cif_path','/home/jp/Dropbox/Main/icepaper/quartz/quartz.cif',
#                 '-slab_size','[9,9,6]','-rundir /home/jp/Dropbox/Main/icepaper/\n','-debug','-align_method','alignpair', '-freeze', '2',
#                 '-place_on_slab','-target_molecule','/home/jp/Dropbox/Main/icepaper/quartz/mol.xyz','-object_align', 'Mn','-expose_type','O',
#                 '-align_dist','1.5','-shave_extra_layers','0','align_distance_method','custom','-surface_atom_type','O','-num_surface_atoms','2'])
#
#subprocess.call(["python","-m",'molSimplify.__main__','-slab_gen','-unit cif_path','molSimplify/Unitcells/anatase_tio2.cif'
#                 '-slab_size','[9,9,9]','-miller_index','[1,0,0]'])

#subprocess.call(["python","-m",'molSimplify.__main__','-place_on_slab','-unit_cell','/home/jp/Downloads/cell.xyz','-cell_vector','[[7.5,0,0],[0,7.5,0],[0,0,19]]',
#                '-target_molecule','/home/jp/Runs/mo5.xyz','align_distance_method','custom','-align_dist','1.5','-align_method','alignpair','-object_align','O',
#                '-surface_atom_type','Ti','-debug','True'])


#subprocess.call(["python","-m",'molSimplify.__main__','slab_gen','-cif_path','/home/jp/Downloads/pd.cif','-slab_size','[10,10,10]','-miller_index','[1,1,1]',
#                '-freeze','2', '-surface_atom_type','Pd','-debug','True'])

#subprocess.call(["python","-m",'molSimplify.__main__','-place_on_slab','-slab_gen','-cif_path','/home/jp/Downloads/pd.cif','-slab_size','[10,10,7]',
#                '-target_molecule','/home/jp/Runs/mo5.xyz','align_distance_method','custom','-align_dist','1.5','-align_method','alignpair','-object_align','O','-freeze','2',
#                '-surface_atom_type','Pd','-debug','True'])


#
#subprocess.call(["python","main.py","-core","Fe","-coord",'6',"-lig","acac,water","-ligocc",'2,2',
#                 '-geometry','oct','-distort','0','checkdirb','True','-ligalign','False',
#3                 '-keepHs','False','-bcharge','0','-rundir','/home/jp/Runs/CLI',
#                 '-oxstate','0','-spin','1'])
print(datetime.now() - start)
