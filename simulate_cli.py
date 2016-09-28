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
subprocess.call(["python",'main.py','-slab_gen','-cif_path','Unitcells/anatase_tio2.cif',
                 '-slab_size','[9,9,9]','-miller_index','[1,0,1]','-place_on_slab','-target_molecule','/home/jp/Runs/co.xyz',
                 '-align_method','alignpair','-surface_atom_type','Ti','-num_surface_atoms','1',
                 '-object_align','C','-align_dist','2.5'])
#subprocess.call(["python",'main.py','-slab_gen','-cif_path','Unitcells/anatase_tio2.cif',
#                 '-slab_size','[9,9,9]','-miller_index','[1,0,1]'])



#
#subprocess.call(["python","main.py","-core","Fe","-coord",'6',"-lig","acac,water","-ligocc",'2,2',
#                 '-geometry','oct','-distort','0','checkdirb','True','-ligalign','False',
#3                 '-keepHs','False','-bcharge','0','-rundir','/home/jp/Runs/CLI',
#                 '-oxstate','0','-spin','1'])
print(datetime.now() - start)
