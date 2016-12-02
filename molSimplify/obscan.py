import os
import subprocess
import shutil
from datetime import datetime
start = datetime.now()



this_distance = 1.5
while (this_distance < 3.0):
    this_name  ='dist_' + str(int(1000*this_distance))
    target_geopath = '/home/jp/test/geos/' + this_name + '.xyz'
    target_jobpath = '/home/jp/test/jobs/' + this_name + '.in' 
    base_jobpath = '/home/jp/test/input.in' 
    subprocess.call(["python","-m",'molSimplify.__main__','-place_on_slab','-unit_cell','/home/jp/test/fakeslab.xyz','-cell_vector','[[5,0,0],[0,5,0],[0,0,5]]',
                '-target_molecule','/home/jp/test/o2.xyz','align_distance_method','custom','-align_dist',str(this_distance),'-align_method','alignpair','-object_align','1',
                '-surface_atom_type','Fe','-debug','True','-rundir','/home/jp/temptest'])
    shutil.move('/home/jp/temptest/loaded_slab/loaded.xyz',target_geopath)
    this_distance += 0.025
    with open(target_jobpath,'w') as new:
        with open(base_jobpath,'r') as old:
            for line in old:
                new.write(line)
            new.write('coordinates '+ target_geopath + '\n')
            new.write('scrdir scr/'+ this_name + '\n')

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
