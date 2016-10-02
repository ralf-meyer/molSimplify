import os
import subprocess
from datetime import datetime
start = datetime.now()
strc =( "molsimplify"+" -core"+" Mn "+"-coord "+' 6 '+"-lig" +" en,en,co,co "+"-ligocc"+' 1,1,1,1,1,1 '+
                 '-rundir ' + '"/home/jp/Dropbox/Main/testingANN\n "' + ' -jobdir "t1" ' + '-name feaci '+
                 '-geometry'+' oct '+'-distort'+' 0 ''checkdirb'+' True '+'-ligalign '+' True '+
                 '-keepHs'+' False '+'-bcharge'+' 0 '+'-ligloc 1 '+
                 '-oxstate'+' 2 '+'-spin'+' 1')
strc =( "molsimplify"+" -core"+" Ni "+"-coord "+' 6 '+"-lig" +" water,water,water,water,water,water "+"-ligocc"+' 1,1,1,1,1,1 '+
                 '-rundir ' + '"/home/jp/Dropbox/Main/testingANN\n "' + ' -jobdir "t1" ' + '-name niwataci '+
                 '-geometry'+' oct '+'-distort'+' 0 ''checkdirb'+' True '+'-ligalign '+' True '+
                 '-keepHs'+' False '+'-bcharge'+' 0 '+'-ligloc 1 '+
                 '-oxstate'+' 2 '+'-spin'+' 1')

print(strc)
subprocess.call(strc,shell=True)
