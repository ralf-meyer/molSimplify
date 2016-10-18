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

#print(strc)
#subprocess.call(strc,shell=True)
subprocess.call(["python","-m",'molSimplify.__main__'," -core"+" Co "+"-coord "+' 6 '+"-lig " +" acac,acac,acac "+
                 '-rundir ' + '/home/jp/Dropbox/Main/testingANN\n ' + ' -jobdir t2 ' + '-name coacacac_NN '+
                 '-geometry'+' oct '+'-distort'+' 0 ''checkdirb'+' True '+'-ligalign '+' True '+
                 '-keepHs'+' yes,yes,yes '+'-bcharge'+' 0 '+'-calccharge yes ' + '-qccode TeraChem '+ '-method UDF ' +
                 '-oxstate'+' III '+'-spin'+' 1 ' + ' -runtyp gradient ' ])

subprocess.call(["python","-m",'molSimplify.__main__'," -core"+" Co "+"-coord "+' 6 '+"-lig " +" acac,acac,acac "+
                 '-rundir ' + '/home/jp/Dropbox/Main/testingANN\n ' + ' -jobdir t2 ' + '-name coacacac_MS '+ ' -skipANN ' +
                 '-geometry'+' oct '+'-distort'+' 0 ''checkdirb'+' True '+'-ligalign '+' True '+
                 '-keepHs'+' yes,yes,yes '+'-bcharge'+' 0 '+ '-calccharge yes ' + '-qccode TeraChem '+ '-method UDF ' + 
                 '-oxstate'+' III '+'-spin'+' 1' + ' -runtyp gradient '])

subprocess.call(["python","-m",'molSimplify.__main__'," -core"+" Co "+"-coord "+' 6 '+"-lig " +" en,en,en "+
                 '-rundir ' + '/home/jp/Dropbox/Main/testingANN\n ' + ' -jobdir t2 ' + '-name coenenen_NN '+
                 '-geometry'+' oct '+'-distort'+' 0 ''checkdirb'+' True '+'-ligalign '+' True '+
                 '-keepHs'+' yes,yes,yes '+'-bcharge'+' 0 '+ '-calccharge yes ' + '-qccode TeraChem '+ '-method UDF ' + 
                 '-oxstate'+' III '+'-spin'+' 1' + ' -runtyp gradient '])


subprocess.call(["python","-m",'molSimplify.__main__'," -core"+" Co "+"-coord "+' 6 '+"-lig " +" en,en,en "+
                 '-rundir ' + '/home/jp/Dropbox/Main/testingANN\n ' + ' -jobdir t2 ' + '-name coenenen_MS '+ ' -skipANN ' +
                 '-geometry'+' oct '+'-distort'+' 0 ''checkdirb'+' True '+'-ligalign '+' True '+
                 '-keepHs'+' yes,yes,yes '+'-bcharge'+' 0 ' + '-calccharge yes ' + '-qccode TeraChem '+ '-method UDF ' + 
                 '-oxstate'+' III '+'-spin'+' 1' + ' -runtyp gradient '])


subprocess.call(["python","-m",'molSimplify.__main__'," -core"+" Fe "+"-coord "+' 6 '+"-lig " +" acac,acac,acac "+
                 '-rundir ' + '/home/jp/Dropbox/Main/testingANN\n ' + ' -jobdir t2 ' + '-name feacacac_NN '+
                 '-geometry'+' oct '+'-distort'+' 0 ''checkdirb'+' True '+'-ligalign '+' True '+
                 '-keepHs'+' yes,yes,yes '+'-bcharge'+' 0 '+ '-calccharge yes ' + '-qccode TeraChem '+ '-method UDF ' +
                 '-oxstate'+' III '+'-spin'+' 6' + ' -runtyp gradient '])


subprocess.call(["python","-m",'molSimplify.__main__'," -core"+" Fe "+"-coord "+' 6 '+"-lig " +" acac,acac,acac "+
                 '-rundir ' + '/home/jp/Dropbox/Main/testingANN\n ' + ' -jobdir t2 ' + '-name feacacac_MS '+ ' -skipANN ' +
                 '-geometry'+' oct '+'-distort'+' 0 ''checkdirb'+' True '+'-ligalign '+' True '+
                 '-keepHs'+' yes,yes,yes '+'-bcharge'+' 0 '+ '-calccharge yes ' + '-qccode TeraChem '+ '-method UDF ' + 
                 '-oxstate'+' III '+'-spin'+' 6' + ' -runtyp gradient '])



subprocess.call(["python","-m",'molSimplify.__main__'," -core"+" Cr "+"-coord "+' 6 '+"-lig " +" bipy,bipy,bipy "+
                 '-rundir ' + '/home/jp/Dropbox/Main/testingANN\n ' + ' -jobdir t2 ' + '-name crbipbipbip_NN '+
                 '-geometry'+' oct '+'-distort'+' 0 ''checkdirb'+' True '+'-ligalign '+' True '+
                 '-keepHs'+' yes,yes,yes '+'-bcharge'+' 0 '+ '-calccharge yes ' + '-qccode TeraChem '+ '-method UDF ' +
                 '-oxstate'+' II '+'-spin'+' 3' + ' -runtyp gradient '])


subprocess.call(["python","-m",'molSimplify.__main__'," -core"+" Cr "+"-coord "+' 6 '+"-lig " +" bipy,bipy,bipy "+
                 '-rundir ' + '/home/jp/Dropbox/Main/testingANN\n ' + ' -jobdir t2 ' + '-name crbipbipbip_MS '+
                 '-geometry'+' oct '+'-distort'+' 0 ''checkdirb'+' True '+'-ligalign '+' True '+ '-skipANN '
                 '-keepHs'+' yes,yes,yes '+'-bcharge'+' 0 '+ '-calccharge yes ' + '-qccode TeraChem '+ '-method UDF ' +
                 '-oxstate'+' II '+'-spin'+' 3' + ' -runtyp gradient '])





subprocess.call(["python","-m",'molSimplify.__main__'," -core"+" Mn "+"-coord "+' 6 '+"-lig " +" misc,misc,misc,misc,misc,misc"+
                 '-rundir ' + '/home/jp/Dropbox/Main/testingANN\n ' + ' -jobdir t2 ' + '-name mnmisc_NN '+
                 '-geometry'+' oct '+'-distort'+' 0 ''checkdirb'+' True '+'-ligalign '+' True '+
                 '-keepHs'+' yes,yes,yes '+'-bcharge'+' 0 '+ '-calccharge yes ' + '-qccode TeraChem '+ '-method UDF ' +
                 '-oxstate'+' II '+'-spin'+' 2 ' + ' -runtyp gradient '])



subprocess.call(["python","-m",'molSimplify.__main__'," -core"+" Mn "+"-coord "+' 6 '+"-lig " +" misc,misc,misc,misc,misc,misc"+
                 '-rundir ' + '/home/jp/Dropbox/Main/testingANN\n ' + ' -jobdir t2 ' + '-name mnmisc_MS '+
                 '-geometry'+' oct '+'-distort'+' 0 ''checkdirb'+' True '+'-ligalign '+' True '+ '-skipANN '
                 '-keepHs'+' yes,yes,yes '+'-bcharge'+' 0 '+ '-calccharge yes ' + '-qccode TeraChem '+ '-method UDF ' +
                 '-oxstate'+' II '+'-spin'+' 2 '+ ' -runtyp gradient '])





