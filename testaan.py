import os
import subprocess
from datetime import datetime
start = datetime.now()

subprocess.call(["python","main.py","-core","Fe","-coord",'6',"-lig","acac,acac,water,water","-ligocc",'1,1,1,1',
                 '-rundir', '/home/jp/Dropbox/Main/testingANN\n','-jobdir t1', '-name feaci \n',
                 '-geometry','oct','-distort','0','checkdirb','True','-ligalign','True',
                 '-keepHs','False','-bcharge','0','-ligloc 1',
                 '-oxstate','0','-spin','1'])
