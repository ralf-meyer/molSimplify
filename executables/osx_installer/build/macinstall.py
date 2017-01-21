#!/usr/bin/env python
# Written by Tim Ioannidis for HJK Group
# Dpt of Chemical Engineering, MIT


import os, sys, subprocess, getpass
from distutils.spawn import find_executable

########################################
### module for running bash commands ###
########################################
def mybash(cmd):
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    stdout = []
    while True:
        line = p.stdout.readline()
        stdout.append(line)        
        if line == '' and p.poll() != None:
            break
    return ''.join(stdout)

# check command exists
def checkexist(cmd):
    exist = True
    try:
        subprocess.call(cmd)
    except OSError as e:
        exist = False
    return exist


if __name__ == '__main__':
    # check if brew is installed
    brew = checkexist('brew')
    # get current user
    user = mybash("logname").replace('\n','')
    # install brew
    instbrew = '/usr/bin/ruby -e '
    instbrew += '"$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"'
    instbrew += ' </dev/null'
    if not brew:
        # install brew as normal user
        cmd = 'sudo -u ' + user +' '+ instbrew
        t = mybash(cmd)
        t = mybash('chown root /usr/local/bin/brew') # change permissions
    # check if openbabel is installed
    try:
        import openbabel
    except:
        t = mybash('/usr/local/bin/brew install open-babel --with-python')
        t = mybash('mkdir -p /Users/'+user+'/Library/Python/2.7/lib/python/site-packages')
        cmd = "echo 'import site; site.addsitedir("+'"/usr/local/lib/python2.7/site-packages")'+"'"
        cmd += '>> /Users/'+user+'/Library/Python/2.7/lib/python/site-packages/homebrew.pth'
        t = mybash(cmd)
    # check if imagemagick is installed
    immag = checkexist('convert')
    if not immag:
        t = mybash('/usr/local/bin/brew install imagemagick')
        t = mybash('/usr/local/bin/brew install ghostscript')
        t = mybash('/usr/local/bin/brew link --overwrite ghostscript')




