# -*- coding: utf-8 -*-
"""
PURPOSE: 
    This script will vary certain LOSALAMOS model parameters within some
    user-defined range to scope out the parameter space and determine the 
    best fit values

Created on Mon Nov 16 13:49:17 2015

@author: cns <Christopher N. Shingledecker>
"""

import math,random,os,subprocess

contFile = 'parameters.f03' #File to be edited
exFile = 'losalamos'        #Binary to run
Nruns = 500                 #Number of simulation runs

#Set the ranges for the parameters to be varied
#(1) Trial frequency
trlnumax = 1E10
trlnumin = 1E11

#(2) Dissociation probability
#NB: Since DISPROB is already between 0,1, there is no need to give max and min

#(3) Number of sub-excitation interactions
nsubexmax = 50
nsubexmin = 0

#(4) Most probable gamma-distribution value
avalmax = 30
avalmin = 3


os.system('rm -rf sim_no*') #Clean up old directories
os.system('rm -rf temp')    #Clean up the temp directory

random.seed(716381)

for j in range(0,Nruns):
    #Mkdir with os package
    temp_path = 'temp'
    os.system('rm -rf temp')
    os.mkdir(temp_path)
    
    #Copy files into new directory
    cmd = 'cp * ./' + temp_path + '/'
    os.system(cmd)

    #Get random numbers
    r1 = random.random()
    r2 = random.random()
    r3 = random.random()
    r4 = random.random()
    r5 = random.random()
    r6 = random.random()
        
    #Calculate new parameter values
    #(1) Trial frequency
    trlnu = trlnumin + r1*(trlnumax - trlnumin)

    #(2) Dissociation probability
    disprob = r2 

    #(3) Sub-excitation number
    nsubex = nsubexmin + r3*(nsubexmax-nsubexmin)

    #(4) A-values
    aval = avalmin + r4*(avalmax-avalmin)

    #(5) Determine fast-reacts
    #NB: fastreact and fragile are strings already
    if r5 <= 0.5:
        fastreact = "(/ 4  /)"
    else:
        fastreact = "(/ 21 /)"

    #(6) Determine fragile list
    if r6 <= 0.5:
        fragile = "(/ 7  /)"
    else:
        fragile = "(/ 21 /)"
    

    #Edit the input files
    infile = './' + temp_path + '/' + contFile
#    print infile
    i = 0
    inf =  open(infile,'r')
    lines = inf.readlines()
    inf.close()
#    outfile = './' + temp_path + '/new_' + contFile
    outf = open(infile,'w')
    for line in lines:
        i = i + 1
        if i == 30:
            newline = line[:14] + '%1.3e' % temp_zeta + line[24:]
            newline = newline[:19] + 'D' + newline[20:]
            outf.write(newline)
        elif i == 73:
            newline = line[:14] + '%.6e' % oH2 + line[26:]
            newline = newline[:22] + 'D' + newline[23:]
            outf.write(newline)
        elif i == 74:
            newline = line[:14] + '%.6e' % pH2 + line[26:]
            newline = newline[:22] + 'D' + newline[23:]
            outf.write(newline)
        else:
            outf.write(line)
    outf.close()
    
    #Change into new directory
    os.chdir(temp_path)
    
    #Recompile with os
    print 'Now recompiling'
    cmd = './' + compFile
    subprocess.call(cmd,shell=True)

    #Run with os
    print 'Starting Nautilus'
    cmd = './'  + exFile
    subprocess.call(cmd,shell=True)

    #Change back into home directory
    os.chdir('..')
    newdir = 'both_run_'+str(j)
    os.mkdir(newdir)
    cmd = 'cp ./temp/output_1D* ./' + newdir + '/'
    os.system(cmd)
    cmd =  'cp ./temp/rates1D.* ./' + newdir + '/'  
    os.system(cmd)
    cmd = 'rm -rf temp/'
    os.system(cmd)

print 'Ending simulation runs'
