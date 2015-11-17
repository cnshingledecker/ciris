# -*- coding: utf-8 -*-
"""
PURPOSE: 
    This script will vary certain LOSALAMOS model parameters within some
    user-defined range to scope out the parameter space and determine the 
    best fit values

Created on Mon Nov 16 13:49:17 2015

@author: cns <Christopher N. Shingledecker>
"""

import sys,math,random,os,subprocess

contFile = 'parameters.f03' #File to be edited
exFile = 'losalamos'        #Binary to run
Nruns = 2                 #Number of simulation runs
baserep = 'sim_no'

#Writing file format parameters
debug=False
lnstart=47
lnend=54
txtend=55

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


cmd = 'rm -rf ' + baserep + '*'
os.system(cmd) #Clean up old directories
os.system('rm -rf temp')    #Clean up the temp directory
os.system('rm *csv')

random.seed(716381)

for j in range(0,Nruns):
    #Mkdir with os package
    temp_path = 'temp'
    os.system('rm -rf temp')
    os.system('rm losalamos*')
    os.mkdir(temp_path)
    
    #Copy files into new directory
    cmd = 'cp * ./' + temp_path + '/'
    os.system(cmd)

    #Change into new directory
    os.chdir(temp_path)
    os.system('echo $PWD')


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
    nsubex = int(nsubexmin + math.floor((nsubexmax+1-nsubexmin)*r3))

    #(4) A-values
    aval = avalmin + r4*(avalmax-avalmin)

    #(5) Determine fast-reacts
    #NB: fastreact and fragile are strings already
    if r5 <= 0.5:
        fastreacts = "(/ 4  /)"
    else:
        fastreacts = "(/ 21 /)"

    #(6) Determine fragile list
    if r6 <= 0.5:
        fragile = "(/ 7  /)"
    else:
        fragile = "(/ 21 /)"

    #Edit the input files
    infile = contFile
    i      = 0
    inf    = open(infile,'r')
    lines  = inf.readlines()
    inf.close()
    outf   = open(infile,'w')
    for line in lines:
        i = i + 1
        # Change trlnu
        if i == 52:
            newline = line[:lnstart] + '%1.1e' % trlnu + line[lnend:]
            if debug==True:
                print trlnu
                print line
                print newline
            else:
                outf.write(newline)
        elif i == 53:
            newline = line[:lnstart] + '%1.2f' % disprob + line[lnend-3:]
            if debug==True:
                print disprob
                print line
                print newline
            else:
                outf.write(newline)
        elif i == 72:
            newline = line[:lnstart] + '%2d' % nsubex + line[lnend-5:]
            if debug==True:
                print nsubex
                print line
                print newline
            else:
                outf.write(newline)
        elif i == 76:
            newline = line[:lnstart] + '%1.1e' % aval + line[lnend:]
            if debug==True:
                print aval
                print line
                print newline
            else:
                outf.write(newline)
        elif i == 89:
            newline = line[:lnstart] + fastreacts + line[txtend:]
            if debug==True:
                print fastreacts
                print line
                print newline
            else:
                outf.write(newline)
        elif i == 90:
            newline = line[:lnstart] + fragile + line[txtend:]
            if debug==True:
                print fragile
                print line
                print newline
            else:
                outf.write(newline)
        else:
            outf.write(line)
    outf.close()

    if debug==False: 
      #Recompile with os
      print 'Now recompiling'
      cmd = 'make'
      subprocess.call(cmd,shell=True)
      cmd = 'make clean'
      subprocess.call(cmd,shell=True)

      #Run with os
      print 'Starting losalamos'
      cmd = './'  + exFile
      subprocess.call(cmd,shell=True)

    #Change back into home directory
    os.chdir('..')
    newdir = baserep + str(j)
    
    os.mkdir(newdir)
    cmd = 'cp ./temp/abundance.csv ./' + newdir + '/'
    os.system(cmd)
    cmd =  'cp ./temp/parameters.f03 ./' + newdir + '/'  
    os.system(cmd)
    cmd = 'cp ./quickplot.p ./' + newdir + '/'
    os.system(cmd)
    cmd = 'rm -rf temp/'
    os.system(cmd)

print 'Ending simulation runs'
