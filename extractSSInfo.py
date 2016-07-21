import sys
import os
import glob
import re
import subprocess

def readDSSPOut(input,input2):

    residueMAP = {'X':'UNW','A':'ALA','C':'CYS','D':'ASP','R':'ARG','N':'ASN','E':'GLU','Q':'GLN','G':'GLY','H':'HIS','I':'ILE','L':'LEU','K':'LYS','M':'MET','F':'PHE','P':'PRO','S':'SER','T':'THR','W':'TRP','Y':'TYR','V':'VAL'}
   
    secondaryStructure = 'NA'
 
    fileName = open(input,'r')
    if fileName:
       for index in range(1,29):
           fileName.next()

       for line in fileName:
           line = line.strip()
          
           if len(line[:-119].split()) == 5 and line[:-119].split()[3].upper() in residueMAP:

              chainInfo = line[:-119].split()[2]

              print line[:-119]

              snpInfo = residueMAP[line[:-119].split()[3].upper()]+line[:-119].split()[1]

    
              if residueMAP[input2[0:1]]+input2[1:] == snpInfo:


                 if line[:-119].split()[4] == 'H' or line[:-119].split()[4] == 'G' or line[:-119].split()[4] == 'I':
                    secondaryStructure = 'helix'

                 elif line[:-119].split()[4] == 'E' or line[:-119].split()[4] == 'B' or line[:-119].split()[4] == 'S':
                      secondaryStructure = 'strand'

                 elif line[:-119].split()[4] == 'T':
                      secondaryStructure = 'turn'

           elif len(line[:-119].split()) == 4 and line[:-119].split()[3].upper() in residueMAP:

                snpInfo = residueMAP[line[:-119].split()[3].upper()]+line[:-119].split()[1]

                if input2 == snpInfo :
                   secondaryStructure = 'coil'


    return secondaryStructure
        

def calSSInfo(input1,input2):

    cmd = '/home2/sk972/resNetwork/dsspcmbi'

    resInfo = input2

    dsspOut = cmd + " " + input1 + " > " + input1+".dssp" 
    p1 = subprocess.Popen(dsspOut,shell=True)
    p1.communicate()

    SS = readDSSPOut(input1+".dssp",resInfo)


    cmd1 = 'rm '+input1+'.dssp'
    p2 = subprocess.Popen(cmd1,shell=True)
    p2.communicate()

    return SS


