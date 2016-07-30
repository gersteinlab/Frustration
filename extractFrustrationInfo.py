"""extractFrustrationInfo

Usage:
  extractFrustrationInfo1.py -I <mappedSNPInfo> -nd <nativePDBDir> -md <mutPDBDir> -F <frstnExecDir> -P <pdbSeqDir> -O <frustrationOutFile>
  extractFrustrationInfo1.py (-h | --help)

Options:
  -h --help     Show this screen.

"""
from docopt import docopt

import sys
import os
import subprocess
import glob
from collections import OrderedDict
 
def extractChainInfo(input):

    """ Read the pdb file to extract the chainId info and map it to each residue belonging to that particular chain"""

    chainId = []
    chainIdnr = []
    keyList = []
    valueList = []

    chainInfoMap = {} 

    pdbFile = input
    
    try:

       fileName = open(pdbFile,"r")
       if fileName:
          for line in fileName:
              if line.split()[0] == 'ATOM':
                 chainId.append(line.split()[4][0])

       chainIdnr = list(OrderedDict.fromkeys(chainId))
 
       for index in range(0,len(chainIdnr)):
           keyList.append(index+1)
           valueList.append(chainIdnr[index])
   
       chainInfoMap = dict(zip(keyList,valueList))

    except IOError: 
           print 'There is no file named', pdbFile

    print chainInfoMap

    return chainInfoMap


def extractChainLength(inputSeq):

    """ extract the sequence of the PDB file and store it in a FASTA format file. """

    sequence = ''
    pdbSeq = arguments["<pdbSeqDir>"]+'pdb_seq.py' 


    pdbSeq = pdbSeq + ' '+inputSeq[:-4]+ ' > ' + inputSeq[:-4]+'.fasta'
    p1 = subprocess.Popen(pdbSeq,shell=True)
    p1.communicate()   

    fileInp = open(inputSeq[:-4]+'.fasta','r')
    if fileInp:
       for line in fileInp:
           line = line.strip()

           if line[0:1] != '>' and line[0:1] !=' ':
              sequence = sequence + line.strip()

    removeFASTA = 'rm '+inputSeq[:-4]+'.fasta'
    p2 = subprocess.Popen(removeFASTA,shell=True)
    p2.communicate()   

    return  len(sequence)

 


def extractResNumFrustration(inputPDB,resIndex,flag):

    """ extract the residule level frustration index for each mutated residue position """

    RMFI = "NA"
    BMFI = "NA"
    AA = "NA"

    ### frustration output file name ###
    furstrationFileName = './'+flag+'/'+inputPDB+'.pdb_res_renum'

    try:
       fileName = open(furstrationFileName,"r")
       if fileName:
          for index in range(1,2):
              fileName.next()

          for line in fileName:
              line = line.strip()

              if not line:
                 RMFI = "NA"
                 BMFI = "NA"
 
              #''' extract relevant info if the residue number matches''' 

              elif int(''.join(filter(lambda x: x.isdigit(), line.split()[0]))) == int(resIndex):
                 RMFI = line.split()[2] #residue level rustration index 
                 BMFI = line.split()[1] #burial frustration index
                 AA = line.split()[3]   # Amino acid identity
    except IOError: 
           print 'There is no file named', furstrationFileName
   
    return RMFI,BMFI,AA




def calcFrustration(inputPDB,chainId,inputFlag,arguments):

    """ run residue level frustration calculation for the given PDB and store the output in relevant directory """

    chainLength = 0

    residueFrustration = 0.0
    burialFrustration = 0.0

    configFrustration = 0.0

    AccessSurfaceArea = ()


    pdbSubset = arguments["<pdbSeqDir>"]+'pdb_subset.py' 
    FrustMeter = arguments["<frstnExecDir>"]

    nativePDBDir = arguments["<nativePDBDir>"]
    mutatedPDBDir = arguments["<mutPDBDir>"]


    ''' for the native PDB frustration output gets stored in the native subdirectory'''

    if inputFlag == 'native':
       
 
       nativeInput = inputPDB.split('/')[2][:-4]+'.pdb_'+chainId+'.pdb' #nativePdb file 

       ''' run pdbTool to extract coOrdinates for the given pdbchain and store in pdb format '''

       pdbFile = inputPDB[0:8]

       if not os.path.exists('./native/'):
         mkNative = 'mkdir  ./native'
         p00 = subprocess.Popen(mkNative,shell=True)
         p00.communicate()

       if pdbFile :
          pdbNativeSubunit = 'python '+ pdbSubset + ' '+inputPDB+ ' -c '+ chainId
          print pdbNativeSubunit
          p0 = subprocess.Popen(pdbNativeSubunit,shell=True)
          p0.communicate() 

          ''' calculate chain length '''

          chainLength = extractChainLength(nativePDBDir+nativeInput)

          if chainLength <= 10000:

             #''' Move the native PDB file from the original native directory '''

             moveNative = 'cp '+nativePDBDir+nativeInput+' .'
             print moveNative
             p0 = subprocess.Popen(moveNative,shell=True)
             p0.communicate() 

             ''' Run frustration code on the native PDB'''
             nativeFrustratometer = 'sh '+FrustMeter + ' ' +nativeInput
             print nativeFrustratometer
             p1 = subprocess.Popen(nativeFrustratometer,shell=True)
             p1.communicate() 

             ''' Move the residue level frutsration output for the given PDB to the native directory'''  

             nativeFindexOut = 'mv ./'+nativeInput+'.done/'+nativeInput+'_res_renum   ./native/'
             print nativeFindexOut
             p2 = subprocess.Popen(nativeFindexOut,shell=True)
             p2.communicate()


             ''' Move the PDB file to the native subdirectory directory'''  

             nativePdbOut  = 'mv '+nativeInput+' ./native/'
             print nativePdbOut
             p4 = subprocess.Popen(nativePdbOut,shell=True)
             p4.communicate()


             ''' Delete the frustration code generated ouput directory to save space and get rid of other files'''  

             removeNatFrustration = 'rm -r ' +nativeInput+'*'
             print removeNatFrustration
             p5 = subprocess.Popen(removeNatFrustration,shell=True)
             p5.communicate()
          

    #''' for the mutated PDB, frustration output gets stored in the mutated subdirectory'''

    elif inputFlag == 'mutated':

         mutPDB = inputPDB+'_'+chainId+'.pdb'
         mutInput = inputPDB.split('/')[-1:][0]+'_'+chainId+'.pdb' 

         if not os.path.exists('./mutated/'):
            mkMutated = 'mkdir  ./mutated'
            p00 = subprocess.Popen(mkMutated,shell=True)
            p00.communicate()

         if mutPDB:

           #''' calculate chain length for mutated PDB '''
           chainLength = extractChainLength(mutatedPDBDir+mutInput)

           if chainLength <= 10000:
     

              pdbSubunit =  'python '+ pdbSubset + ' '+inputPDB+ ' -c '+ chainId
              print pdbSubunit
              p01 = subprocess.Popen(pdbSubunit,shell=True)
              p01.communicate() 

              moveMut = 'cp '+mutPDB+' .'
              print moveMut
              p0 = subprocess.Popen(moveMut,shell=True)
              p0.communicate() 


              ''' Run frustration code on the mutated PDB'''

              runFrustratometer = FrustMeter + ' ' +mutInput
              print runFrustratometer 
              p02 = subprocess.Popen(runFrustratometer,shell=True)
              p02.communicate() 


              ''' Move the residue level frutsration output for the mutated PDB to the mutated directory'''  

              mutFindexOut = 'mv ./'+mutInput+'.done/'+mutInput+'_res_renum   ./mutated/'
              print mutFindexOut
              p03 = subprocess.Popen(mutFindexOut,shell=True)
              p03.communicate()



              ''' Move the mutated PDB file to the mutated subdirectory directory'''  

              mutPdbOut  = 'cp '+mutInput+' ./mutated/'
              print mutPdbOut
              p05 = subprocess.Popen(mutPdbOut,shell=True)
              p05.communicate()

              ''' Delete the frustration code generated ouput directory to save space and get rid of other files'''  

              removeMutFrustration = 'rm -r ' +mutInput+'*'
              print removeMutFrustration
              p06 = subprocess.Popen(removeMutFrustration,shell=True)
              p06.communicate()


 
def extractFIndexVal(arguments):

    ''' extract out the residue level frustration index along with other info '''


    residueMAP = {'NA':'UNW1','X':'UNW','A':'ALA','C':'CYS','D':'ASP','R':'ARG','N':'ASN','E':'GLU','Q':'GLN','G':'GLY','H':'HIS','I':'ILE','L':'LEU','K':'LYS','M':'MET','F':'PHE','P':'PRO','S':'SER','T':'THR','W':'TRP','Y':'TYR','V':'VAL'}


    nativeAA = 'NA'
    mutAA = 'NA'
    MAF = 0.0

    inputInfoFile = arguments["<mappedSNPInfo>"]
    frstnOutFile = arguments["<frustrationOutFile>"]
    nativePDBDir = arguments["<nativePDBDir>"]
    mutatedPDBDir = arguments["<mutPDBDir>"]

    ### open output file in the append mode ###

    outputFile = open(frstnOutFile,'a')


    ### open the ENST file to read mutated and Native PDB along with snpInfo ###
    fileName = open(inputInfoFile,"r")
    if fileName:
       for line in fileName:
           line = line.strip()

           inputFileMutated = mutatedPDBDir+'.'.join(line.split()[2].split('.')[0:2])+line.split()[0][4:]+'.pdb'  #mutated pdb file name

           inputFileNative = nativePDBDir+line.split()[2][0:4]+'.pdb' #native pdb file name

           mutResIndex = line.split()[2].split('.')[1][6:] #residue number of the mutated residue
           print inputFileMutated 

           chainIdMap = extractChainInfo(inputFileMutated)

           if chainIdMap:

              nativeInput = inputFileNative.split('/')[-1:][0]+'_'+chainIdMap[int(inputFileMutated.split('.')[2])] #frustation output filename for native
              mutInput = inputFileMutated.split('/')[-1:][0]+'_'+chainIdMap[int(inputFileMutated.split('.')[2])]  # frustration output file name after mutation

              ### Extract the residue level frustration index, burial frustration and native amoni acid identity###
              nativeResFrustration,nativeBuriedFrustration,nativeAA = extractResNumFrustration(nativeInput,mutResIndex,'native')
  
              ### Extract the residue level frustration index, burial frustration and mutated amino acid identity###
              mutResFrustration,mutBuriedFrustration,mutAA = extractResNumFrustration(mutInput,mutResIndex,'mutated')
              
 
              outputFile.write(line.split()[0]+'\t'+line.split()[2]+'\t'+line.split()[1]+'\t'+residueMAP[nativeAA]+'\t'+residueMAP[mutAA]+'\t'+str(nativeResFrustration)+'\t'+str(mutResFrustration)+'\t'+str(nativeBuriedFrustration)+'\t'+str(mutBuriedFrustration))
              outputFile.write('\n')

          
def returnPDBMap(mappedSNVInp,arguments):

    ''' return pdb map with nativePdbId,mutatedPdbId along with the chainInfo'''

    chainIdMap = {}
    nativePdbMap = {}
    mutatedPdbMap = {}

    nativePDB = []
    mutatedPDB = []

    chainIdList = []
 
    nativePDBDir = arguments["<nativePDBDir>"]
    mutatedPDBDir = arguments["<mutPDBDir>"]


    fileName = open(mappedSNVInp,"r")
    if fileName:
       for line in fileName:
           line = line.strip()

           ''' list of mutated pdbIds'''
           inputFileMutated = mutatedPDBDir+'.'.join(line.split()[2].split('.')[0:2])+line.split()[0][4:]+'.pdb' 
           mutatedPDB.append(inputFileMutated)

           
           ''' list of native pdbIds'''
           inputFileNative = nativePDBDir+line.split()[2][0:4]+'.pdb'
           nativePDB.append(inputFileNative)

           chainIdMap = extractChainInfo(inputFileMutated)

           if chainIdMap: 

               ''' list of chainIds mapping to different pdbs'''
               chainInfo = chainIdMap[int(inputFileMutated.split('.')[2])]
               chainIdList.append(chainInfo)

    
    ''' map native PdbIds with corresponding chain'''
    nativePdbMap = dict(zip(nativePDB,chainIdList))

    ''' map mutated PdbIds with corresponding chain'''
    mutatedPdbMap = dict(zip(mutatedPDB,chainIdList))

    return nativePdbMap,mutatedPdbMap


def main(arguments):

    nativePdb = {}
    mutatedPdb = {}


    mappedSNPInfo = arguments["<mappedSNPInfo>"]
    nativePDBDir = arguments["<nativePDBDir>"]
    mutatedPDBDir = arguments["<mutPDBDir>"]

    print mappedSNPInfo+'\t'+nativePDBDir+'\t'+mutatedPDBDir


    ### obtain native & mutated pdbId along with chain Info #### 
    nativePdb,mutatedPdb = returnPDBMap(mappedSNPInfo,arguments)

    #### For each native Pdb run frustratometer and store the frustation files in native sub directory#### 
    for key in  nativePdb.keys():
        calcFrustration(key,nativePdb[key],'native',arguments)

    #### For each mutated Pdb run frustratometer and store the frustation files in mutated sub directory#### 
    for key in  mutatedPdb.keys():
        calcFrustration(key,mutatedPdb[key],'mutated',arguments)

    #### Extract out frustartion value along with other info ####
    extractFIndexVal(arguments)


           
if __name__ == '__main__':

   arguments = docopt(__doc__, version='extractFrustrationInfo 2.0')
   main(arguments)
