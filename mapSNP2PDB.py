"""mapSNP

Usage:
  mapSNP2PDB.py -p <pdbIdList> -b <bioMartFile> -I <snpSummaryFile> -B <blastPDir> -M <modellerDir> -P <pbdSeqDir> -O <outLogFile>
  mapSNP2PDB.py (-h | --help)

Options:
  -h --help     Show this screen.

"""
from docopt import docopt

import sys
import os
import time
import subprocess
import read
import operator
import timeit

def filterPDBs(pdbFileList):

    """ extract PDB Ids from the filtered PDB list (filtering based on resolution and R-factor of crystal structures only """
    
    filteredIds = []

    fileName = open(pdbFileList,'r')
    if fileName:
       for line in fileName:
           line = line.strip()
           filteredIds.append(line.split()[0])

    return filteredIds 


def extractPdbId(bioMartInput,pdbL):

    """ read the biomart file and  extract GeneId, transcriptId & PdbId for all filtered PDB dataset """

    pdbTups = ()
    pdbIdList = []

    fileName = open(bioMartInput,"r")
    if fileName:
       for i in range(1,2):
           fileName.next()
       for line in fileName:
           line = line.strip()

           if line.split(',')[3] in pdbL:
              pdbTups = (line.split(',')[0],line.split(',')[1],line.split(',')[3])
              pdbIdList.append(pdbTups)

    return set(pdbIdList)     



def generatePdbChain(inputPDBFile):

    """ extract the pdbChain Info """

    lineNum = 0
    lineNumList = []

    pdbId = inputPDBFile[:-6]

    fileInp = open(inputPDBFile,'r')
    if fileInp:
       for line in fileInp:
           line = line.strip()
           lineNum = lineNum + 1
 
           if line[0:1] == '>':              
              lineNumList.append(lineNum)

   
    lineNumList.append(lineNum)      
    return lineNumList  
    

def runBlastP(pdbId,subjectStr,transcriptId,arguments):

    """ run blastP on the transcript fasta sequence against PDB sequence. Map the SNP on the PDB file."""

    chainTups = ()

    chainTupsList = []
    chainInfoList = []

    chainStartPos = []
    bestAlignIndex = []

    identityMap = {}

    chainInfoMap = {}

    chainSequence = ''

    residueMAP = {'A':'ALA','C':'CYS','D':'ASP','R':'ARG','N':'ASN','E':'GLU','Q':'GLN','G':'GLY','H':'HIS','I':'ILE','L':'LEU','K':'LYS','M':'MET','F':'PHE','P':'PRO','S':'SER','T':'THR','W':'TRP','Y':'TYR','V':'VAL'} 
    
    index = 1
    lineNo = 0
    fileInpName = ''
    sequence = ''


    blastP = arguments["<blastPDir>"]+'blastp'
    pdbSeq = arguments["<pbdSeqDir>"]+'pdb_seq.py'

    pdbFile = pdbId .lower()

    """ extract the sequence of the PDB file and store it in a FASTA format file. """

    pdbSeqRun = pdbSeq+ ' ' +pdbFile+ ' > ' + pdbFile[:-4]+'.fasta' 
    p0 = subprocess.Popen(pdbSeqRun,shell=True)
    p0.communicate()


    """ extract the line number corresponding to the start of sequence of each chain in the FASTA file """

    chainStartPos = generatePdbChain(pdbFile[:-4]+'.fasta')



    """ write the sequence of native protein in a FASTA format. """

    fileSubject = open('subject'+pdbFile[:-4]+'.fasta','w')
    fileSubject.write(">sequence \n")
    fileSubject.write(subjectStr.strip('*')) 
    fileSubject.close()



    """ For each chain of the PDB file generate a separate FASTA file. """ 

    while index < len(chainStartPos):

               fileInp = open(pdbFile[:-4]+'.fasta','r')
               if fileInp:
                  for line in fileInp:
                      line = line.strip()
                      lineNo = lineNo + 1

                      if lineNo >= int(chainStartPos[index-1]) and lineNo < int(chainStartPos[index]) and line.strip():

                         chainSequence = chainSequence + line.strip()

                         fileInpName = str(pdbFile[:-4])+"."+str(index)+".fasta"
                         fileOut = open(fileInpName,"a")
                         fileOut.write(line.strip())
                         fileOut.write("\n")

               fileOut.close() 

 
               if  len(chainSequence) <= 10000: 

 
                  """ run blastP on each PDB hcain sequence against the  native protein sequence of the transcript. """  
    
                  blastSeq = blastP +' -query '+pdbFile[:-4]+"."+str(index)+'.fasta' + "  -subject subject"+pdbFile[:-4]+".fasta >  " +pdbFile[:-4]+'.'+str(index)+'.'+transcriptId+'.fasta' 
                  p = subprocess.Popen(blastSeq,shell=True)
                  p.communicate()

                  #### extract Alignment index , chainInfo for the chain with best alignment score  with respect to native sequence present in biomart file #### 
                  bestAlignIndex,identityMap,chainInfo = read.extractBestAlign(pdbFile[:-4]+'.'+str(index)+'.'+transcriptId+'.fasta')


                  if len(identityMap) == 1 and 0.0 in identityMap.values():
                     return

                  else:

                      ### Store the linumber of the alignment file, alignment identity score, last line of alignment and chainIndex for the best alignment in a tuple
                      chainTups = ( max(identityMap.iteritems(), key=operator.itemgetter(1))[0] , max(identityMap.iteritems(), key=operator.itemgetter(1))[1], bestAlignIndex[bestAlignIndex.index(max(identityMap.iteritems(), key=operator.itemgetter(1))[0])+1], index )
                      
                      ### compile the above tuple info in a list
                      chainTupsList.append(chainTups) 

                      ### STore the PdbId and chainId  
                      chainInfoList.append(chainInfo)

               ### If the chain length is greater than the prescribed cutoff don't bother about it !!!
               else:
                   return
           
               lineNo = 0
               chainSequence = ''

               index = index + 1


    ### Map alignment description of each chain with its chain info(pdbId & chainId) ###
    chainInfoMap = dict(zip(chainInfoList,chainTupsList))

    return chainInfoMap



def generatePDBModel(pdbInfoMap,pdbId,snpInfo,transcriptId,snpId,arguments):

    residueIndexMap = {}
    transcriptIndexMap = {}

    modelerSH = arguments["<modellerDir>"]+'modpy.sh'
    mutateModel = arguments["<modellerDir>"]+'mutate_model.py'
    outFile = arguments["<outLogFile>"]
 
    residueMAP = {'A':'ALA','C':'CYS','D':'ASP','R':'ARG','N':'ASN','E':'GLU','Q':'GLN','G':'GLY','H':'HIS','I':'ILE','L':'LEU','K':'LYS','M':'MET','F':'PHE','P':'PRO','S':'SER','T':'THR','W':'TRP','Y':'TYR','V':'VAL'} 
   
    currDir = os.getcwd()

    if not os.path.exists(currDir+'/nativePDB/'):
       os.mkdir('nativePDB')

    if not os.path.exists(currDir+'/mutatedPDB/'):
       os.mkdir('mutatedPDB')

    nativePDBDir = currDir+'/nativePDB/'
    mutPDBDir = currDir+'/mutatedPDB/'

    fileLog = open(outFile,'a') #store the geneId, transcriptId and pdbId for each SNP record in the log file

    pdbFile = pdbId.lower()

    if pdbInfoMap:


       #### extract the chainInfo along with identity score for the best aligned chain of a given PDB
       pdbKey = max(pdbInfoMap.items(), key=lambda(k,v):v[1])[0]
       chainIndex = pdbInfoMap[max(pdbInfoMap.items(), key=lambda(k,v):v[1])[0]][3]

       residueIndexMap,transcriptIndexMap = read.extractAlignInfo(pdbFile[:-4]+'.'+str(chainIndex)+'.'+transcriptId+'.fasta',pdbInfoMap[pdbKey],pdbKey)

               
       if residueIndexMap and transcriptIndexMap:

          if int(snpInfo.split("_")[2]) in transcriptIndexMap.keys() and snpInfo.split("_")[2] in residueIndexMap.keys():


             if residueIndexMap[snpInfo.split("_")[2]] == snpInfo.split("_")[3].split("->")[0]:

                #### Run Modeler to obtain mutated PDB file ####
                mutateResidue = modelerSH + ' python '+mutateModel+ ' '+pdbFile+' '+snpInfo.split("_")[2]+' '+residueMAP[snpInfo.split("_")[3].split("->")[1][0:1]]+ ' '+pdbKey[4:5]
                p1 = subprocess.Popen(mutateResidue,shell=True)
                p1.communicate()

                     
                #### copy the original  PDB file ####
                moveNatPDB = 'cp '+pdbFile[:-4].lower()+'.pdb'+' '+ nativePDBDir
                print moveNatPDB
                p2 = subprocess.Popen(moveNatPDB,shell=True)
                p2.communicate()


                #### copy the mutated  PDB file ####
                moveMutPDB = 'mv '+pdbFile[:-4].lower()+'.pdb'+residueMAP[snpInfo.split("_")[3].split("->")[1][0:1]]+snpInfo.split("_")[2]+'.pdb'+ ' '+ mutPDBDir+pdbFile[:-4].lower()+'.pdb'+residueMAP[snpInfo.split("_")[3].split("->")[1][0:1]]+snpInfo.split("_")[2]+'.'+str(chainIndex)+'.'+transcriptId+'.pdb'
                print moveMutPDB
                p3 = subprocess.Popen(moveMutPDB,shell=True)
                p3.communicate()

                ### Write the description output file ###
                fileLog.write(pdbFile[:-4].lower()+'.'+str(chainIndex)+'.'+transcriptId+"\t"+snpId+'\t'+pdbFile[:-4].lower()+'.pdb'+residueMAP[snpInfo.split("_")[3].split("->")[1][0:1]]+snpInfo.split("_")[2]+'.pdb'+snpInfo.split("_")[2]+'.pdb')
                fileLog.write("\n")

    ### remove the fasta file ###
    removeFile = 'rm *.fasta'
    p4 = subprocess.Popen(removeFile,shell=True)
    p4.communicate()
               


def mapSNP2PDB(snpSummary,biomartPDBList,arguments):

    """ read the SNP summary file and map each filtered PdbId with the corresponding SNP record """ 

    pdbChainMap = {}
    
    fileName = open(snpSummary,"r")
    if fileName:
       for line in fileName:
           line = line.strip()


           for item in biomartPDBList:

               ### If the transcriptId and geneId for the given transcript matches with biomart entry obtain corresponding PDBId
               if line.split()[0][0:15]  == item[0] and line.split()[4][0:15] == item[1]:

                  pdbId = item[2]+'.pdb' ## Matching transcript's pdbId
                  transcriptInfo = line.split()[0]+'_'+line.split()[4]  ## concatanate geneId & transcriptId for the given transcript
                 
                  snvInfo = line.split()[5] ## SNP info
                  sequence = line.split()[6] ## AA sequence
                  snvId = line.split()[1]    ## SNP Id


                  pdbChainMap = runBlastP(pdbId,sequence,transcriptInfo,arguments) 

                  generatePDBModel(pdbChainMap,pdbId,snvInfo,transcriptInfo,snvId,arguments)
                 

def main(arguments):

    #start = timeit.timeit()

    filterPDB = []
    pdbIdInfo = set()

    pdbList = arguments["<pdbIdList>"]
    bioMartFile = arguments["<bioMartFile>"]
    snpSummaryFile = arguments["<snpSummaryFile>"] 

    blastP = arguments["<blastPDir>"]+'blastp'
    modelerSH = arguments["<modellerDir>"]+'modpy.sh'
    mutateModel = arguments["<modellerDir>"]+'mutate_model.py'
    pdbSeq = arguments["<pbdSeqDir>"]+'pdb_seq.py'
    outFile = arguments["<outLogFile>"]


    filterPDB = filterPDBs(pdbList)

    ''' Retrun the list of tuple(geneId,transcriptId,PdbId) '''
    pdbIdInfo = extractPdbId(bioMartFile,filterPDB)

    ''' Map SNP onto PdbId for the given trancript'''
    mapSNP2PDB(snpSummaryFile,pdbIdInfo,arguments)

    #end = timeit.timeit()
    #print end-start

if __name__ == '__main__':

    arguments = docopt(__doc__, version='mapSNP2PDB 1.0')
    main(arguments)

