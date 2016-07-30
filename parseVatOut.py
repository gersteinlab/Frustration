"""parseVatOut

Usage:
  parseVatOut.py -d <dataResource> -v <vatOut> -b <bioMartFile> -type <snpType>  
  parseVatOut.py (-h | --help)

Options:
  -h --help     Show this screen.

"""
from docopt import docopt
import sys
import os
import collections

##################################################################################
# parseVatOut.py is a  python script written to extract the coding SNV annotation#  
# Written by Sushant Kumar, Yale University, 01/06/2016                         #
##################################################################################


def extractFastaSeq(input1,input2,input3):

    """ obtain the FASTA sequence for the given PDB """
    lineNum = 0

    protSequence = ''
   
    for index in range(input2+1,input3): 
        protSequence = protSequence+input1[index]

    return protSequence            



def extractFASTAInfo(input):

    """ extract the information for each entry in the biomart FASTA file """

    lineNum = 0
  
    keyList = []
    valList = []

    seqList = []
    seqKeyList = []
    sequence = []
    fastaTupleList = []

    fastaTuple = ()

    fastaInfoMap = {}
    fastaSeqMap = {}
    fastaSeqIndx = {}

    tupVal = () 
 

    fileName = open(input,"r")
    if fileName:
       for line in fileName:
           line = line.strip()

           lineNum = lineNum + 1

           if line[0:1] == '>':
              keyList.append(lineNum)

              tupVal = (line.split('|')[1],line.split('|')[0][1:])
              valList.append(tupVal)
 
           elif line[0:1] != '>':

              seqKeyList.append(lineNum)
              sequence.append(line)

    ''' GeneId and TranscritpId info from biomart FASTA file'''
    fastaInfoMap = dict(zip(keyList,valList))

    ''' Fasta Sequence of the translated protein from the biomart FASTA file '''
    fastaSeqIndx = dict(zip(seqKeyList,sequence))

    #### extract sequence from the gencode translated fasta file##############
    for index in range(0,len(keyList)-1):
        seqList.append(extractFastaSeq(fastaSeqIndx,keyList[index],keyList[index+1]))

    seqList.append(extractFastaSeq(fastaSeqIndx,keyList[len(keyList)-1],lineNum))
    
 
    fastaSeqMap = dict(zip(keyList,seqList)) 

    sorted_fastaInfoMap = collections.OrderedDict(sorted(fastaInfoMap.items()))

    ####### generate list of tuple, which holds geneId, transcriptId & protein sequence #######
    for key in sorted_fastaInfoMap.keys():
        fastaTuple=(fastaInfoMap[key][0],fastaInfoMap[key][1],fastaSeqMap[key])
        fastaTupleList.append(fastaTuple)

    return fastaTupleList


def extractGeneID(inputVatFile,bioMartInfo,dataSource,snpType):

    """ Read in the vatoutput file and etxract the details of each non-synonymous SNP along with the original and mutated Sequence """    


    index = 0
    sequence = ''

    vIndex = []
    residueIndex = []

    fileOut = open(inputVatFile+'.'+snpType+".summary","a")

    fileName = open(inputVatFile,"r")
    if fileName:

       ##### skip headers in the vatOutput file #######

       #for indx in range(1,2):
       #    fileName.next()

       for line in fileName:
           line = line.strip()

           if not line.startswith("#"):

              ''' List of index for a given type of variant in the vat ouput file''' 

              if snpType in line.split()[7].split(':'):
                 vIndex = [i for i, x in enumerate(line.split()[7].split(':')) if x == snpType]

               
              for element in vIndex:
    
               
                  offSet = int(line.split()[7].split(':')[element+1].split('/')[0])   ### number of transcript for the given gene listed in the vatOutput file
                  transcriptIndex = element + 3 ## index for the first transcript listed in the vatOutput file
       
 
                  ''' Check whether ensembel geneId for the variant is present in the gencode19 annotation file'''

                  if line.split()[7].split(':')[element - 2] in [item[0] for item in bioMartInfo]:

                  
                     ''' list out indices for each transcript corresponding to the given geneId in the vatoutput which is also present in gencode19''' 

                     residueIndex = [i for i, x in enumerate([item[0] for item in bioMartInfo]) if x == line.split()[7].split(':')[element - 2]]

               
                     ''' For each element in the trancript list for the given geneId reported by gencode19, check if that particular transcript is present in the vat output file'''

                     for member in residueIndex:
                  
                         while index < 3*offSet-2:
                                         
                               ''' If the transcript is listed in the vat ouput file, get relevant information like residue index for mutation, mutated and original residue and sequence of the translated form of the transcript'''

                               if line.split()[7].split(':')[transcriptIndex+index] ==  [item[1] for item in bioMartInfo][member]:

                                  resIndex = [item[1] for item in bioMartInfo].index(line.split()[7].split(':')[transcriptIndex+index])
                                  sequence = [item[2] for item in bioMartInfo ][resIndex]

                                  mutationIndex = int(line.split()[7].split(':')[transcriptIndex+index+1].split('_')[2])-1
                                  origResidue = line.split()[7].split(':')[transcriptIndex+index+1].split('_')[3].split('->')[0]
                                  mutResidue = line.split()[7].split(':')[transcriptIndex+index+1].split('_')[3].split('->')[1]

                                  try:    
                                     ''' If the original listed reported by VAT matches with original residue in the FASTA sequence , mutate that position and report relevant info'''  
                  
                                     if sequence[mutationIndex] == origResidue:
                                        mutSequence = sequence[0:mutationIndex]+mutResidue+sequence[mutationIndex+1:] # mutated protein sequence
                                        SNPINFO = line.split()[0]+':'+line.split()[1]+':'+line.split()[3]+':'+line.split()[4] # snpInfo = chromNum:snvPos:refAllele:altAllele

                                        snpAttribute = line.split()[7].split(':')[0].split(';')[2]+':'+':'.join(line.split()[7].split(':')[2].split(';')[1:]) # Additional SNP attributes: VAF,geneName etc. 

                                        #print line.split()[7].split(':')[element-2]+"\t"+SNPINFO+'\t'+snpAttribute+'\t'+line.split()[7].split(':')[element]+'\t'+line.split()[7].split(':')[transcriptIndex+index]+"\t"+line.split()[7].split(':')[transcriptIndex+index+1]


                                        fileOut.write(line.split()[7].split(':')[element-2]+"\t"+SNPINFO+'\t'+snpAttribute+'\t'+line.split()[7].split(':')[element]+"\t"+line.split()[7].split(':')[transcriptIndex+index]+"\t"+line.split()[7].split(':')[transcriptIndex+index+1]+"\t"+str(sequence)+"\t"+str(mutSequence))

                                        fileOut.write("\n")

                                  except Exception:
                                         print "exception found\n"

                               index = index + 3


                         index = 0 
     
              vIndex = []
  
def main(arguments):


    dataResource = arguments["<dataResource>"]
    vatOutputFile = arguments["<vatOut>"]
    biomartFile = arguments["<bioMartFile>"]
    codingSNPType = arguments["<snpType>"]

    ensembleOutList = extractFASTAInfo(dataResource)
        
    extractGeneID(vatOutputFile,ensembleOutList,biomartFile,codingSNPType)


if __name__ == '__main__':

    arguments = docopt(__doc__, version='parseVatOut 1.0')
    main(arguments)

