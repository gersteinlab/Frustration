# Frustration Project

This repository contains source code for workflow evaluating the change in Localized frustration index of a protein residue
upon mutation. The input to this workflow is VAT output file, which is ran on the user provided
list of single nucleotide variants(SNVs). More details about running VAT can be found here (http://vat.gersteinlab.org/).

This workflow consist of three steps for evaluating changes in frustration:

1) Parsing VAT output of all SNVs to etxract residue position and residue identity for the mutated residue on protein sequence --- parseVatOut.py
Usage:
  parseVatOut.py -d <dataResource> -v <vatOut> -b <bioMartFile> -type <snpType>
  parseVatOut.py (-h | --help)

2) Mapping each SNV onto user-provided list of PDB strcuture --- mapSNP2PDB.py
Usage:
  mapSNP2PDB.py -p <pdbIdList> -b <bioMartFile> -I <snpSummaryFile> -B <blastPDir> -M <modellerDir> -P <pbdSeqDir> -O <outLogFile>
  mapSNP2PDB.py (-h | --help)

3) Evaluating Frustration changes of residues  --- extractFrustrationInfo.py
Usage:
  extractFrustrationInfo1.py -I <mappedSNPInfo> -nd <nativePDBDir> -md <mutPDBDir> -F <frstnExecDir> -P <pdbSeqDir> -O <frustrationOutFile>
  extractFrustrationInfo1.py (-h | --help)
