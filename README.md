# Frustration Project

This repository contains source code for workflow evalauting change in Localized frustration index of protein residues
upon mutation. The input to this workflow is VAT (http://vat.gersteinlab.org/) output, which is ran on the user provided
list of single nucleotide variants(SNVs).
This workflow consist of three steps for evaluating changes in frustration:

1) Parsing VAT output of all SNVs to etxract residue position and residue identity for the mutated residue on protein sequence --- parseVatOut.py

2) Mapping each SNV onto user-provided list of PDB strcuture --- mapSNP2PDB.py

3) Evaluating Frustration changes of residues  --- extractFrustrationInfo.py
