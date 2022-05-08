# Oncomerge Pathways

This program takes somatic, CNA, and fusion inputs from GISTIC and MutSig2CV andoutputs a combined matrix that  integrates the different mutation types. It does the same for pathway mutations, classifying them into Activation, Loss of Function, and pure old Protein Affecting Mutations. To quantify the significance of this classifcation, a permutation test is perform by creating fake pathways from random genes and comparing the mutation frequency of those fake pathway to real ones. The pvalues are then correcting with a Benjamini-Hochberg procedure.

The Oncomerge.py file only contains code for performing oncomerge on individual genes. OncomergePathways.py contains code for both individual genes and pathways.

The OncomergePathways code was run on all 33 cancer types in the the TCGA dataset. The folder named output contains the outputof this run when including kegg, GO and human cyc pathways. The foler named output alternative pathways runs the code with PID, OncoSig, and cancer hallmark pathways.

Have Fun!
