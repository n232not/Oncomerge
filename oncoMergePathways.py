# %%
##########################################################
## OncoMerge:  oncoMerge.py                             ##
##  ______     ______     __  __                        ##
## /\  __ \   /\  ___\   /\ \/\ \                       ##
## \ \  __ \  \ \___  \  \ \ \_\ \                      ##
##  \ \_\ \_\  \/\_____\  \ \_____\                     ##
##   \/_/\/_/   \/_____/   \/_____/                     ##
## @Developed by: Plaisier Lab                          ##
##   (https://plaisierlab.engineering.asu.edu/)         ##
##   Arizona State University                           ##
##   242 ISTB1, 550 E Orange St                         ##
##   Tempe, AZ  85281                                   ##
## @github: https://github.com/plaisier-lab/OncoMerge   ##
## @Author:  Chris Plaisier                             ##
## @License:  GNU GPLv3                                 ##
##                                                      ##
## If this program is used in your analysis please      ##
## mention who built it. Thanks. :-)                    ##
##########################################################


#####################
## Import packages ##
#####################

import json
import argparse
import pandas as pd
import numpy as np
from pyparsing import null_debug_action
from statsmodels.stats.multitest import multipletests # pip install statsmodels==0.10.0
import itertools
import sys
import os
from tqdm import tqdm
from multiprocessing import Pool
from statsmodels.stats.multitest import fdrcorrection
# %%

#import mygene
#mg = mygene.MyGeneInfo() # pull in function to map genes

##################################
## Read in command line options ##
##################################
###############################################################################################################

parser = argparse.ArgumentParser(description='OncoMerge merges patient Protein Affecting Mutations (PAMs) and Copy Number Alterations (CNAs) into a unified mutation matrix.')
parser.add_argument('-cf', '--config_file', help='Path to JSON encoded configuration file, overrides command line parameters', type = str)
parser.add_argument('-gp', '--gistic_path', help='Path to GISTIC output folder', type = str, default='C:/Users/nihan/Dropbox/OncoMerge_pathways/TCGA/GISTIC/BRCA')
parser.add_argument('-df', '--del_file', help='Path to GISTIC deletion file (default = del_genes.conf_99.txt)', type = str, default = 'del_genes.conf_99.txt')
parser.add_argument('-af', '--amp_file', help='Path to the GISTIC amplification file (default = amp_genes.conf_99.txt)', type = str, default = 'amp_genes.conf_99.txt')
parser.add_argument('-gdf', '--gene_data_file', help='Path to the GISTIC gene data file (default = all_data_by_genes.txt)', type = str, default = 'all_data_by_genes.txt')
parser.add_argument('-aaf', '--alternate_annotation_file', help='Supply alternate annotation file to convert gene symbols to Entrez IDs (default does not import alternate annotation file and instead uses conversion embedded in GISTIC output files).', type = str, default='C:/Users/nihan/Dropbox/OncoMerge_pathways/TCGA/OncoMerge_input_g2e_converter.csv')
parser.add_argument('-ln', '--label_name', help='Label for Entrez ID column in GISTIC gene data file (default = \'Gene ID\')', type = str, default = 'Locus ID')
parser.add_argument('-tf', '--thresh_file', help='Path to the GISTIC all_thresholded file (default = all_thresholded.by_genes.txt)', type = str, default = 'all_thresholded.by_genes.txt')
parser.add_argument('-pam', '--pam_file', help='Path to the protein affecting mutation (PAM) file (CSV matrix where columns are patients and genes are rows) [0 = not mutated, and 1 = mutated]', type = str, default='C:/Users/nihan/Dropbox/OncoMerge_pathways/TCGA/PAM/BRCA_somMutMC3.csv')
parser.add_argument('-mscv', '--mutsig2cv_file', help='Path to a MutSig2CV output file', type = str, default='C:/Users/nihan/Dropbox/OncoMerge_pathways/TCGA/MutSig2cv/BRCA_sig2cv.csv')
parser.add_argument('-op', '--output_path', help='Path you would like to output OncoMerged files (default = current directory)', type = str, default = '.')
parser.add_argument('-fus', '--fusions_file', help='Path to the gene fusions file (CSV matrix where columns are patients and genes are rows) [0 = not fused, and 1 = fused]', type = str, default='C:/Users/nihan/Dropbox/OncoMerge_pathways/TCGA/FUSIONS/BRCA_fusions.csv')
parser.add_argument('-mmf', '--min_mut_freq', help='Minimum frequency of mutation (range = 0-1; default = 0.05)', type = float, default = 0.05)
parser.add_argument('-pq', '--perm_qv', help='Permuted p-value FDR BH corrected cutoff (default = 0.1)', type = float, default = 0.1)
parser.add_argument('-sp', '--save_permutation', help='Run and save out permutation analysis to be used for comparability in another OncoMerge run (default off)', action='store_true')
parser.add_argument('-lp', '--load_permutation', help='Do not run permutation anlaysis and load permutation anlaysis from previous run (default off)', type = str, default = None)
parser.add_argument('-mlg', '--min_loci_genes', help='Minimum number of genes in loci to apply maximum final frequency filter (default = 10)', type = int, default = 10)
parser.add_argument('-mpf', '--min_pam_freq', help='Minimum PAM frequency (default = 0.01)', type = float, default = 0.01)
parser.add_argument('-tcga', '--tcga', help='Clip gistic TCGA names.', type = bool, default = True)
parser.add_argument('-bl', '--blacklist', help='List of patients (one per line) to exclude for frequency calculations.', type = str, default='C:/Users/nihan/Dropbox/OncoMerge_pathways/TCGA/blacklist/blacklist_29850653_29625053.csv')
parser.add_argument('-g2g', '--gene2go', help='Gene2Go pathways', type = str, default = 'pathways/gene2go.hsa.csv')
parser.add_argument('-kegg', '--kegg', help='Kegg pathways', type = str, default = 'pathways/keggPathwayGenes_hsa_3_10_2019.csv')
parser.add_argument('-huam', '--huamncyc', help='huamncyc pathways', type = str, default='pathways/huamncyc_PC2_3_11_2019.csv')
parser.add_argument('-pid', '--pid', help='PID pathways', type = str, default = 'pathways/PID.csv')
parser.add_argument('-os', '--oncosig', help='Oncogenic Signature pathways', type = str, default = 'pathways/OncoSig.csv')
parser.add_argument('-hall', '--hallmark', help='Hallmark pathways', type = str, default='pathways/hall.csv')
args = parser.parse_args()

def pathwayPool(i):
    perms = []
    tmp1 = [np.array(i[1].loc[np.random.choice(i[1].index, i[0])].sum(axis=0)) for j in range(i[4])]
    tmp2 = [np.array(i[3].loc[np.random.choice(i[3].index, i[0])].sum(axis=0)) for j in range(i[4])]
    if type(i[2])!=str:
        tmp3 = [np.array(i[2].loc[np.random.choice(i[2].index, i[0])].sum(axis=0)) for j in range(i[4])]
        temp = []
        for j in range(i[4]): 
            temp.append((tmp1[j] + tmp2[j] + tmp3[j]).clip(0,1).mean())
        perms = [str(i[0]), temp]
    else:
        temp = []
        for j in range(i[4]):
            temp.append((tmp1[j] + tmp2[j]).clip(0,1).mean())
        perms = [str(i[0]), temp]
    return perms
#%%
if __name__ == "__main__":



    #######################
    ## Define parameters ##
    #######################

    params = args.__dict__
    #params['fusions_file'] = 'C:/Users/nihan/Dropbox/OncoMerge_pathways/TCGA/FUSIONS/BRCA_fusions.csv'
    if args.config_file:
        with open(args.config_file, "r") as cfg:
            tmp = json.loads(cfg.read())
            for i in tmp:
                params[i] = tmp[i]
    if (not params['gistic_path']) or (not params['pam_file']) or (not params['mutsig2cv_file']) or (params['save_permutation'] and not params['load_permutation']==None):
        parser.print_help()
        sys.exit(1)


    ##################
    ## Load up data ##
    ##################

    print('Loading data...')

    # Create conversion series for gene symbol to Entrez ID
    if not params['alternate_annotation_file']:
        n1 = pd.read_csv(params['gistic_path']+'/'+params['gene_data_file'],index_col=0,sep='\t',usecols=[0,1])
        n1.index = [i.split('|')[0] for i in n1.index]
        n1 = n1.drop_duplicates()
        n1 = n1[params['label_name']]
    else:
        n1 = pd.read_csv(params['alternate_annotation_file'],index_col=0)[params['label_name']].apply(int)

    # load up significantly mutated genes
    mutSig2CV = pd.read_csv(params['mutsig2cv_file'],index_col=1)
    mutSig2CV = mutSig2CV.loc[mutSig2CV.index.map(lambda x: x in n1.index)]
    mutSig2CV.index = mutSig2CV.index.map(lambda x: n1.loc[x])
    mutSig2CV = mutSig2CV.loc[~mutSig2CV.index.duplicated(keep='first')]
    sigPAMs = list(mutSig2CV.index[mutSig2CV['q']<=0.1]) # Set at 0.1 based on Lawerence et al., Nature 2013 (PMID = 23770567)

    # Get list of significantly CNA amplified genes
    ampLoci = {}
    ampLoci_qv = {}
    amp1 = pd.read_csv(params['gistic_path']+'/'+params['amp_file'],index_col=0,sep='\t')
    for col1 in amp1.columns:
        if float(amp1[col1]['residual q value'])<=0.05 and not (col1[0]=='X' or col1=='Y'):
            ampLoci[col1] = list(set([n1.loc[i] for i in [i.lstrip('[').rstrip(']').split('|')[0] for i in list(amp1[col1].dropna()[3:])] if i in n1.index and n1.loc[i]>0]))
            ampLoci_qv[col1] = float(amp1[col1]['residual q value'])

    # Get list of significantly CNA deleted genes
    delLoci = {}
    delLoci_qv = {}
    del1 = pd.read_csv(params['gistic_path']+'/'+params['del_file'],index_col=0,sep='\t')
    for col1 in del1.columns:
        if float(del1[col1]['residual q value'])<=0.05 and not (col1[0]=='X' or col1=='Y'):
            delLoci[col1] = list(set([n1.loc[i] for i in [i.lstrip('[').rstrip(']').split('|')[0] for i in list(del1[col1].dropna()[3:])] if i in n1.index and n1.loc[i]>0]))
            delLoci_qv[col1] = float(del1[col1]['residual q value'])

    # Set up background gene numbers for gold-standard enrichment analysis
    backgrounds = {'Activating':[], 'LossOfFunction':[], 'Aggregate':[]}
    actTmp = [gene for locus in ampLoci for gene in ampLoci[locus]] #question on what this creates
    backgrounds['Activating'] += actTmp
    backgrounds['Aggregate'] += actTmp
    delTmp = [gene for locus in delLoci for gene in delLoci[locus]]
    backgrounds['LossOfFunction'] += delTmp
    backgrounds['Aggregate'] += delTmp
    backgrounds['Aggregate'] += sigPAMs

    # Load up somatically mutated genes
    somMuts = pd.read_csv(params['pam_file'],index_col=0,header=0)
    somMuts[somMuts != 0] = 1
    if not somMuts.index.dtype=='int64':
        somMuts = somMuts.loc[somMuts.index.map(lambda x: x in n1.index)]
        somMuts.index = somMuts.index.map(lambda x: n1.loc[x])

    somMuts = somMuts.loc[~somMuts.index.duplicated(keep='first')]


    # Load up gene fusions
    
    if not params['fusions_file']=='none':
        fusions = pd.read_csv(params['fusions_file'],index_col=0,header=0)
        fusions[fusions != 0] = 1
        if not fusions.index.dtype=='int64':
            fusions = fusions.loc[somMuts.index.map(lambda x: x in n1.index)]
            fusions.index = fusions.index.map(lambda x: n1.loc[x])
        fusions = fusions.loc[~fusions.index.duplicated(keep='first')]

    # Read in gistic2 all_data_by_genes file
    with open(params['gistic_path']+'/'+params['thresh_file'],'r') as inFile:
        tmp = inFile.readline().strip().split('\t')
        numCols1 = len(tmp)

    d1 = pd.read_csv(params['gistic_path']+'/'+params['thresh_file'],index_col=0,sep='\t').drop(tmp[1], axis = 1)
    if params['tcga']:
        d1.columns = [i[:12] for i in d1.columns]

    d1 = d1.loc[d1.index.map(lambda x: x.split('|')[0] in n1.index)]
    d1.index = d1.index.map(lambda x: n1.loc[x.split('|')[0]])
    d1.index.name = 'Locus ID'

    # Removing sex chromosomes (issues in CNA analysis) from d1
    lociThresh = d1['Cytoband']
    include = []
    for i in lociThresh:
        if not(i[0]=='X' or i[0]=='Y'):
            include.append(True)
        else:
            include.append(False)

    d1 = d1.loc[lociThresh[include].index].drop('Cytoband', axis = 1)

    # Make sure somMuts and gistic have same samples
    if params['fusions_file']=='none':
        pats = list(set(d1.columns).intersection(somMuts.columns))
    else:
        pats = list((set(d1.columns).intersection(somMuts.columns)).intersection(fusions.columns))
        fusions = fusions[pats]

    somMuts = somMuts[pats]
    d1 = d1[pats]

    # Get rid of duplicated rows
    d1 = d1[~d1.index.duplicated(keep='first')]

    ## Fill out summary matrix
    cols = ['Symbol', 'PAM_freq', 'MutSig2CV_qvalue']
    if not params['fusions_file'] == 'none':
        cols += ['Fusion_freq']
    cols += ['CNA_freq', 'CNA_locus', 'CNA_type', 'GISTIC_residual_q_value', 'Act_freq', 'LoF_freq', 'OM_type_selected', 'OM_empirical_p_value', 'OM_empirical_q_value', 'Genes_in_locus',  'Final_mutation_type', 'Final_freq', 'Delta_over_PAM']
    summaryMatrix = pd.DataFrame(index= list(set([gene for locus in ampLoci.values() for gene in locus] + [gene for locus in delLoci.values() for gene in locus] + list(somMuts.index))), columns = cols)

    # Add gene symbols
    toMatch = [i for i in summaryMatrix.index if i in n1.values]
    summaryMatrix.loc[toMatch,'Symbol'] = [n1.index[n1==i][0] for i in toMatch]

    # Add MutSig2CV q-values
    toMatch = [i for i in summaryMatrix.index if i in mutSig2CV.index]
    summaryMatrix.loc[toMatch,'MutSig2CV_qvalue'] = mutSig2CV.loc[toMatch,'q']

    # Load up black list of patients if given one
    blacklist = []
    if not params['blacklist']=='none':
        with open(params['blacklist'], 'r') as inFile:
            while 1:
                line = inFile.readline()
                if not line:
                    break
                blacklist.append(line.strip().split(',')[0])

    #Print out some useful information
    print('\tSize of somatic mutation matrix: '+str(somMuts.shape))
    if not params['fusions_file']=='none':
        print('\tSize of fusions matrix: '+str(fusions.shape))
    print('\tSize of CNA matrix: '+str(d1.shape))
    print('Finished loading data.')

    # Make the output directory if it doesn't exists already
    if not os.path.exists(params['output_path']):
        os.makedirs(params['output_path'])


    ########################
    ## Begin OncoMerging ##
    ########################


    # Precompute positive and negative dichotomized matrices
    print('Precomputing dichotomized matrices...')
    posdicot = (lambda x: 1 if x>=2 else 0)
    posD1 = d1.applymap(posdicot)
    posFreq = posD1[posD1.columns.difference(blacklist)].mean(axis=1)
    ampGenes = {int(j):i for i in ampLoci for j in ampLoci[i]}
    negdicot = (lambda x: 1 if x<=(-2) else 0)
    negD1 = d1.applymap(negdicot)
    negFreq = negD1[negD1.columns.difference(blacklist)].mean(axis=1)
    delGenes = {int(j):i for i in delLoci for j in delLoci[i]}
    print('Finished precomputing dichotomized matrices.')

    # %%
    # Add back genes that are high frequency amplification or deletion (20%)
    for gene1 in list([i for i in posD1.index if posFreq[i]>=0.2]):
        if (not gene1 in ampGenes) and isinstance(lociThresh[gene1], str) and (lociThresh[gene1] in ampLoci_qv) and (gene1 in posD1.index):
            ampGenes[gene1] = lociThresh[gene1]
            ampLoci[lociThresh[gene1]].append(gene1)
            #print('Pos', gene1, posFreq[gene1], lociThresh[gene1])

    for gene1 in list([i for i in negD1.index if negFreq[i]>=0.2]):
        if (not gene1 in delGenes) and isinstance(lociThresh[gene1], str) and (lociThresh[gene1] in delLoci_qv) and (gene1 in negD1.index):
            delGenes[gene1] = lociThresh[gene1]
            delLoci[lociThresh[gene1]].append(gene1)
            #print('Neg', gene1, negFreq[gene1], lociThresh[gene1])

    # %%
    #Add Pathways
    #gene2go = pd.read_csv(params['gene2go'], header=None).iloc[: , 1:]
    #kegg = pd.read_csv(params['kegg'], header=0).iloc[: , 1:]
    #huam = pd.read_csv(params['huamncyc'], header=None, sep=",").iloc[: , 1:3]

    #gene2go.iloc[:,0] ='GO::' + gene2go.iloc[:,0]
    #kegg.iloc[:,0] ='kegg::' + kegg.iloc[:,0]
    #huam.iloc[:,0] ='humanCyc::' + huam.iloc[:,0]
    #kegg.columns = pd.Index([1,2])

    PID = pd.read_csv(params['pid'], header=0, index_col=0)
    OncoSig = pd.read_csv(params['oncosig'], header=0, index_col=0)
    hall = pd.read_csv(params['hallmark'], header=0, index_col=0)

    #pathways =pd.DataFrame(pd.concat([PID, OncoSig, hall], axis=0).iloc[:,0:2]).dropna()
    pathways = pd.concat([PID, OncoSig, hall], axis=0)
    pathways.columns = ['Pathway', 'Genes']
    #pathways['Genes'] = pathways['Genes'].apply(lambda x: list(map(int, str(x).split(' '))))
    pathways['Genes'] = pathways['Genes'].apply( lambda x: list(map(int, list(str(x[1:-1]).split(", ")))))
    pathways = pathways.loc[pathways['Genes'].map(len) > 1]
    genes = pathways['Genes'].map(set)
    pathways['Genes'] = genes.map(list)
    pathways = pathways[~pathways.duplicated(subset = "Pathway", keep='last')]
    pathways.index = range(len(pathways.index))
    pathways["PathwayID"] = pd.Series(pathways.index).astype(int) + 1000000001
    pathways["Pathway"] = pathways["Pathway"].str.replace('_',"-")

    # %%
    inSigPAMs=[]
    for i in pathways.index:
        inSigPAMs.append(any(x in sigPAMs for x in pathways["Genes"][i]))

    # %%
    inAmpGenes=[]
    ampList= [int(x) for x in list(ampGenes.keys())]
    for i in pathways.index:
        inAmpGenes.append(any(x in ampList for x in pathways["Genes"][i]))

    # %%
    inDelGenes=[]
    delList= [int(x) for x in list(delGenes.keys())]
    for i in pathways.index:
        inDelGenes.append(any(x in delList for x in pathways["Genes"][i]))

    # %%
    genes = pathways['Genes'].map(set)

    # %%
    m1 = somMuts.sum(axis=1).clip(0,1)
    m1set = set(m1[m1>0].index)
    mutgenes = [list(m1set & i) for i in genes]
    pathways["Somatically Mutated Genes"] = mutgenes  

    # %%
    pathways["Number of Somatically Mutated Genes"] = [len(i) for i in mutgenes]

    # %%
    m1 = posD1.sum(axis=1).clip(0,1)
    m1set = set(m1[m1>0].index)
    mutgenes = [list(m1set & i) for i in genes]
    pathways["Amplified Genes"] = mutgenes  

    # %%
    m1 = negD1.sum(axis=1).clip(0,1)
    m1set = set(m1[m1>0].index)
    mutgenes = [list(m1set & i) for i in genes]
    pathways["Deleted Genes"] = mutgenes  
    # %%
    pathways.to_csv("pathways.csv")
    # %%
    def pathwiseAddition(origDf): 
        indices = set(origDf.index)
        genes = pathways['Genes'].map(set)
        rows= [pd.DataFrame(origDf.loc[list((indices & i))].sum(axis=0)) for i in genes]
        return pd.concat(rows,axis=1).T.set_index(pathways['PathwayID']).clip(0,1)

    # Include pathway name when sum
    #Instead of subsettign w/ somMutPoints, filter based on freq later
    #dump out new mutation frequency for pathway and test w/ amp loci del loci, SomMutPoitn etc
    # Histogram of freqs, seperated by mut type and whether or not pathways qould be included in amp loci, del loci, etc

    somMutsPathways = pathwiseAddition(somMuts)
    posD1Pathways = pathwiseAddition(posD1)
    negD1Pathways = pathwiseAddition(negD1)
    fusionPathways= pathwiseAddition(fusions)

    # %%
    inSigPAMs = list(np.array(inSigPAMs) & np.array(somMutsPathways.mean(axis=1) > params['min_mut_freq']))
    inAmpGenes = list(np.array(inAmpGenes) & np.array(posD1Pathways.mean(axis=1) > params['min_mut_freq']))
    inDelGenes = list(np.array(inDelGenes) & np.array(negD1Pathways.mean(axis=1) > params['min_mut_freq']))
    # %%
    somMutsPathways.index = somMutsPathways.index.astype(int)
    posD1Pathways.index = posD1Pathways.index.astype(int)
    negD1Pathways.index = negD1Pathways.index.astype(int)
    fusionPathways.index = fusionPathways.index.astype(int)

    # %%
    ampLoci.update(dict(zip(list(pathways[inAmpGenes]["Pathway"]), list(posD1Pathways[inAmpGenes].index.map(lambda x: list([x]))))))
    ampGenes.update(dict(zip(list(posD1Pathways[inAmpGenes].index), list(pathways[inAmpGenes]["Pathway"]))))

    # %%
    delLoci.update(dict(zip(list(pathways[inDelGenes]["Pathway"]), list(negD1Pathways[inDelGenes].index.map(lambda x: list([x]))))))
    delGenes.update(dict(zip(list(negD1Pathways[inDelGenes].index), list(pathways[inDelGenes]["Pathway"]))))

    # %%
    sigPAMs = sigPAMs + list(somMutsPathways[inSigPAMs].index)

    # %%
    somMuts = pd.concat([somMuts, somMutsPathways], axis = 0)
    posD1 = pd.concat([posD1, posD1Pathways], axis = 0)
    negD1 = pd.concat([negD1, negD1Pathways], axis = 0)
    fusions = pd.concat([fusions, fusionPathways], axis = 0)

    # %%
    summaryMatrix1 = summaryMatrix.copy()

    # %%
    summaryMatrix = pd.DataFrame(index= list(set([gene for locus in ampLoci.values() for gene in locus] + [gene for locus in delLoci.values() for gene in locus] + list(somMuts.index) + list(pathways["PathwayID"]))), columns = cols)

    # %%
    summaryMatrix.loc[summaryMatrix1.index] = summaryMatrix1

    # %%
    del summaryMatrix1

    # %%
    summaryMatrix = summaryMatrix[~summaryMatrix.index.duplicated(keep='first')]

    # %%
    freq1 = somMuts[somMuts.columns.difference(blacklist)].sum(axis=1)/len(list(somMuts.columns.difference(blacklist)))
    somMutPoint = freq1[freq1>=params['min_mut_freq']].index
    tmp1 = [i for i in summaryMatrix.index if i in freq1.index]
    summaryMatrix.loc[tmp1, 'PAM_freq'] = freq1.loc[tmp1]

    # %%
    summaryMatrix["Pathway_Length"] = int(1)

    # %%
    somMutInd = set(somMuts.index)
    pathways["length"] = pathways["Genes"].apply(lambda x: len(set(x) & somMutInd))

    # %%
    summaryMatrix.loc[pathways['PathwayID'], "Pathway_Length"] = list(pathways["length"])

    # %%
    summaryMatrix["Pathway_Name"] = ""

    # %%
    summaryMatrix.loc[pathways['PathwayID'], "Pathway_Name"] = list(pathways["Pathway"])
    # %%
    summaryMatrix["Somatically_Mutated_Genes"] = float("NaN")
    summaryMatrix["Amplified_Genes"] = float("NaN")
    summaryMatrix["Deleted_Genes"] = float("NaN")

    # %%
    pathways['PathwayID']

    # %%
    summaryMatrix.loc[pathways['PathwayID'], "Somatically_Mutated_Genes"] = list(pathways["Somatically Mutated Genes"])
    summaryMatrix.loc[pathways['PathwayID'], "Amplified_Genes"] = list(pathways["Amplified Genes"])
    summaryMatrix.loc[pathways['PathwayID'], "Deleted_Genes"] = list(pathways["Deleted Genes"])

    # %%
    # Cutoff fusions based on the minimum mutation frequency (mf)
    if params['fusions_file']:
        freq2 = fusions[fusions.columns.difference(blacklist)].sum(axis=1)/len(list(fusions.columns.difference(blacklist)))
        fusionsMut = freq2[freq2>=params['min_mut_freq']].index
        tmp1 = [i for i in summaryMatrix.index if i in freq2.index]
        summaryMatrix.loc[tmp1, 'Fusion_freq'] = freq2.loc[tmp1]

    posFreq = posD1[posD1.columns.difference(blacklist)].mean(axis=1)
    negFreq = negD1[negD1.columns.difference(blacklist)].mean(axis=1)

    # %%
    n1Path = pd.Series(pathways['PathwayID'])

    # %%
    n1Path.index = pathways['Pathway']

    # %%
    n1 = pd.concat([n1,n1Path])

    # %%
    print('Generating CNA summaryMatrix data...')
    somMuts1 = list(set(list(ampGenes.keys())+list(delGenes.keys())).intersection(list(set(list(posD1.index)) | set(list(negD1.index)) | set(list(d1.index)))))
    delGenesSet = set(delGenes.keys())
    ampGenesSet = set(ampGenes.keys())
    bothGenesSet = (ampGenesSet.intersection(delGenesSet)).intersection(somMuts1)
    onlyDel = delGenesSet.difference(ampGenesSet).intersection(somMuts1)
    onlyAmp = ampGenesSet.difference(delGenesSet).intersection(somMuts1)
    # %%
    with tqdm(total=len(bothGenesSet)+len(onlyDel)+len(onlyAmp)) as pbar:
        for both1 in bothGenesSet:
            # Choose amplification if is higher frequency, otherwise choose deletion
            if (both1 < 1000000000):
                if posFreq[both1]>negFreq[both1]:
                    summaryMatrix.loc[both1,['CNA_type','CNA_locus','GISTIC_residual_q_value','CNA_freq']] = ['Amp', ampGenes[both1], ampLoci_qv[ampGenes[both1]], posFreq[both1]]
                # Otherwise choose deletion (deletion becomse default for equal frequency; unlikely to happen)
                else:
                    summaryMatrix.loc[both1,['CNA_type','CNA_locus','GISTIC_residual_q_value','CNA_freq']] = ['Del', delGenes[both1], delLoci_qv[delGenes[both1]], negFreq[both1]]
            else:
                if posFreq[both1]>negFreq[both1]:
                    summaryMatrix.loc[both1,['CNA_type','CNA_freq']] = ['Amp', posFreq[both1]]
                # Otherwise choose deletion (deletion becomse default for equal frequency; unlikely to happen)
                else:
                    summaryMatrix.loc[both1,['CNA_type','CNA_freq']] = ['Del', negFreq[both1]]
            pbar.update(1)
        # Straight up deletion
        for del1 in onlyDel:
            if (del1 < 1000000000):
                summaryMatrix.loc[del1,['CNA_type','CNA_locus','GISTIC_residual_q_value','CNA_freq']] = ['Del', delGenes[del1], delLoci_qv[delGenes[del1]], negFreq[del1]]
            else:
                summaryMatrix.loc[del1,['CNA_type','CNA_freq']] = ['Del', negFreq[del1]]
            pbar.update(1)
        # Straight up amplification
        for amp1 in onlyAmp:
            if (amp1 < 1000000000):
                summaryMatrix.loc[amp1,['CNA_type','CNA_locus','GISTIC_residual_q_value','CNA_freq']] = ['Amp', ampGenes[amp1], ampLoci_qv[ampGenes[amp1]], posFreq[amp1]]
            else: 
                summaryMatrix.loc[amp1,['CNA_type','CNA_freq']] = ['Amp', posFreq[amp1]]

            pbar.update(1)
# %%
    # Bundle together loci for deletions and amplifications that are synonymous
    lociCNAgenes = {}
    lociCNA = pd.DataFrame(columns=d1.columns)
    print('Bundling amplification loci...')
    for loci1 in ampLoci:
        # Get matrix of CNAs for genes in loci
        dt = posD1.loc[list(set(posD1.index).intersection(ampLoci[loci1]))]
        dt = dt[dt.columns.difference(blacklist)]
        # Get unique rows
        dedup = dt.drop_duplicates(keep='first', ignore_index=True)
        # Get genes which match and add to output dictionaries
        for i in range(len(dedup.index)):
            cnaName = loci1+'_'+str(i)+'_CNAamp'
            lociCNA.loc[cnaName] = dedup.iloc[i]
            lociCNAgenes[cnaName] = [j for j in dt.index if dedup.iloc[i].equals(dt.loc[j])]
# %%
    print('Bundling deletion loci...')
    for loci1 in delLoci:
        # Get matrix of CNAs for genes in loci
        dt = negD1.loc[list(set(posD1.index).intersection(delLoci[loci1]))]
        dt = dt[dt.columns.difference(blacklist)]
        # Get unique rows
        dedup = dt.drop_duplicates(keep='first', ignore_index=True)
        # Get genes which match and add to output dictionaries
        for i in range(len(dedup.index)):
            cnaName = loci1+'_'+str(i)+'_CNAdel'
            lociCNA.loc[cnaName] = dedup.iloc[i]
            lociCNAgenes[cnaName] = [j for j in dt.index if dedup.iloc[i].equals(dt.loc[j])]
# %%

    # Make combined matrix
    # LoF = deletions + somatic point mutations (+ fusions if have data)
    # Act = amplifications + somatic point mutations (+ fusions if have data)
    print('Starting somatic mutations...')
    pamLofAct = {}
    freq = {}
    potMuts = somMutPoint
    if not params['fusions_file'] == 'none':
        potMuts = potMuts.append(fusionsMut)
    for s1 in potMuts:
        if s1>0:
            if not str(s1) in pamLofAct:
                pamLofAct[str(s1)] = {}
            if s1 in somMuts.index:
                tmpSom = somMuts.loc[s1]
                tmpSomMean = tmpSom[tmpSom.index.difference(blacklist)].mean()
                # If potential PAM, store PAM
                if not (str(s1)+'_PAM' in pamLofAct[str(s1)] or sum(tmpSom)==0):
                    pamLofAct[str(s1)][str(s1)+'_PAM'] = tmpSom
            if (not params['fusions_file'] == 'none') and (s1 in fusions.index):
                tmpFusion = fusions.loc[s1]
                if not (str(s1)+'_Fusion' in pamLofAct[str(s1)] or sum(tmpFusion)==0):
                    pamLofAct[str(s1)][str(s1)+'_Fusion'] = tmpFusion
                if s1 in somMuts.index:
                    tmpSom = tmpSom.add(tmpFusion)
                    tmpSom[tmpSom > 1] = 1
                else:
                    tmpSom = tmpFusion
                    tmpSom[tmpSom > 1] = 1
            if (s1 in negD1.index and s1 in posD1.index):
                tmpNeg = negD1.loc[s1]
                tmpLoF = tmpSom.add(tmpNeg)[tmpNeg.index]
                tmpLoF[tmpLoF > 1] = 1
                tmpPos = posD1.loc[s1]
                tmpAct = tmpSom.add(tmpPos)[tmpPos.index]
                tmpAct[tmpAct > 1] = 1
                if not s1 in freq:
                    if (not params['fusions_file'] == 'none') and (s1 in fusions.index):
                        freq[str(s1)] = {'PAM':tmpSomMean,'Fusion':tmpFusion[tmpFusion.index.difference(blacklist)].mean(),'CNAdel':tmpNeg[tmpNeg.index.difference(blacklist)].mean(),'CNAamp':tmpPos[tmpPos.index.difference(blacklist)].mean(),'LoF':tmpLoF[tmpLoF.index.difference(blacklist)].mean(),'Act':tmpAct[tmpAct.index.difference(blacklist)].mean()}
                    else:
                        freq[str(s1)] = {'PAM':tmpSomMean,'CNAdel':tmpNeg[tmpNeg.index.difference(blacklist)].mean(),'CNAamp':tmpPos[tmpPos.index.difference(blacklist)].mean(),'LoF':tmpLoF[tmpLoF.index.difference(blacklist)].mean(),'Act':tmpAct[tmpAct.index.difference(blacklist)].mean()}
            else:
                if not s1 in freq:
                    if (params['fusions_file'] == 'none') and (s1 in fusions.index):
                        freq[str(s1)] = {'PAM':tmpSomMean,'CNAdel':0,'CNAamp':0,'LoF':0,'Act':0}
                    else:
                        freq[str(s1)] = {'PAM':tmpSomMean,'Fusion':0,'CNAdel':0,'CNAamp':0,'LoF':0,'Act':0}
# %%
    print('Starting amplifications...')
    for loci1 in ampLoci:
        for s1 in  set(ampLoci[loci1]).intersection(somMuts.index):
            if s1>0:
                # If potential Act
                if s1 in somMuts.index and s1 in posD1.index:
                    tmpSom = somMuts.loc[s1]
                    tmpSomMean = tmpSom[tmpSom.index.difference(blacklist)].mean()
                    if not params['fusions_file'] == 'none' and (s1 in fusions.index):
                        tmpFusion = fusions.loc[s1]
                        tmpSom = tmpSom.add(tmpFusion)
                        tmpSom[tmpSom > 1] = 1
                    tmpNeg = negD1.loc[s1]
                    tmpLoF = (tmpSom.add(tmpNeg)[tmpNeg.index]).clip(0,1)
                    tmpPos = posD1.loc[s1]
                    tmpAct = (tmpSom.add(tmpPos)[tmpPos.index]).clip(0,1)
                    if not s1 in freq:
                        if params['fusions_file'] == 'none':
                            freq[str(s1)] = {'PAM':tmpSomMean,'Fusion':tmpFusion[tmpFusion.index.difference(blacklist)].mean(),'CNAdel':tmpNeg[tmpNeg.index.difference(blacklist)].mean(),'CNAamp':tmpPos[tmpPos.index.difference(blacklist)].mean(),'LoF':tmpLoF[tmpLoF.index.difference(blacklist)].mean(),'Act':tmpAct[tmpAct.index.difference(blacklist)].mean()}
                        else:
                            freq[str(s1)] = {'PAM':tmpSomMean,'CNAdel':tmpNeg[tmpNeg.index.difference(blacklist)].mean(),'CNAamp':tmpPos[tmpPos.index.difference(blacklist)].mean(),'LoF':tmpLoF[tmpLoF.index.difference(blacklist)].mean(),'Act':tmpAct[tmpAct.index.difference(blacklist)].mean()}
                    # Store Act
                    if not str(s1) in pamLofAct:
                        pamLofAct[str(s1)] = {}
                    if not (str(s1)+'_Act' in pamLofAct[str(s1)] or tmpAct[tmpAct.index.difference(blacklist)].sum() == tmpSom[tmpSom.index.difference(blacklist)].sum() or tmpAct[tmpAct.index.difference(blacklist)].sum() == tmpPos[tmpPos.index.difference(blacklist)].sum()):
                        pamLofAct[str(s1)][str(s1)+'_Act'] = tmpAct
                        pamLofAct[str(s1)][str(s1)+'_CNAamp'] = tmpPos
                    elif not (str(s1)+'_CNAamp' in pamLofAct[str(s1)] or tmpAct[tmpAct.index.difference(blacklist)].sum() == tmpSom[tmpSom.index.difference(blacklist)].sum()): 
                        pamLofAct[str(s1)][str(s1)+'_CNAamp'] = tmpPos
        for s1 in set(ampLoci[loci1]).difference(somMuts.index):
            if s1>0:
                tmpNeg = negD1.loc[s1]
                tmpPos = posD1.loc[s1]
                if not s1 in freq:
                    if (params['fusions_file'] == 'none') and (s1 in fusions.index):
                        freq[str(s1)] = {'PAM':np.nan,'CNAdel':tmpNeg[tmpNeg.index.difference(blacklist)].mean(),'CNAamp':tmpPos[tmpPos.index.difference(blacklist)].mean(),'LoF':np.nan,'Act':np.nan}
                    else:
                        freq[str(s1)] = {'PAM':np.nan,'Fusion':np.nan,'CNAdel':tmpNeg[tmpNeg.index.difference(blacklist)].mean(),'CNAamp':tmpPos[tmpPos.index.difference(blacklist)].mean(),'LoF':np.nan,'Act':np.nan}
                # Store Amp
                if not str(s1) in pamLofAct:
                    pamLofAct[str(s1)] = {}
                if not str(s1)+'_CNAamp' in pamLofAct[str(s1)]:
                    pamLofAct[str(s1)][str(s1)+'_CNAamp'] = tmpPos
    # %%
    print('Starting deletions...')
    for loci1 in delLoci:
        for s1 in set(delLoci[loci1]).intersection(somMuts.index):
            if s1>0:
                # If potential Lof
                if s1 in somMuts.index and s1 in negD1.index:
                    tmpSom = somMuts.loc[s1]
                    tmpSomMean = tmpSom[tmpSom.index.difference(blacklist)].mean()
                    if not params['fusions_file'] == 'none' and (s1 in fusions.index):
                        tmpFusion = fusions.loc[s1]
                        tmpSom = tmpSom.add(tmpFusion)
                        tmpSom[tmpSom > 1] = 1
                    tmpNeg = negD1.loc[s1]
                    tmpLoF = tmpSom.add(tmpNeg)[tmpNeg.index]
                    tmpLoF[tmpLoF > 1] = 1
                    tmpPos = posD1.loc[s1]
                    tmpAct = tmpSom.add(tmpPos)[tmpPos.index]
                    tmpAct[tmpAct > 1] = 1
                    if not s1 in freq:
                        if params['fusions_file'] == 'none':
                            freq[str(s1)] = {'PAM':tmpSomMean,'Fusion':tmpFusion[tmpFusion.index.difference(blacklist)].mean(),'CNAdel':tmpNeg[tmpNeg.index.difference(blacklist)].mean(),'CNAamp':tmpPos[tmpPos.index.difference(blacklist)].mean(),'LoF':tmpLoF[tmpLoF.index.difference(blacklist)].mean(),'Act':tmpAct[tmpAct.index.difference(blacklist)].mean()}
                        else:
                            freq[str(s1)] = {'PAM':tmpSomMean,'CNAdel':tmpNeg[tmpNeg.index.difference(blacklist)].mean(),'CNAamp':tmpPos[tmpPos.index.difference(blacklist)].mean(),'LoF':tmpLoF[tmpLoF.index.difference(blacklist)].mean(),'Act':tmpAct[tmpAct.index.difference(blacklist)].mean()}
                    # Store LoF
                    if not str(s1) in pamLofAct:
                        pamLofAct[str(s1)] = {}
                    if not (str(s1)+'_LoF' in pamLofAct[str(s1)] or tmpLoF[tmpLoF.index.difference(blacklist)].sum() == tmpSom[tmpSom.index.difference(blacklist)].sum() or tmpLoF[tmpLoF.index.difference(blacklist)].sum() == tmpNeg[tmpNeg.index.difference(blacklist)].sum()):
                        pamLofAct[str(s1)][str(s1)+'_LoF'] = tmpLoF
                        pamLofAct[str(s1)][str(s1)+'_CNAdel'] = tmpNeg
        for s1 in set(delLoci[loci1]).difference(somMuts.index):
            if s1>0:
                tmpNeg = negD1.loc[s1]
                tmpPos = posD1.loc[s1]
                if not s1 in freq:
                    if (params['fusions_file'] == 'none') and (s1 in fusions.index):
                        freq[str(s1)] = {'PAM':np.nan,'CNAdel':tmpNeg[tmpNeg.index.difference(blacklist)].mean(),'CNAamp':tmpPos[tmpPos.index.difference(blacklist)].mean(),'LoF':np.nan,'Act':np.nan}
                    else:
                        freq[str(s1)] = {'PAM':np.nan,'Fusion':np.nan,'CNAdel':tmpNeg[tmpNeg.index.difference(blacklist)].mean(),'CNAamp':tmpPos[tmpPos.index.difference(blacklist)].mean(),'LoF':np.nan,'Act':np.nan}
                # Store Del
                if not str(s1) in pamLofAct:
                    pamLofAct[str(s1)] = {}
                if not str(s1)+'_CNAdel' in pamLofAct[str(s1)]:
                    pamLofAct[str(s1)][str(s1)+'_CNAdel'] = tmpNeg
    # %%
    # Decide which mutations will be tested in permutation analysis
    print('Screening for frequency...')
    keepPAM = []
    keepFusion = []
    keepDel = []
    keepAmp = []
    keepers = {}
    calcSig = []
    for s1 in pamLofAct:
        if s1 in freq:
            freqPAM = freq[s1]['PAM']
            if 'Fusion' in freq[s1]:
                freqFusion = freq[s1]['Fusion']
            freqPos = freq[s1]['CNAamp']
            freqNeg = freq[s1]['CNAdel']
            #if summaryMatrix.loc[int(s1),'CNA_type']=='Del':
            #    summaryMatrix.loc[int(s1), 'CNA_freq'] = freq[s1]['CNAdel']
            #elif summaryMatrix.loc[int(s1),'CNA_type']=='Amp':
            #    summaryMatrix.loc[int(s1), 'CNA_freq'] = freq[s1]['CNAamp']
            freqAct = freq[s1]['Act']
            summaryMatrix.loc[int(s1), 'Act_freq'] = freq[s1]['Act']
            freqLoF = freq[s1]['LoF']
            summaryMatrix.loc[int(s1), 'LoF_freq'] = freq[s1]['LoF']
            if freqLoF>=params['min_mut_freq'] or freqAct>=params['min_mut_freq'] or freqPAM>=params['min_mut_freq'] or freqPos>=params['min_mut_freq'] or freqNeg>=params['min_mut_freq'] or ('Fusion' in freq[s1] and freqFusion>=params['min_mut_freq']):
                name1 = 'Unknown'
                if sum(n1.isin([int(s1)]))==1:
                    name1 = n1.index[n1==int(s1)][0]
                if 'Fusion' in freq[s1]:
                    print('\t'+''.join([str(i) for i in [name1+' ('+str(s1),') - FreqPAM: ', round(freqPAM,3), ' | FreqFusion: ', round(freqFusion,3), ' | FreqNeg: ', round(freqNeg,3), ' | FreqLoF: ', round(freqLoF,3), ' | FreqPos: ', round(freqPos,3),' | FreqAct: ', round(freqAct,3)]]))
                else:
                    print('\t'+''.join([str(i) for i in [name1+' ('+str(s1),') - FreqPAM: ', round(freqPAM,3), ' | FreqNeg: ', round(freqNeg,3), ' | FreqLoF: ', round(freqLoF,3), ' | FreqPos: ', round(freqPos,3),' | FreqAct: ', round(freqAct,3)]]))
            # Add PAM
            if freqPAM>0 and freqPAM>=params['min_mut_freq'] and int(s1) in somMutPoint and int(s1) in sigPAMs:
                keepers[str(s1)+'_PAM'] = pamLofAct[str(s1)][str(s1)+'_PAM']
                keepPAM.append(str(s1)+'_PAM')
                summaryMatrix.loc[int(s1),'OM_type_selected'] = 'PAM'
            # Add Fusion
            if 'Fusion' in freq[s1]:
                if freqFusion>0 and freqFusion>=params['min_mut_freq'] and int(s1) in fusionsMut:
                    keepers[str(s1)+'_Fusion'] = pamLofAct[str(s1)][str(s1)+'_Fusion']
                    keepFusion.append(str(s1)+'_Fusion')
                    summaryMatrix.loc[int(s1),'OM_type_selected'] = 'Fusion'
            # Add Act
            if str(s1)+'_Act' in pamLofAct[str(s1)] and freqAct>freqPAM and freqAct>=params['min_mut_freq'] and freqAct>freqLoF:
                if freqPAM>=params['min_pam_freq']:
                    keepers[str(s1)+'_Act'] = pamLofAct[str(s1)][str(s1)+'_Act']
                    calcSig.append(str(s1)+'_Act')
                    summaryMatrix.loc[int(s1),'OM_type_selected'] = 'Act'
            # Add LoF
            if str(s1)+'_LoF' in pamLofAct[str(s1)] and freqLoF>freqPAM and freqLoF>=params['min_mut_freq'] and freqLoF>freqAct:
                if freqPAM>=params['min_pam_freq']:
                    keepers[str(s1)+'_LoF'] = pamLofAct[str(s1)][str(s1)+'_LoF']
                    calcSig.append(str(s1)+'_LoF')
                    summaryMatrix.loc[int(s1),'OM_type_selected'] = 'LoF'
            # Add CNAamp
            if (str(s1)+'_CNAamp' in pamLofAct[str(s1)]) and (freqPos>=params['min_mut_freq']) and (not ((summaryMatrix.loc[int(s1),'OM_type_selected']=='Act') or (summaryMatrix.loc[int(s1),'OM_type_selected']=='LoF'))):
                #print(str(s1)+'_CNAamp')
                keepers[str(s1)+'_CNAamp'] = pamLofAct[str(s1)][str(s1)+'_CNAamp']
                keepAmp.append(str(s1)+'_CNAamp')
                summaryMatrix.loc[int(s1),'OM_type_selected'] = 'CNAamp'
            # Add CNAdel
            if ((str(s1)+'_CNAdel' in pamLofAct[str(s1)]) and (freqNeg>=params['min_mut_freq'])) and (not ((summaryMatrix.loc[int(s1),'OM_type_selected']=='Act') or (summaryMatrix.loc[int(s1),'OM_type_selected']=='LoF'))):
                #print(str(s1)+'_CNAdel')
                keepers[str(s1)+'_CNAdel'] = pamLofAct[str(s1)][str(s1)+'_CNAdel']
                keepDel.append(str(s1)+'_CNAdel')
                summaryMatrix.loc[int(s1),'OM_type_selected'] = 'CNAdel'

    # %%
    ##################################
    ## Conduct permutation analysis ##
    ##################################
    print('Permutation anlaysis...')

    # %%
    numPermutes = 1000
    # Permute to get frequency
    def singlePermute(somMutsMF, somFusionMF, somCNAsMF):
        perms = []
        tmp1 = np.array(somMutsMF.loc[np.random.choice(somMutsMF.index, 1)])
        tmp2 = np.array(somCNAsMF.loc[np.random.choice(somCNAsMF.index, 1)])
        if type(somFusionMF)!=str:
            tmp3 = np.array(somFusionMF.loc[np.random.choice(somFusionMF.index, 1)])
            temp = (tmp1 + tmp2 + tmp3).clip(0,1)
            perms = temp.mean()
        else:
            temp = (tmp1 + tmp2).clip(0,1)
            perms = temp.mean()
        return perms        
    # %%
    permMF_neg = []
    permMF_pos = []
    if params['load_permutation']==None:
        ## Compute permutations if not loading from previous run
        # Deletions
        print('\tPermuting deletions...')
        permutedGenes = (set(sigPAMs) | ampGenesSet | delGenesSet) & set(list(posD1.index)) & set(list(negD1.index)) & set(somMuts.index)
        if not params['fusions_file'] == 'none':
            permutedGenes = permutedGenes & set(fusions.index)
        permutedGenes = np.array(list(permutedGenes))
        somMutsMF = somMuts[somMuts.index<1000000000][somMuts.columns.difference(blacklist)]
        somMutsMF = somMutsMF.loc[permutedGenes[permutedGenes<1000000000],:]
        if not params['fusions_file'] == 'none':
            somFusionMF = fusions[fusions.index<1000000000][fusions.columns.difference(blacklist)]
            somFusionMF = somFusionMF.loc[permutedGenes[permutedGenes<1000000000],:]
        else:
            somFusionMF = 'none'
        somCNAsNegMF = negD1[negD1.index<1000000000][negD1.columns.difference(blacklist)]
        somCNAsNegMF = somCNAsNegMF.loc[permutedGenes[permutedGenes<1000000000],:]
        with tqdm(total=numPermutes*10) as pbar:
            for i in range(numPermutes*10):
                permMF_neg.append(singlePermute(somMutsMF, somFusionMF, somCNAsNegMF))
                pbar.update(1)

        # Amplifications
        print('\tPermuting amplifications...')
        somCNAsPosMF = posD1[posD1.index<1000000000][posD1.columns.difference(blacklist)]
        somCNAsPosMF = somCNAsPosMF.loc[permutedGenes[permutedGenes<1000000000],:]
        with tqdm(total=numPermutes*10) as pbar:
            for i in range(numPermutes*10):
                permMF_pos.append(singlePermute(somMutsMF, somFusionMF, somCNAsPosMF))
                pbar.update(1)

        # Change precision so that is the same as what will be written out,
        # - Fixes bug where precision change leads to different behavior from freshly run (then saved) and loaded permutation data
        permMF_neg = [float(str(i)) for i in permMF_neg]
        permMF_pos = [float(str(i)) for i in permMF_pos]

        # If requested to save out permutations
        if params['save_permutation']==True:
            print('\tSaving permutations...')
            np.save(params['output_path']+'/oncomerge_delPerm', permMF_neg)
            np.save(params['output_path']+'/oncomerge_ampPerm', permMF_pos)
    elif not params['load_permutation']==None:
        ## Load up permutations from previous run
        print('\tLoading previous permutations...')
        permMF_neg = np.load(params['load_permutation']+'/oncomerge_delPerm.npy')
        permMF_pos = np.load(params['load_permutation']+'/oncomerge_ampPerm.npy')

    # Write Permutation Analysis file Lof_Act_sig
    lofActSig = pd.DataFrame(columns = ['Symbol', 'Type','Freq','Emp.p_value'], index = calcSig)
    for sig1 in calcSig:
        if sig1.find('LoF')>0:
            lofActSig['Symbol'].loc[sig1] = n1.index[n1==int(sig1.rstrip('_LoF'))][0]
            lofActSig['Type'].loc[sig1] = 'LoF'
            lofActSig['Freq'].loc[sig1] = freq[sig1.rstrip('_LoF')]['LoF']
        elif sig1.find('Act')>0:
            lofActSig['Symbol'].loc[sig1] = n1.index[n1==int(sig1.rstrip('_Act'))][0]
            lofActSig['Type'].loc[sig1] = 'Act'
            lofActSig['Freq'].loc[sig1] = freq[sig1.rstrip('_Act')]['Act']


##############################################################################################################################
    # %%
    # Pathway Permutations
    def singlePathwayPermute(somMuts, somFusions, somCNAs, numPermutes):
        with Pool(processes=7) as p: 
            permuts =  list(tqdm(p.imap(pathwayPool, [[i, somMuts, somFusions, somCNAs, numPermutes] for i in pathways["length"].unique()]), total=len(pathways["length"].unique())))
        perm = dict()
        for i in permuts:
            perm[i[0]] = i[1]
        p.terminate()
        return perm
    print('\tPermuting pathway deletions...')  
    permPath_neg = singlePathwayPermute(somMutsMF, somFusionMF, somCNAsNegMF, numPermutes)
    
    # %%
    print('\tPermuting pathway amplifications...')
    permPath_pos = singlePathwayPermute(somMutsMF, somFusionMF, somCNAsPosMF, numPermutes)
    # %%
    permPath_neg = {j:[float(str(i)) for i in permPath_neg[j]] for j in permPath_neg}
    permPath_pos = {j:[float(str(i)) for i in permPath_pos[j]] for j in permPath_pos}

    # %%
    print('\tCalculating p_values...')

    # %%
    # Precalculate the permuted p-values for each frequency
    permdict1 = {}
    permdict1['LoF'] = {}
    permdict1['Act'] = {}
    oMfreqs = {sig1:freq[sig1.split('_')[0]][sig1.split('_')[1]] for sig1 in calcSig}

    # %%
    pathways["PathwayID"].index = pathways.index
    # %%
    with tqdm(total=len(oMfreqs)) as pbar:
        for f in oMfreqs.keys():
            if int(f.split('_')[0])<1000000000:
                permdict1['LoF'][f] = float(len([i for i in permMF_neg if i >= oMfreqs[f]]))/len(permMF_neg)
                permdict1['Act'][f] = float(len([i for i in permMF_pos if i >= oMfreqs[f]]))/len(permMF_pos)
            else:
                perm_neg = permPath_neg[str(pathways.loc[pathways["PathwayID"] == int(f.split('_')[0]), 'length'].iloc[0])]
                perm_pos = permPath_pos[str(pathways.loc[pathways["PathwayID"] == int(f.split('_')[0]), 'length'].iloc[0])]
                permdict1['LoF'][f] = float(len([i for i in perm_neg if i >= oMfreqs[f]]))/len(perm_neg)
                permdict1['Act'][f] = float(len([i for i in perm_pos if i >= oMfreqs[f]]))/len(perm_pos)
            pbar.update(1)

    # %%
    # Add permuted p-values to summary matrix and permutation summary
    for sig1 in calcSig:
        if sig1.find('LoF')>0:
            lofActSig.loc[sig1, 'Emp.p_value'] = permdict1['LoF'][sig1]
            summaryMatrix.loc[int(sig1.rstrip('_LoF')),'OM_empirical_p_value'] = permdict1['LoF'][sig1]
        elif sig1.find('Act')>0:
            lofActSig.loc[sig1, 'Emp.p_value'] = permdict1['Act'][sig1]
            summaryMatrix.loc[int(sig1.rstrip('_Act')),'OM_empirical_p_value'] = permdict1['Act'][sig1]
    # %%
    # Filter LoF and Act based on permuted p-values
    lofActSig["q_value"]=np.empty((lofActSig.shape[0],1))
    if len(lofActSig)>0:
        lofActSig.loc[np.array(lofActSig.index.str.split('_').str[0]).astype(int)<1000000000, 'q_value'] = multipletests(lofActSig[np.array(lofActSig.index.str.split('_').str[0]).astype(int)<1000000000]['Emp.p_value'], 0.05, method='fdr_bh')[1]
        lofActSig.loc[np.array(lofActSig.index.str.split('_').str[0]).astype(int)>1000000000, 'q_value'] = multipletests(lofActSig[np.array(lofActSig.index.str.split('_').str[0]).astype(int)>1000000000]['Emp.p_value'], 0.05, method='fdr_bh')[1]
        summaryMatrix.loc[[int(i.split('_')[0]) for i in lofActSig.index],'OM_empirical_q_value'] = list(lofActSig['q_value'])
        lofActSig.sort_values('q_value').to_csv(params['output_path']+'/oncoMerge_ActLofPermPV.csv')
        # Screen out LoF and Act that don't meet significance cutoffs
        keepLofAct0 = list(lofActSig.index[np.logical_or(lofActSig['q_value']<=float(params['perm_qv']), lofActSig['Emp.p_value']<=0.01)])
    else:
        # No LoF or Act to filter
        lofActSig.to_csv(params['output_path']+'/oncoMerge_ActLofPermPV.csv')
        keepLofAct0 = []
    # %%
    # Function to map mutation to locus
    def findLoci(mutation, ampLoci, delLoci):
        gene, mutType = mutation.split('_')
        if (mutType == 'Act') or (mutType == 'CNAamp'):
            loci = ampLoci
        if (mutType == 'LoF') or (mutType == 'CNAdel'):
            loci = delLoci
        return [locus1 for locus1 in loci.keys() if int(gene) in loci[locus1]]

    # Tabulate the number of genes per locus Lof Act
    lofActLoci = {}
    # %%
    for mut in keepLofAct0:
        mutLoci = findLoci(mut, ampLoci, delLoci)
        for locus1 in mutLoci:
            if locus1 not in lofActLoci.keys():
                lofActLoci[locus1] = []
            lofActLoci[locus1].append(mut)
    # %%
    # Tabulate the number of genes per locus CNAs
    combinedLoci = lofActLoci.copy()
    # %%
    for mut in keepDel+keepAmp:
        mutLoci = findLoci(mut, ampLoci, delLoci)
        for locus1 in mutLoci:
            if not locus1 in lofActLoci:
                if locus1 not in combinedLoci.keys():
                    combinedLoci[locus1] = []
                combinedLoci[locus1].append(mut)
    # %%
    ## Decide whether to apply the maximum final frequency filter
    keepLofAct1 = []
    for locus in combinedLoci.keys():
        # Don't filter further with maximum final frequency filter
        if len(combinedLoci[locus]) < params['min_loci_genes']:
            keepLofAct1 += combinedLoci[locus]
            for mut in combinedLoci[locus]:
                gene, mutType = mut.split('_')
                summaryMatrix.loc[int(gene), 'Final_mutation_type'] = mutType
                if mutType=='Act' or mutType=='LoF':
                    summaryMatrix.loc[int(gene), 'Final_freq'] = summaryMatrix.loc[int(gene), mutType+'_freq']
                    summaryMatrix.loc[int(gene), 'Delta_over_PAM'] = summaryMatrix.loc[int(gene), mutType+'_freq'] - summaryMatrix.loc[int(gene), 'PAM_freq']
                else:
                    summaryMatrix.loc[int(gene), 'Final_freq'] = summaryMatrix.loc[int(gene), 'CNA_freq']
                summaryMatrix.loc[int(gene), 'Genes_in_locus'] = len(combinedLoci[locus])
        # Filter with maximum final frequency filter
        else:
            LociFF = [summaryMatrix.loc[int(mut.split('_')[0]), mut.split('_')[1]+'_freq'] if (mut.split('_')[1]=='Act' or mut.split('_')[1]=='LoF') else summaryMatrix.loc[int(mut.split('_')[0]), 'CNA_freq'] for mut in combinedLoci[locus]]
            maxFF = max(LociFF)
            gl1 = pd.Series(np.array(LociFF) == maxFF).sum()
            for mut in combinedLoci[locus]:
                gene, mutType = mut.split('_')
                keepLofAct1.append(mut)
                if ((mutType=='Act' or mutType=='LoF') and (summaryMatrix.loc[int(gene), mutType+'_freq']==maxFF)) or summaryMatrix.loc[int(gene), 'CNA_freq']==maxFF:
                    summaryMatrix.loc[int(gene), 'Genes_in_locus'] = gl1
                    summaryMatrix.loc[int(gene), 'Final_mutation_type'] = mutType
                    if mutType=='Act' or mutType=='LoF':
                        summaryMatrix.loc[int(gene), 'Final_freq'] = summaryMatrix.loc[int(gene), mutType+'_freq']
                        summaryMatrix.loc[int(gene), 'Delta_over_PAM'] = summaryMatrix.loc[int(gene), mutType+'_freq'] - summaryMatrix.loc[int(gene), 'PAM_freq']
                    else:
                        summaryMatrix.loc[int(gene), 'Final_freq'] = summaryMatrix.loc[int(gene), 'CNA_freq']

    keepLofAct = list(set(keepLofAct1))
    # %%
    # Screen out PAMs that are LoF/Act
    newKeepPAM = []
    for pam1 in keepPAM:
        found = 0
        tmp1 = pam1.split('_')[0]
        for lofAct in keepLofAct:
            if tmp1==lofAct.split('_')[0]:
                found = 1
        if found==0:
            newKeepPAM.append(pam1)
            summaryMatrix.loc[int(tmp1), 'Final_mutation_type'] = 'PAM'
            summaryMatrix.loc[int(tmp1), 'Final_freq'] = summaryMatrix.loc[int(tmp1), 'PAM_freq']
            summaryMatrix.loc[int(tmp1), 'Delta_over_PAM'] = float(0)

    # Screen out Fusions that are LoF/Act
    for fus1 in keepFusion:
        found = 0
        tmp1 = fus1.split('_')[0]
        for lofAct in keepLofAct:
            if tmp1==lofAct.split('_')[0]:
                found = 1
        if found==0:
            newKeepPAM.append(fus1)
            summaryMatrix.loc[int(tmp1), 'Final_mutation_type'] = 'Fusion'
            summaryMatrix.loc[int(tmp1), 'Final_freq'] = summaryMatrix.loc[int(tmp1), 'Fusion_freq']
            summaryMatrix.loc[int(tmp1), 'Delta_over_PAM'] = summaryMatrix.loc[int(tmp1), 'Fusion_freq']-summaryMatrix.loc[int(tmp1), 'PAM_freq']

    ## Screen out loci that have a representative gene
    # Mutations that are at or above minimum mutation frequency cutoff
    highFreqLoci = lociCNA.loc[lociCNA.mean(axis=1)>=params['min_mut_freq']]
    highFreqLoci = highFreqLoci[~highFreqLoci.index.str.split('_').str[0].isin(pathways["Pathway"])]
    # %%
    # Figure out what loci are explained by current Act or LoF genes
    explainedLoc = []
    keepLoc = []
    AmpDelLoci = {**ampLoci, **delLoci}
    # %%
    for locus1 in highFreqLoci.index:
        genesInLocus = [i for i in lociCNAgenes[locus1] if (str(i)+'_Act' in keepLofAct or str(i)+'_LoF' in keepLofAct or str(i)+'_CNAamp' in keepLofAct or str(i)+'_CNAdel' in keepLofAct)]
        # ActLoFInLocus = [i for i in AmpDelLoci[locus1.split('_')[0]] if (str(i)+'_Act' in keepLofAct or str(i)+'_LoF' in keepLofAct)]
        if len(genesInLocus)>0:
            explainedLoc.append(locus1.split('_')[0])
        else:
            keepLoc.append(locus1)
            
    # %%
    keepLoc_dict = {}
    for locus1 in set(['_'.join([i.split('_')[0],i.split('_')[-1]]) for i in keepLoc]):
        locus2 = locus1.split('_')[0]
        if locus2 in ampLoci.keys():
            keepLoc_dict[locus1] = posD1.loc[ampLoci[locus2]]
        if locus2 in delLoci.keys():
            keepLoc_dict[locus1] = negD1.loc[delLoci[locus2]]

    def mode2(pat1col):
        tmpser = pat1col.mode()
        if len(tmpser)>1:
            tmp10 = 1
        else:
            tmp10 = tmpser.iloc[0]
        return tmp10

    # Use most common value to determine the mutational profile of a loci
    keepLoc_df = pd.DataFrame(columns = d1.columns)
    tmpMode = pd.Series(index = d1.columns)
    for locus1 in keepLoc_dict.keys():
        for pat1 in keepLoc_dict[locus1].columns:
            tmpMode.loc[pat1] = mode2(keepLoc_dict[locus1][pat1])
        if tmpMode.mean()>=0.05:
            keepLoc_df.loc[locus1] = tmpMode

    keepLoc_df = keepLoc_df.applymap(int)
    # %%
    ####################################
    ## Compile OncoMerge output files ##
    ####################################
    finalMutFile = pd.concat([pd.DataFrame(keepers).transpose().loc[newKeepPAM].sort_index(), pd.DataFrame(keepers).transpose().loc[keepLofAct].sort_index(), keepLoc_df.sort_index()], sort=True)
    # %%
    # Rename all loci with only one gene
    ind_list = finalMutFile.index.tolist()
    for locus1 in keepLoc_df.index:
        splitUp = locus1.split('_')
        if splitUp[-1]=='CNAamp' and len(ampLoci[splitUp[0]])==1:
            idx = ind_list.index(locus1)
            ind_list[idx] = str(ampLoci[splitUp[0]][0]) + '_CNAamp'
            summaryMatrix.loc[int(ampLoci[splitUp[0]][0]), 'Final_mutation_type'] = 'CNAamp'
            summaryMatrix.loc[int(ampLoci[splitUp[0]][0]), 'Final_freq'] = summaryMatrix.loc[int(ampLoci[splitUp[0]][0]), 'CNA_freq']
            summaryMatrix.loc[int(ampLoci[splitUp[0]][0]), 'Delta_over_PAM'] = summaryMatrix.loc[int(ampLoci[splitUp[0]][0]), 'CNA_freq']

        if splitUp[-1]=='CNAdel' and len(delLoci[splitUp[0]])==1:
            idx = ind_list.index(locus1)
            ind_list[idx] = str(delLoci[splitUp[0]][0]) + '_CNAdel'
            summaryMatrix.loc[int(delLoci[splitUp[0]][0]), 'Final_mutation_type'] = 'CNAdel'
            summaryMatrix.loc[int(delLoci[splitUp[0]][0]), 'Final_freq'] = summaryMatrix.loc[int(delLoci[splitUp[0]][0]), 'CNA_freq']
            summaryMatrix.loc[int(delLoci[splitUp[0]][0]), 'Delta_over_PAM'] = summaryMatrix.loc[int(delLoci[splitUp[0]][0]), 'CNA_freq']
    # %%
    ####################################
    ## ME Analysis                    ##
    ####################################
    print("Starting Mutual Exclusivity Analysis...")
    import rpy2
    import rpy2.robjects as ro
    from rpy2.robjects.vectors import DataFrame
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.conversion import localconverter

    # add fusions

    ro.r('''
        # create a function `f`
        MEanalysis <- function(somMuts, posD1, negD1, Act, LoF, pathways){
            library("Rediscover")
            library("tidyverse")
            library("discover")
            out = data.frame(names=character(), ME_pvalues=double())
            somPM <- getPM(data.matrix(somMuts))
            posPM <- getPM(data.matrix(posD1))
            negPM <- getPM(data.matrix(negD1))
            actPM <- getPM(data.matrix(Act))
            lofPM <- getPM(data.matrix(LoF))
            if (nrow(pathways) > 0){
                for (row in 1:nrow(pathways)) {
                    if (pathways[row, "Final_mutation_type"] == "PAM"){
                        muts = somMuts
                        PM = somPM
                    }
                    if (pathways[row, "Final_mutation_type"] == "CNAamp"){
                        muts = posD1
                        PM = posPM
                    }
                    if (pathways[row, "Final_mutation_type"] == "CNADel"){
                        muts = negD1
                        PM = negPM
                    }
                    if (pathways[row, "Final_mutation_type"] == "Act"){
                        muts = Act
                        PM = actPM
                    }
                    if (pathways[row, "Final_mutation_type"] == "LoF"){
                        muts = LoF
                        PM = lofPM
                    }
                    if (nchar(pathways[row, "Somatically_Mutated_Genes"])>2){
                        genes <-as.vector(strsplit(substring(pathways[row, "Somatically_Mutated_Genes"], 2,nchar(pathways[row, "Somatically_Mutated_Genes"])-1), ", "))[[1]]
                        genes <- genes[genes %in% rownames(muts)]
                        A = data.matrix(muts)
                        x= data.frame(rownames(muts))
                        x$indices = rownames(x)
                        p_val = getMutexGroup(data.matrix(A[genes,]), PM[as.numeric(c(filter(x, rownames.muts. %in% genes)$indices)),], "Coverage")
                        out[nrow(out) + 1,] = c(pathways[row, "Pathway_Name"], p_val)
                    }
                }
            }
            return(out)
            }
        ''')
    MEanalysis_r = ro.globalenv['MEanalysis']
    MEmuts = (somMuts.loc[list(set(somMuts.index).intersection(set(fusions.index))),:] + fusions.loc[list(set(somMuts.index).intersection(set(fusions.index))),:]).clip(0,1).astype(int)
    MEpos = posD1.copy().astype(int)
    MEneg = negD1.copy().astype(int)
    MEact = (MEmuts.loc[list(set(MEmuts.index).intersection(set(posD1.index))),:] + posD1.loc[list(set(MEmuts.index).intersection(set(posD1.index))),:]).clip(0,1).astype(int)
    MElof = (MEmuts.loc[list(set(MEmuts.index).intersection(set(negD1.index))),:] + negD1.loc[list(set(MEmuts.index).intersection(set(negD1.index))),:]).clip(0,1).astype(int)
    #with localconverter(ro.default_converter + pandas2ri.converter):
    #    MEmuts_r = ro.conversion.py2rpy(MEmuts)
    #print("Converted Somatic Mutations.")
    # %%
    MEpathways = summaryMatrix.copy()
    MEpathways = MEpathways.loc[MEpathways["Pathway_Length"]>1,]
    MEpathways = MEpathways.loc[~MEpathways["OM_empirical_p_value"].isna(),]
    MEpathways = MEpathways.loc[MEpathways["OM_empirical_p_value"]<0.01]
    MEpathways = MEpathways.loc[~MEpathways["Final_mutation_type"].isna(),]
    MEpathways = MEpathways.loc[MEpathways["Final_mutation_type"] != "",] 
    MEpathways = MEpathways.loc[:,["Somatically_Mutated_Genes", "Pathway_Name", "Final_mutation_type"]]
    MEpathways["Somatically_Mutated_Genes"] = MEpathways["Somatically_Mutated_Genes"].astype(str)
    #with localconverter(ro.default_converter + pandas2ri.converter):
    #    pathways_r = ro.conversion.py2rpy(MEpathways)
    #print("Converted Pathways.")
    # %%
    with localconverter(ro.default_converter + pandas2ri.converter):
        MEpathways_r = ro.conversion.py2rpy(MEpathways)
        MEmuts_r = ro.conversion.py2rpy(MEmuts)
        MEpos_r = ro.conversion.py2rpy(MEpos)
        MEneg_r = ro.conversion.py2rpy(MEneg)
        MEact_r = ro.conversion.py2rpy(MEact)
        MElof_r = ro.conversion.py2rpy(MElof)
    # %%
    print("Finished Converting")
    MEpvals_r = MEanalysis_r(MEmuts_r, MEpos_r, MEneg_r, MEact_r, MElof_r, MEpathways_r)
    #MEpvals_r = MEanalysis_r(MEact_r, MElof_r, MEpathways_r)
    print("Finished Running Analysis.")
    with localconverter(ro.default_converter + pandas2ri.converter):
        MEpvals = ro.conversion.rpy2py(MEpvals_r)
    # %%
    summaryMatrixInd = summaryMatrix.index
    summaryMatrix = summaryMatrix.merge(MEpvals, how="left", left_on="Pathway_Name", right_on="names")
    summaryMatrix = summaryMatrix.drop("names", axis = 1)
    summaryMatrix.index = summaryMatrixInd
    # %%
    ####################################
    ## Compile OncoMerge output files ##
    ####################################
    finalMutFile.index = ind_list

    # Write out final mutation matrix to csv file
    finalMutFile.to_csv(params['output_path']+'/oncoMerge_mergedMuts.csv')

    # Write out summary matrix to csv file
    summaryMatrix.to_csv(params['output_path']+'/oncoMerge_summaryMatrix.csv')

    ## Write out loci information
    # Prepare for writing out
    writeLoci = ['Locus_name,Genes']
    for locus1 in keepLoc_df.index:
        splitUp = locus1.split('_')
        if splitUp[-1]=='CNAamp':
            writeLoci.append(locus1+','+' '.join([str(i) for i in ampLoci[splitUp[0]]]))
        if splitUp[-1]=='CNAdel':
            writeLoci.append(locus1+','+' '.join([str(i) for i in delLoci[splitUp[0]]]))

    # Write out loci file
    with open(params['output_path']+'/oncoMerge_CNA_loci.csv','w') as outFile:
        outFile.write('\n'.join(writeLoci))

    ## Write out information for hypergeometric analysis
    # Can be used to load in gene information from final mutation file:
    #    [[int(j.split('_')[0]),j.split('_')[1]] for j in finalMutFile.index[[not i for i in list(finalMutFile.index.str.contains('p|q'))]]]

    # Write out background information for hypergeometric analysis
    backgrounds = {i:[int(j) for j in backgrounds[i]] for i in backgrounds}
    with open(params['output_path']+'/background.json', 'w') as outFile:
        json.dump(backgrounds, outFile)

    print('Done.')
# %%