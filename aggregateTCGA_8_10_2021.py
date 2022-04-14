##########################################################
## aggregateTCGA.py                                     ##
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

# %%
#####################
## Import packages ##
#####################

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42

import pandas as pd
import json
from scipy.stats import hypergeom
from scipy.stats import pearsonr
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import matplotlib.figure as figs
from matplotlib_venn import venn2, venn3
import mygene as mg
import numpy as np
import copy

###############
## Functions ##
###############

def enrichment(overlap1, background1, set1, set2):
    x = len(overlap1)
    M = len(background1)
    n = len(set1)
    N = len(set2)
    return [x, n, N, M, 1-hypergeom.cdf(x, M, n, N)]


######################
## Common variables ##
######################

colors1 = {'PAM':'tab:blue', 'Fusion':'#f46d43', 'CNAamp':'gold', 'Act':'tab:green', 'CNAdel':'tab:red', 'LoF':'tab:purple', 'Overlap':'tab:orange'}
colors2 = {'NA': 0, 'PAM':1, 'Fusion':2, 'CNAamp':3, 'Act':4, 'CNAdel':5, 'LoF':6, 'Overlap':7}
colors3 = ['white', 'tab:blue', '#f46d43', 'gold', 'tab:green', 'tab:red', 'tab:purple'] #, 'tab:orange']

tumors = ['ACC', 'BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 'DLBC', 'ESCA', 'GBM', 'HNSC', 'KICH', 'KIRC', 'KIRP', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'MESO', 'OV', 'PAAD', 'PCPG', 'PRAD', 'READ', 'SARC', 'STAD', 'TGCT', 'THCA', 'THYM', 'UCEC', 'UCS', 'UVM']
#tumors = ['ACC', 'BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD']
analyses = ['pq_mff']
#analyses = ['no_filter','pq','mff','pq_mff','mpf','pq_mpf','pq_mff_mpf']
#analyses = ['no_filter','pq','mpf','pq_mpf']
toLoad = { "TCGA_Consensus": "C:/Users/nihan/Dropbox/OncoMerge_pathways/gold_standards/TCGA_Consensus/TCGA_consensus_gold_standard.csv", # Pubmed = 29625053
           "CGC": "C:/Users/nihan/Dropbox/OncoMerge_pathways/gold_standards/CGC/CGC_gold_standard.csv",
           "Moonlight": "C:/Users/nihan/Dropbox/OncoMerge_pathways/gold_standards/Moonlight/Moonlight_gold_standard.csv",
           "OncoROLE": "C:/Users/nihan/Dropbox/OncoMerge_pathways/gold_standards/OncoROLE/OncoROLE_gold_standard.csv",
           "ONGene": "C:/Users/nihan/Dropbox/OncoMerge_pathways/gold_standards/ONGene/ONGene_gold_standard.csv",
           "TSGene": "C:/Users/nihan/Dropbox/OncoMerge_pathways/gold_standards/TSGene/TSGene_gold_standard.csv",
           "Tokheim": "C:/Users/nihan/Dropbox/OncoMerge_pathways/gold_standards/Tokheim/Tokheim_gold_standard.csv",
           "Vogelstein": "C:/Users/nihan/Dropbox/OncoMerge_pathways/gold_standards/Vogelstein/Vogelstein_gold_standard.csv" }

# %%
###############
## Load data ##
###############

print('Loading data...')
# Load up enrichment sets
enrichmentSets = {}
for set1 in toLoad:
    enSet1 = pd.read_csv(toLoad[set1], header=0)
    enrichmentSets[set1] = {}
    enrichmentSets[set1]['Aggregate'] = {'Activating': [], 'LossOfFunction': [], 'All': []}
    for cancer in set(enSet1['Cancer']):
        enrichmentSets[set1][cancer] = {}
        for type1 in ['Activating', 'LossOfFunction', 'All']:
            enrichmentSets[set1][cancer][type1] = list(enSet1.loc[(enSet1['Cancer']==cancer) & (enSet1['Type']==type1),'Entrez'])
            enrichmentSets[set1]['Aggregate'][type1] += list(enSet1.loc[(enSet1['Cancer']==cancer) & (enSet1['Type']==type1),'Entrez'])

'''
with PdfPages('enrichmentSet_Overlap_8_10_2021_by5.pdf') as pp:
    fig = plt.figure(figsize=(9.5,4))
    grid = plt.GridSpec(1,2, wspace=0.4, hspace=0.3)
    ax = plt.subplot(grid[0,1]
    activating = [set([j for i in enrichmentSets['TCGA_Consensus'] for j in enrichmentSets['TCGA_Consensus'][i]['Activating']]), set([j for i in enrichmentSets['CGC'] for j in enrichmentSets['CGC'][i]['Activating']]), set([j for i in enrichmentSets['Vogelstein'] for j in enrichmentSets['Vogelstein'][i]['Activating']]), set([j for i in enrichmentSets['OncoROLE'] for j in enrichmentSets['OncoROLE'][i]['Activating']]), set([j for i in enrichmentSets['Tokheim'] for j in enrichmentSets['Tokheim'][i]['Activating']])]
    act_labels = venn.get_labels(activating, fill=['number'])
    fig, ax = venn.venn5(act_labels, names=['TCGA Consensus','CGC','Vogelstein','OncodriveROLE','Tokheim'])
    fig.savefig(pp, format='pdf')
    plt.close()

    lossOfFunction = [set([j for i in enrichmentSets['TCGA_Consensus'] for j in enrichmentSets['TCGA_Consensus'][i]['LossOfFunction']]), set([j for i in enrichmentSets['CGC'] for j in enrichmentSets['CGC'][i]['LossOfFunction']]), set([j for i in enrichmentSets['Vogelstein'] for j in enrichmentSets['Vogelstein'][i]['LossOfFunction']]), set([j for i in enrichmentSets['OncoROLE'] for j in enrichmentSets['OncoROLE'][i]['LossOfFunction']]), set([j for i in enrichmentSets['Tokheim'] for j in enrichmentSets['Tokheim'][i]['LossOfFunction']])]
    lof_labels = venn.get_labels(lossOfFunction, fill=['number'])
    fig, ax = venn.venn5(lof_labels, names=['TCGA Consensus','CGC','Vogelstein','OncodriveROLE','Tokheim'])
    fig.savefig(pp, format='pdf')
    plt.close()
'''

'''
with PdfPages('enrichmentSet_Overlap_8_10_2021.pdf') as pp:
    fig = plt.figure(figsize=(9.5,4))
    grid = plt.GridSpec(1,2, wspace=0.4, hspace=0.3)
    ax = plt.subplot(grid[0,1])
    activating = [set([j for i in enrichmentSets['TCGA_Consensus'] for j in enrichmentSets['TCGA_Consensus'][i]['Activating']]), set([j for i in enrichmentSets['CGC'] for j in enrichmentSets['CGC'][i]['Activating']]), set([j for i in enrichmentSets['Vogelstein'] for j in enrichmentStes['Vogelstein'][i]['Activating']]), set([j for i in enrichmentSets['OncoROLE'] for j in enrichmentStes['OncoROLE'][i]['Activating']]), set([j for i in enrichmentSets['Tokheim'] for j in enrichmentStes['Tokheim'][i]['Activating']])]



    ax = plt.subplot(grid[0,1])
    venn2([set([j for i in enrichmentSets['TCGA_Consensus'] for j in enrichmentSets['TCGA_Consensus'][i]['Activating']]),set([j for i in enrichmentSets['TCGA_Consensus'] for j in enrichmentSets['TCGA_Consensus'][i]['LossOfFunction']])], ('TCGA Consensus\n(Activating)','TCGA Consensus\n(Loss of Function)'), ax=ax)
    #venn3([set([j for i in enrichmentSets['TCGA_Consensus'] for j in enrichmentSets['TCGA_Consensus'][i]['Activating']]),set([j for i in enrichmentSets['TCGA_Consensus'] for j in enrichmentSets['TCGA_Consensus'][i]['LossOfFunction']]),set(enrichmentSets['TSGene']['Any']['LossOfFunction'])], ('TCGA Consensus\n(Activating)','TCGA Consensus\n(Loss of Function)','TSGene'), ax=ax)


    ax = plt.subplot(grid[0,0])
    #venn3([set([j for i in enrichmentSets['TCGA_Consensus'] for j in enrichmentSets['TCGA_Consensus'][i]['Activating']]),set([j for i in enrichmentSets['TCGA_Consensus'] for j in enrichmentSets['TCGA_Consensus'][i]['LossOfFunction']]),set(enrichmentSets['Vogelstein']['Any']['All'])], ('TCGA Consensus\n(Activating)','TCGA Consensus\n(Loss of Function)','Vogelstein'), ax=ax)
    venn3([set(enrichmentSets['TCGA_Consensus']['Any']['All']),set(enrichmentSets['CGC']['Any']['All']), set(enrichmentSets['Vogelstein']['Any']['All'])], ('TCGA consensus','CGC','Vogelstein'), ax=ax)
    #ax.save_fig()
    fig.savefig(pp, format='pdf')
    plt.close()
'''

# Load up background information
background = {}
for tumor in tumors:
    with open('C:/Users/nihan/Dropbox/OncoMerge_pathways/output/'+tumor+'/pq_mff/background.json', 'r') as bkgd:
            background[tumor] = json.loads(bkgd.read())

# Load up gene sets for each analysis
merged = {}
oncomergeSummary = {}
geneSets = {}
summaryTable = { 'pq_mff':pd.DataFrame(index=tumors, columns=['Patients','PAM','Fusion','CNAamp','Act','CNAdel','LoF'])}
                 #'mpf':pd.DataFrame(index=tumors, columns=['Patients','PAM','CNAamp','Act','CNAdel','LoF']),
                 #'pq_mpf':pd.DataFrame(index=tumors, columns=['Patients','PAM','CNAamp','Act','CNAdel','LoF']),
                 #'pq_mff_mpf':pd.DataFrame(index=tumors, columns=['Patients','PAM','CNAamp','Act','CNAdel','LoF'])}

# %%
deltaOverPAM = []
mutFreq = pd.DataFrame(columns = ['tumor','freq'])
for tumor in tumors:
    merged[tumor] = {}
    oncomergeSummary[tumor] = {}
    geneSets[tumor] = {}
    for analysis in analyses:
        merged[tumor][analysis] = pd.read_csv('C:/Users/nihan/Dropbox/OncoMerge_pathways/output/'+tumor+'/'+analysis+'/oncoMerge_mergedMuts.csv', header = 0, index_col = 0)
        oncomergeSummary[tumor][analysis] = pd.read_csv('C:/Users/nihan/Dropbox/OncoMerge_pathways/output/'+tumor+'/'+analysis+'/oncoMerge_summaryMatrix.csv', header = 0, index_col = 0)
        geneSets[tumor][analysis] = [[int(j.split('_')[0]),j.split('_')[1]] for j in merged[tumor][analysis].index[[not i for i in list(merged[tumor][analysis].index.str.contains('p|q'))]]]
        #if analysis=='pq_mff':
        #    mutFreq = mutFreq.append(pd.DataFrame({'tumor':[tumor]*merged[tumor]['pq_mff'].shape[0], 'freq':list(merged[tumor]['pq_mff'].T.mean())}))
        tmp1 = [i.split('_') for i in merged[tumor][analysis].index]
        summaryTable[analysis].loc[tumor, 'Patients'] = merged[tumor][analysis].shape[1]
        for mutType in ['PAM','Fusion','CNAamp','Act','CNAdel','LoF']:
            summaryTable[analysis].loc[tumor, mutType] = len([i for i in tmp1 if i[1]==mutType])

print('Done.')

# %%
# Mutations with Delta_over_PAM greater than 0.35
# {tumor:oncomergeSummary[tumor]['pq_mff'][oncomergeSummary[tumor]['pq_mff']['Delta_over_PAM']>=0.3][['Symbol','Delta_over_PAM']] for tumor in oncomergeSummary.keys()}

# Get top mutaions by frequency
{tumor:oncomergeSummary[tumor]['pq_mff'][oncomergeSummary[tumor]['pq_mff']['Final_freq']>=0.5][['Symbol','Final_freq']] for tumor in oncomergeSummary.keys()}

"""
# Maximum final frequency with ties
removeThem = {}
for tumor in tumors:
    oncomergeSummary[tumor]['pq_mff'] = oncomergeSummary[tumor]['pq'].copy(deep=True)
    geneSets[tumor]['pq_mff'] = geneSets[tumor]['pq']
    tmp1 = oncomergeSummary[tumor]['pq'].loc[oncomergeSummary[tumor]['pq']['Genes_in_locus'].notna()]
    removeThem[tumor] = []
    for locus in list(set(tmp1.loc[tmp1['Genes_in_locus']>=10]['CNA_locus'].dropna())):
        tmp2 = tmp1.loc[tmp1['CNA_locus']==locus]
        maxFF = tmp2['Final_freq'].max()
        for gene1 in tmp2.loc[tmp2['Final_freq']<maxFF].index:
            for column1 in ['OM_type_selected', 'OM_empirical_p_value', 'OM_empirical_q_value', 'Genes_in_locus', 'Final_mutation_type', 'Final_freq', 'Delta_over_PAM']:
                oncomergeSummary[tumor]['pq_mff'].loc[gene1,column1] = np.nan
                removeThem[tumor].append(gene1)
        tmp3 = tmp2.loc[tmp2['Final_freq']==maxFF].index
        for gene1 in tmp3:
            oncomergeSummary[tumor]['pq_mff'].loc[gene1,'Genes_in_locus'] = len(tmp3)

for tumor in tumors:
    geneSets[tumor]['pq_mff'] = [i for i in geneSets[tumor]['pq_mff'] if not i[0] in removeThem[tumor]]
"""

#analyses = ['no_filter','pq','mff','pq_mff','mpf','pq_mpf']

# Genes per loci
#unfiltered_loci = [j for i in oncomergeSummary for j in list(oncomergeSummary[i]['no_filter'][['CNA_locus','Genes_in_locus']].dropna().drop_duplicates()['Genes_in_locus'])]+[j for i in oncomergeSummary for i in [1.0]*sum(oncomergeSummary[i]['no_filter']['Final_mutation_type']=='PAM')]

#filtered_loci = [j for i in oncomergeSummary for j in list(oncomergeSummary[i]['pq_mcr'][['CNA_locus','Genes_in_locus']].dropna().drop_duplicates()['Genes_in_locus'])]+[j for i in oncomergeSummary for i in [1.0]*sum(oncomergeSummary[i]['no_filter']['Final_mutation_type']=='PAM')]

print('Loci counts...')
loci_counts = []
for tumor in tumors:
    loci_counts.append(oncomergeSummary[tumor]['pq_mff'][['CNA_locus','Final_mutation_type']].dropna().groupby('CNA_locus').count()
        #oncomergeSummary[tumor]['mpf'][['CNA_locus','Final_mutation_type']].dropna().groupby('CNA_locus').count(),
        #oncomergeSummary[tumor]['pq_mpf'][['CNA_locus','Final_mutation_type']].dropna().groupby('CNA_locus').count(),
        #oncomergeSummary[tumor]['pq_mff_mpf'][['CNA_locus','Final_mutation_type']].dropna().groupby('CNA_locus').count()],
        )

forPlotting = pd.concat(loci_counts,axis=0)
forPlotting.columns = ['PQ & MFF'] #,'MPF','PQ & MPF', 'PQ & MFF & MPF']
loci_counts_long = []
for tumor in tumors:
        tmp5 = oncomergeSummary[tumor]['pq_mff'][['CNA_locus','Final_mutation_type']].dropna().groupby('CNA_locus').count()
        tmp5['filter'] = 'PQ & MFF'
        '''
        tmp5 = oncomergeSummary[tumor]['mpf'][['CNA_locus','Final_mutation_type']].dropna().groupby('CNA_locus').count()
        tmp5['filter'] = 'Minimum PAM Freq (MPF)'
        tmp6 = oncomergeSummary[tumor]['pq_mpf'][['CNA_locus','Final_mutation_type']].dropna().groupby('CNA_locus').count()
        tmp6['filter'] = 'PQ & MPF'
        tmp7 = oncomergeSummary[tumor]['pq_mff_mpf'][['CNA_locus','Final_mutation_type']].dropna().groupby('CNA_locus').count()
        tmp7['filter'] = 'PQ & MFF & MPF'
        '''
        tmp8 = pd.concat([tmp5], axis=0)
        tmp8['tumor'] = tumor
        loci_counts_long.append(tmp8)

forPlotting_violin = pd.concat(loci_counts_long,axis=0)
print('Done.')

'''
ax = forPlotting.plot.scatter("PQ & MCR", "No filter", color=[1,0,0,0.25])
ax.set_yscale('log')
ax.set_xscale('log')
diag_line, = ax.plot((0,300), (0,300), ls="--", c=".3")
#plt.ylim(-100,250)
#plt.xlim(-100,250)

#ax = forPlotting.plot.density()
#ax = np.log10(forPlotting[['no_filter','pq_mcr']]).plot.hist(density=True, alpha=0.5)
ax = forPlotting[['No filter','PQ & MCR']].plot.hist(density=True, alpha=0.5, bins=50)
#ax.set_yscale('log')
#ax.set_xscale('log')
#diag_line, = ax.plot((0,300), (0,300), ls="--", c=".3")
#plt.ylim(-100,250)
#plt.xlim(-100,250)
plt.show()

#ax = forPlotting.plot.density()
#ax = np.log10(forPlotting[['no_filter','pq_mcr']]).plot.hist(density=True, alpha=0.5)
ax = forPlotting[['No filter','PQ & MCR']].plot.hist(density=True, alpha=0.5, bins=50)
#ax.set_yscale('log')
#ax.set_xscale('log')
#diag_line, = ax.plot((0,300), (0,300), ls="--", c=".3")
#plt.ylim(-100,250)
#plt.xlim(-100,250)
plt.show()
'''
# %%
print('Effect of filter plot...')
with PdfPages('effectOfFiltersOnGenesPerLoci_8_10_2021.pdf') as pp:
    fig, ax = plt.subplots(2, 1, sharex=True, figsize=(4,4))
    plt.xlabel('Genes per locus')
    fig.text(0, 0.5, 'Frequency', va='center', rotation='vertical')

    ax[0].set_ylim(0.125, 0.20)  # outliers only
    ax[1].set_ylim(-0.0025, .025)  # most of the data
    ax[0].set_ylabel("")  # outliers only
    ax[1].set_ylabel("")  # most of the data

    ax[0].spines['bottom'].set_visible(False)
    ax[1].spines['top'].set_visible(False)
    ax[0].xaxis.tick_top()
    ax[0].tick_params(labeltop='off')  # don't put tick labels at the top
    ax[1].xaxis.tick_bottom()

    d = .015  # how big to make the diagonal lines in axes coordinates
    kwargs = dict(transform=ax[0].transAxes, color='k', clip_on=False)
    ax[0].plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
    ax[0].plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal
    kwargs.update(transform=ax[1].transAxes)  # switch to the bottom axes
    ax[1].plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
    ax[1].plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal
    fig.savefig(pp, format='pdf')

    '''fig, ax = plt.subplots(figsize=(4,4))
    forPlotting.boxplot(ax=ax)
    ax.set_yscale('log')
    fig.savefig(pp, format='pdf')
    '''

    fig, ax = plt.subplots(figsize=(4,4))
    sns.boxplot(x='filter', y='Final_mutation_type', data=forPlotting_violin, ax=ax)
    #ax.set_ylim(1,500)
    ax.set_yscale('log')
    ax.set_ylabel('Genes per locus')
    ax.set_xlabel('')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45)
    fig.tight_layout()
    fig.savefig(pp, format='pdf')

    fig, ax = plt.subplots(figsize=(4,4))
    sns.boxenplot(x='filter', y='Final_mutation_type', data=forPlotting_violin, ax=ax)
    #sns.swarmplot(x='filter', y='Final_mutation_type', data=forPlotting_violin, ax=ax, size=2, color=[0.5,0.5,0.5,0.2])
    #ax.set_ylim(1,500)
    ax.set_yscale('log')
    ax.set_ylabel('Genes per locus')
    ax.set_xlabel('')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    fig.tight_layout()
    fig.savefig(pp, format='pdf')

    fig, ax = plt.subplots(figsize=(4,4))
    plt.bar(x=np.arange(5), height=[6193,5463,1489,1399,1085])
    #sns.swarmplot(x='filter', y='Final_mutation_type', data=forPlotting_violin, ax=ax, size=2, color=[0.5,0.5,0.5,0.2])
    #ax.set_ylim(1,500)
    #ax.set_yscale('log')
    ax.set_ylabel('Somatically mutated genes')
    ax.set_xlabel('')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    fig.tight_layout()
    fig.savefig(pp, format='pdf')

    fig, ax = plt.subplots(figsize=(4,4))
    sns.violinplot(x='filter', y='Final_mutation_type', data=forPlotting_violin, ax=ax)
    #sns.swarmplot(x='filter', y='Final_mutation_type', data=forPlotting_violin, ax=ax, size=1, color=[0.5,0.5,0.5,0.2])
    #ax.set_ylim(1,500)
    ax.set_yscale('log')
    ax.set_ylabel('Genes per locus')
    ax.set_xlabel('')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    fig.tight_layout()
    fig.savefig(pp, format='pdf')

print ('Done.')

# %%

print('Summary table...')
# Write out summary table
summaryTable['pq_mff'].to_csv('summaryTable_pq_mff_8_10_2021.csv')
#summaryTable['mpf'].to_csv('summaryTable_pq_mpf_8_10_2021.csv')
#summaryTable['pq_mpf'].to_csv('summaryTable_pq_mpf_8_10_2021.csv')
#summaryTable['pq_mff_mpf'].to_csv('summaryTable_pq_mpf_8_10_2021.csv')

with PdfPages('summaryPlots_pq_mff.pdf') as pp:
    fig, ax = plt.subplots(2,1,figsize=(8,6)) #sharex='col'
    # Plot final mutation counts
    #fig, ax = plt.subplots(figsize=(8,4))
    tmp = copy.copy(summaryTable['pq_mff'])
    del tmp['Patients']
    #summaryTable['pq_mff'].T.reindex(summaryTable['pq_mff'].T.sum().sort_values().index, axis=1).T.plot(kind='bar', stacked=True, ax=ax[0], color=[colors1[i] for i in summaryTable['pq_mff'].columns if not i=='Patients'])
    tmp.T.reindex(tmp.T.sum().sort_values().index, axis=1).T.plot(kind='bar', stacked=True, ax=ax[0], color=[colors1[i] for i in tmp.columns])
    ax[0].set_ylabel('Number of somatic mutations')
    #fig.tight_layout()
    #fig.savefig(pp, format='pdf')

    # Plot final mutation frequencies by tumor type
    #mutFreq = pd.DataFrame(column=['tumor','freq'])

    o1 = tmp.T.sum().sort_values().index
    #fig, ax = plt.subplots(figsize=(8,4))
    #sns.stripplot(x='tumor',y='freq', data=mutFreq, ax=ax, order=o1, size=3, jitter=True)
    sns.swarmplot(x='tumor',y='freq', data=mutFreq, ax=ax[1], order=o1, size=2, edgecolor='gray',linewidth=0.01, dodge=True) #, palette= sns.color_palette('deep', 33))
    ax[1].set_xticklabels(ax[1].get_xticklabels(),rotation=90)
    ax[1].set_ylabel('Somatic mutation frequency')
    ax[1].set_xlabel('Tumor type')
    ax[1].set(ylim=(0,1))
    fig.tight_layout()
    fig.savefig(pp, format='pdf')

    '''# Plot total mutation counts versus number of patients
    fig, ax = plt.subplots(figsize=(4,4))
    comparisonTable = pd.DataFrame(index=tumors, columns=['mut_counts','num_patients'])
    comparisonTable['mut_counts'] = tmp.T.sum().sort_values()
    for tumor in tumors:
        comparisonTable.loc[tumor,'num_patients'] = len(merged[tumor]['pq_mff_mhc'].columns)
    ax = sns.regplot(x='mut_counts', y='num_patients', data = comparisonTable)
    ax.set_ylabel('Number of patients in study')
    ax.set_xlabel('Number of somatic mutations')
    pc1 = pearsonr(comparisonTable['mut_counts'],comparisonTable['num_patients'])
    ax.set_title('R = '+'%.2f'%pc1[0]+'; p-value = '+'%.3f'%pc1[1])
    for tumor in tumors:
        ax.text(comparisonTable.loc[tumor,'mut_counts'], comparisonTable.loc[tumor,'num_patients'], tumor, horizontalalignment='center', fontsize=6)
    fig.tight_layout()
    fig.savefig(pp, format='pdf')
    '''
print('Done.')
# %%
####################
## Hypergeometric ##
####################
"""
print('Hypergeometric enrichment...')
# Compute hypergeometric for cancer-specific
#analyses = ['no_filter','pq','pq_mff']
csp_hypgeoResults = []
for analysis in analyses+['mhc','pq_mff_mhc']:
    print(analysis)
    gs_cancerSpecific_Act = [[i,j[0]] for i in geneSets for j in geneSets[i][analysis] if j[1]=='Act' or j[1]=='CNAamp']
    gs_cancerSpecific_LoF = [[i,j[0]] for i in geneSets for j in geneSets[i][analysis] if j[1]=='LoF' or j[1]=='CNAdel']
    gs_Act = list(set([j[0] for i in geneSets for j in geneSets[i][analysis] if j[1]=='Act' or j[1]=='CNAamp']))
    gs_LoF = list(set([j[0] for i in geneSets for j in geneSets[i][analysis] if j[1]=='LoF' or j[1]=='CNAdel']))

    # TCGA_Consensus (Activating)
    eSet_cancerSpecific_Act = [[i,j] for i in enrichmentSets['TCGA_Consensus'] for j in enrichmentSets['TCGA_Consensus'][i]['Activating'] if not i in ['PANCAN','Aggregate']]
    ovlp_cancerSpecific_Act = [i for i in gs_cancerSpecific_Act if i in eSet_cancerSpecific_Act]
    bkgd_cancerSpecific_Act = [[i,j] for i in background for j in background[i]['Activating']]

    # TCGA_Consensus Loss of function
    eSet_cancerSpecific_LoF = [[i,j] for i in enrichmentSets['TCGA_Consensus'] for j in enrichmentSets['TCGA_Consensus'][i]['LossOfFunction'] if not i in ['PANCAN','Aggregate']]
    ovlp_cancerSpecific_LoF = [i for i in gs_cancerSpecific_LoF if i in eSet_cancerSpecific_LoF]
    bkgd_cancerSpecific_LoF = [[i,j] for i in background for j in background[i]['LossOfFunction']]

    # CGC (Activating)
    eSet_cancerSpecific_cgcAct = list(set([j for i in enrichmentSets['CGC'] for j in enrichmentSets['CGC'][i]['Activating'] if not i in ['PANCAN','Aggregate']]))
    ovlp_cancerSpecific_cgcAct = [i for i in gs_Act if i in eSet_cancerSpecific_cgcAct]
    bkgd_cancerSpecific_cgcAct = list(set([j for i in background for j in background[i]['Activating']]))

    # CGC Loss of function
    eSet_cancerSpecific_cgcLoF = list(set([j for i in enrichmentSets['CGC'] for j in enrichmentSets['CGC'][i]['LossOfFunction'] if not i in ['PANCAN','Aggregate']]))
    ovlp_cancerSpecific_cgcLoF = [i for i in gs_LoF if i in eSet_cancerSpecific_cgcLoF]
    bkgd_cancerSpecific_cgcLoF = list(set([j for i in background for j in background[i]['LossOfFunction']]))

    # Vogelstein (Activating)
    eSet_cancerSpecific_volAct = list(set([j for i in enrichmentSets['Vogelstein'] for j in enrichmentSets['Vogelstein'][i]['Activating'] if not i in ['PANCAN','Aggregate']]))
    ovlp_cancerSpecific_volAct = [i for i in gs_Act if i in eSet_cancerSpecific_volAct]
    bkgd_cancerSpecific_volAct = list(set([j for i in background for j in background[i]['Activating']]))

    # Vogelstein Loss of function
    eSet_cancerSpecific_volLoF = list(set([j for i in enrichmentSets['Vogelstein'] for j in enrichmentSets['Vogelstein'][i]['LossOfFunction'] if not i in ['PANCAN','Aggregate']]))
    ovlp_cancerSpecific_volLoF = [i for i in gs_LoF if i in eSet_cancerSpecific_volLoF]
    bkgd_cancerSpecific_volLoF = list(set([j for i in background for j in background[i]['LossOfFunction']]))

    # OncoROLE (Activating)
    eSet_cancerSpecific_onrAct = list(set([j for i in enrichmentSets['OncoROLE'] for j in enrichmentSets['OncoROLE'][i]['Activating'] if not i in ['PANCAN','Aggregate']]))
    ovlp_cancerSpecific_onrAct = [i for i in gs_Act if i in eSet_cancerSpecific_onrAct]
    bkgd_cancerSpecific_onrAct = list(set([j for i in background for j in background[i]['Activating']]))

    # OncoROLE Loss of function
    eSet_cancerSpecific_onrLoF = list(set([j for i in enrichmentSets['OncoROLE'] for j in enrichmentSets['OncoROLE'][i]['LossOfFunction'] if not i in ['PANCAN','Aggregate']]))
    ovlp_cancerSpecific_onrLoF = [i for i in gs_LoF if i in eSet_cancerSpecific_onrLoF]
    bkgd_cancerSpecific_onrLoF = list(set([j for i in background for j in background[i]['LossOfFunction']]))

    # Tokheim (Activating)
    eSet_cancerSpecific_tokAct = list(set([j for i in enrichmentSets['Tokheim'] for j in enrichmentSets['Tokheim'][i]['Activating'] if not i in ['PANCAN','Aggregate']]))
    ovlp_cancerSpecific_tokAct = [i for i in gs_Act if i in eSet_cancerSpecific_tokAct]
    bkgd_cancerSpecific_tokAct = list(set([j for i in background for j in background[i]['Activating']]))

    # Tokheim Loss of function
    eSet_cancerSpecific_tokLoF = list(set([j for i in enrichmentSets['Tokheim'] for j in enrichmentSets['Tokheim'][i]['LossOfFunction'] if not i in ['PANCAN','Aggregate']]))
    ovlp_cancerSpecific_tokLoF = [i for i in gs_LoF if i in eSet_cancerSpecific_tokLoF]
    bkgd_cancerSpecific_tokLoF = list(set([j for i in background for j in background[i]['LossOfFunction']]))

    '''# Moonlight (Activating)
    eSet_cancerSpecific_mlAct = list(set([j for i in enrichmentSets['Moonlight'] for j in enrichmentSets['Moonlight'][i]['Activating'] if not i in ['PANCAN','Aggregate']]))
    ovlp_cancerSpecific_mlAct = [i for i in gs_Act if i in eSet_cancerSpecific_mlAct]
    bkgd_cancerSpecific_mlAct = list(set([j for i in background for j in background[i]['Activating']]))

    # Moonlight Loss of function
    eSet_cancerSpecific_mlLoF = list(set([j for i in enrichmentSets['Moonlight'] for j in enrichmentSets['Moonlight'][i]['LossOfFunction'] if not i in ['PANCAN','Aggregate']]))
    ovlp_cancerSpecific_mlLoF = [i for i in gs_LoF if i in eSet_cancerSpecific_mlLoF]
    bkgd_cancerSpecific_mlLoF = list(set([j for i in background for j in background[i]['LossOfFunction']]))
    '''

    # Oncogene (http://ongene.bioinfo-minzhao.org/ongene_human.txt)
    eSet_cancerSpecific_og = list(set([i for i in enrichmentSets['ONGene']['Any']['Activating']]))
    ovlp_cancerSpecific_og = [i for i in gs_Act if i in eSet_cancerSpecific_og]
    bkgd_cancerSpecific_og = list(set([j for i in background for j in background[i]['Activating']]))

    # Tumor suppressors (https://bioinfo.uth.edu/TSGene/Human_TSGs.txt)
    eSet_cancerSpecific_ts = list(set([i for i in enrichmentSets['TSGene']['Any']['LossOfFunction']]))
    ovlp_cancerSpecific_ts = [i for i in gs_LoF if i in eSet_cancerSpecific_ts]
    bkgd_cancerSpecific_ts = list(set([j for i in background for j in background[i]['LossOfFunction']]))

    # All = CGC + TCGA_Consensus + Vogelstein
    gs_cancerSpecific_All = list(set([j[0] for i in geneSets for j in geneSets[i][analysis]]))
    eSet_cancerSpecific_All = list(set(enrichmentSets['CGC']['Any']['All']+enrichmentSets['TCGA_Consensus']['Any']['All']+enrichmentSets['Vogelstein']['Any']['All']))
    ovlp_cancerSpecific_All = [i for i in gs_cancerSpecific_All if i in eSet_cancerSpecific_All]
    bkgd_cancerSpecific_All = list(set([j for i in background for j in background[i]['Aggregate']]))

    ## Compute hypergeo
    # ONGene
    csp_hypgeoResults.append([analysis,'ONGene','Activating']+enrichment(ovlp_cancerSpecific_og, bkgd_cancerSpecific_og, gs_Act, eSet_cancerSpecific_og)+[' '.join([str(i) for i in ovlp_cancerSpecific_og])])
    csp_hypgeoResults.append([analysis,'ONGene','LossOfFunction']+enrichment([i for i in gs_LoF if i in eSet_cancerSpecific_og], bkgd_cancerSpecific_og, gs_LoF, eSet_cancerSpecific_og)+[' '.join([str(i) for i in ovlp_cancerSpecific_ts])])

    # TCGA Consensus
    csp_hypgeoResults.append([analysis,'TCGA_Consensus (Act)','Activating']+enrichment(ovlp_cancerSpecific_Act, bkgd_cancerSpecific_Act, gs_cancerSpecific_Act, eSet_cancerSpecific_Act)+[' '.join([';'.join([str(j) for j in i]) for i in ovlp_cancerSpecific_Act])])
    csp_hypgeoResults.append([analysis,'TCGA_Consensus (Act)','LossOfFunction']+enrichment([k for k in gs_cancerSpecific_LoF if k in eSet_cancerSpecific_Act], bkgd_cancerSpecific_Act, gs_cancerSpecific_LoF, eSet_cancerSpecific_Act)+[' '.join([';'.join([str(j) for j in i]) for i in [k for k in gs_cancerSpecific_LoF if k in eSet_cancerSpecific_Act]])])
    csp_hypgeoResults.append([analysis,'TCGA_Consensus (LoF)','Activating']+enrichment([k for k in gs_cancerSpecific_Act if k in eSet_cancerSpecific_LoF], bkgd_cancerSpecific_LoF, gs_cancerSpecific_Act, eSet_cancerSpecific_LoF)+[' '.join([';'.join([str(j) for j in i]) for i in [i for i in gs_cancerSpecific_Act if i in [k for k in gs_cancerSpecific_Act if k in eSet_cancerSpecific_LoF]]])])
    csp_hypgeoResults.append([analysis,'TCGA_Consensus (LoF)','LossOfFunction']+enrichment(ovlp_cancerSpecific_LoF, bkgd_cancerSpecific_LoF, gs_cancerSpecific_LoF, eSet_cancerSpecific_LoF)+[' '.join([';'.join([str(j) for j in i]) for i in ovlp_cancerSpecific_LoF])])

    # CGC
    csp_hypgeoResults.append([analysis,'CGC (Act)','Activating']+enrichment(ovlp_cancerSpecific_cgcAct, bkgd_cancerSpecific_cgcAct, gs_Act, eSet_cancerSpecific_cgcAct)+[' '.join([str(i) for i in ovlp_cancerSpecific_cgcAct])])
    csp_hypgeoResults.append([analysis,'CGC (Act)','LossOfFunction']+enrichment(list(set(gs_LoF).intersection(eSet_cancerSpecific_cgcAct)), bkgd_cancerSpecific_cgcAct, gs_LoF, eSet_cancerSpecific_cgcAct)+[' '.join([str(i) for i in list(set(gs_LoF).intersection(eSet_cancerSpecific_cgcAct))])])
    csp_hypgeoResults.append([analysis,'CGC (LoF)','Activating']+enrichment(list(set(gs_Act).intersection(eSet_cancerSpecific_cgcLoF)), bkgd_cancerSpecific_cgcLoF, gs_Act, eSet_cancerSpecific_cgcLoF)+[' '.join([str(i) for i in gs_Act if i in list(set(gs_Act).intersection(eSet_cancerSpecific_cgcLoF))])])
    csp_hypgeoResults.append([analysis,'CGC (LoF)','LossOfFunction']+enrichment(ovlp_cancerSpecific_cgcLoF, bkgd_cancerSpecific_cgcLoF, gs_LoF, eSet_cancerSpecific_cgcLoF)+[' '.join([str(i) for i in ovlp_cancerSpecific_cgcLoF])])

    # Volgelstein
    csp_hypgeoResults.append([analysis,'Vogelstein (Act)','Activating']+enrichment(ovlp_cancerSpecific_volAct, bkgd_cancerSpecific_volAct, gs_Act, eSet_cancerSpecific_volAct)+[' '.join([str(i) for i in ovlp_cancerSpecific_volAct])])
    csp_hypgeoResults.append([analysis,'Vogelstein (Act)','LossOfFunction']+enrichment(list(set(gs_LoF).intersection(eSet_cancerSpecific_volAct)), bkgd_cancerSpecific_volAct, gs_LoF, eSet_cancerSpecific_volAct)+[' '.join([str(i) for i in list(set(gs_LoF).intersection(eSet_cancerSpecific_volAct))])])
    csp_hypgeoResults.append([analysis,'Vogelstein (LoF)','Activating']+enrichment(list(set(gs_Act).intersection(eSet_cancerSpecific_volLoF)), bkgd_cancerSpecific_volLoF, gs_Act, eSet_cancerSpecific_volLoF)+[' '.join([str(i) for i in gs_Act if i in list(set(gs_Act).intersection(eSet_cancerSpecific_volLoF))])])
    csp_hypgeoResults.append([analysis,'Vogelstein (LoF)','LossOfFunction']+enrichment(ovlp_cancerSpecific_volLoF, bkgd_cancerSpecific_volLoF, gs_LoF, eSet_cancerSpecific_volLoF)+[' '.join([str(i) for i in ovlp_cancerSpecific_volLoF])])

    # OncoROLE
    csp_hypgeoResults.append([analysis,'OncoROLE (Act)','Activating']+enrichment(ovlp_cancerSpecific_onrAct, bkgd_cancerSpecific_onrAct, gs_Act, eSet_cancerSpecific_onrAct)+[' '.join([str(i) for i in ovlp_cancerSpecific_onrAct])])
    csp_hypgeoResults.append([analysis,'OncoROLE (Act)','LossOfFunction']+enrichment(list(set(gs_LoF).intersection(eSet_cancerSpecific_onrAct)), bkgd_cancerSpecific_onrAct, gs_LoF, eSet_cancerSpecific_onrAct)+[' '.join([str(i) for i in list(set(gs_LoF).intersection(eSet_cancerSpecific_onrAct))])])
    csp_hypgeoResults.append([analysis,'OncoROLE (LoF)','Activating']+enrichment(list(set(gs_Act).intersection(eSet_cancerSpecific_onrLoF)), bkgd_cancerSpecific_onrLoF, gs_Act, eSet_cancerSpecific_onrLoF)+[' '.join([str(i) for i in gs_Act if i in list(set(gs_Act).intersection(eSet_cancerSpecific_onrLoF))])])
    csp_hypgeoResults.append([analysis,'OncoROLE (LoF)','LossOfFunction']+enrichment(ovlp_cancerSpecific_onrLoF, bkgd_cancerSpecific_onrLoF, gs_LoF, eSet_cancerSpecific_onrLoF)+[' '.join([str(i) for i in ovlp_cancerSpecific_onrLoF])])

    # Tokheim
    csp_hypgeoResults.append([analysis,'Tokheim (Act)','Activating']+enrichment(ovlp_cancerSpecific_tokAct, bkgd_cancerSpecific_tokAct, gs_Act, eSet_cancerSpecific_tokAct)+[' '.join([str(i) for i in ovlp_cancerSpecific_tokAct])])
    csp_hypgeoResults.append([analysis,'Tokheim (Act)','LossOfFunction']+enrichment(list(set(gs_LoF).intersection(eSet_cancerSpecific_tokAct)), bkgd_cancerSpecific_tokAct, gs_LoF, eSet_cancerSpecific_tokAct)+[' '.join([str(i) for i in list(set(gs_LoF).intersection(eSet_cancerSpecific_tokAct))])])
    csp_hypgeoResults.append([analysis,'Tokheim (LoF)','Activating']+enrichment(list(set(gs_Act).intersection(eSet_cancerSpecific_tokLoF)), bkgd_cancerSpecific_tokLoF, gs_Act, eSet_cancerSpecific_tokLoF)+[' '.join([str(i) for i in gs_Act if i in list(set(gs_Act).intersection(eSet_cancerSpecific_tokLoF))])])
    csp_hypgeoResults.append([analysis,'Tokheim (LoF)','LossOfFunction']+enrichment(ovlp_cancerSpecific_tokLoF, bkgd_cancerSpecific_tokLoF, gs_LoF, eSet_cancerSpecific_tokLoF)+[' '.join([str(i) for i in ovlp_cancerSpecific_tokLoF])])

    '''# Moonlight
    csp_hypgeoResults.append([analysis,'Moonlight (Act)','Activating']+enrichment(ovlp_cancerSpecific_mlAct, bkgd_cancerSpecific_mlAct, gs_Act, eSet_cancerSpecific_mlAct)+[' '.join([str(i) for i in ovlp_cancerSpecific_mlAct])])
    csp_hypgeoResults.append([analysis,'Moonlight (Act)','LossOfFunction']+enrichment(list(set(gs_LoF).intersection(eSet_cancerSpecific_mlAct)), bkgd_cancerSpecific_mlAct, gs_LoF, eSet_cancerSpecific_mlAct)+[' '.join([str(i) for i in list(set(gs_LoF).intersection(eSet_cancerSpecific_mlAct))])])
    csp_hypgeoResults.append([analysis,'Moonlight (LoF)','Activating']+enrichment(list(set(gs_Act).intersection(eSet_cancerSpecific_mlLoF)), bkgd_cancerSpecific_mlLoF, gs_Act, eSet_cancerSpecific_mlLoF)+[' '.join([str(i) for i in gs_Act if i in list(set(gs_Act).intersection(eSet_cancerSpecific_mlLoF))])])
    csp_hypgeoResults.append([analysis,'Moonlight (LoF)','LossOfFunction']+enrichment(ovlp_cancerSpecific_mlLoF, bkgd_cancerSpecific_mlLoF, gs_LoF, eSet_cancerSpecific_mlLoF)+[' '.join([str(i) for i in ovlp_cancerSpecific_mlLoF])])
    '''

    # All
    csp_hypgeoResults.append([analysis,'All','All']+enrichment(ovlp_cancerSpecific_All, bkgd_cancerSpecific_All, gs_cancerSpecific_All, eSet_cancerSpecific_All)+[' '.join([str(i) for i in ovlp_cancerSpecific_All])])

with open('csp_hypgeoResults_8_10_2021.csv','w') as outFile:
    outFile.write('Analysis,Enrichment_Source,OncoMerge_Direction,Overlap,OncoMerge_Discovered,Enrichment_Set,Background,P-Value,Combos\n')
    outFile.write('\n'.join([','.join([str(j) for j in i]) for i in csp_hypgeoResults]))

print('Done.')
"""
#%%
######################
## Enrichment plots ##
######################

print('Enrichment plots...')

csp = pd.read_csv('csp_hypgeoResults_8_10_2021.csv',header=0)
tmp1 = csp.pivot(index=['Enrichment_Source','OncoMerge_Direction'],columns='Analysis', values='P-Value')
nlpv = np.log10(tmp1).replace(-np.inf,-16)

with PdfPages('enrichment_negLogPV_8_10_2021_heatmap.pdf') as pp:
     fig, ax = plt.subplots(1,1,figsize=(3,4))
     sns.set(font_scale=0.4)
     sns.heatmap(nlpv[['pq_mff']].loc[[('All','All'), ('TCGA_Consensus (Act)','Activating'), ('TCGA_Consensus (Act)', 'LossOfFunction'), ('TCGA_Consensus (LoF)', 'Activating'), ('TCGA_Consensus (LoF)', 'LossOfFunction'), ('CGC (Act)', 'Activating'), ('CGC (Act)', 'LossOfFunction'), ('CGC (LoF)', 'Activating'), ('CGC (LoF)', 'LossOfFunction'), ('Vogelstein (Act)', 'Activating'), ('Vogelstein (Act)', 'LossOfFunction'), ('Vogelstein (LoF)','Activating'),
         ('Vogelstein (LoF)', 'LossOfFunction'), ('OncoROLE (Act)','Activating'), ('OncoROLE (Act)', 'LossOfFunction'), ('OncoROLE (LoF)', 'Activating'), ('OncoROLE (LoF)', 'LossOfFunction'), ('Tokheim (Act)', 'Activating'), ('Tokheim (Act)', 'LossOfFunction'), ('Tokheim (LoF)','Activating'), ('Tokheim (LoF)', 'LossOfFunction')]], xticklabels=True, yticklabels=True, cmap = sns.color_palette('blend:#084594,#eff3ff', as_cmap=True))
     fig.tight_layout()
     fig.savefig(pp, format='pdf')


#%%
#################
## Value Added ##
#################

print('Value added...')

#mj = mg.MyGeneInfo() # pull in function to map genes
tmp0 = pd.read_csv("C:/Users/nihan/Dropbox/OncoMerge_pathways/TCGA/OncoMerge_input_g2e_converter.csv",index_col=0)
convertMe = tmp0['Locus ID']

pams = {}
cnvs = {}
pamTumors = {}
allThreshs = {}
for tumor in tumors:
    print('\t'+tumor)

    ## PAMs
    print('\t\tPAM')
    pams[tumor] = {}
    # Need to characterize PAMs prior to oncomerging
    #  - Somatic mutation freqeuncy >= 5%
    #  - MutSigCV2 q-value <= 0.05
    mutSigCV2 = pd.read_csv('C:/Users/nihan/Dropbox/OncoMerge_pathways/TCGA/MutSig2cv/'+tumor+'_sig2cv.csv', header=0, index_col=0)
    genes = list(convertMe[set(convertMe.index).intersection(list(mutSigCV2.loc[mutSigCV2['q']<=0.05,'gene']))])
    #genes = list(mj.querymany(list(mutSigCV2.loc[mutSigCV2['q']<=0.05,'gene']), scopes='symbol', fields='entrezgene', species='human', as_dataframe=True, df_index=True, entrezonly=True, verbose=False)['entrezgene'].dropna())
    pamTumor = pd.read_csv('C:/Users/nihan/Dropbox/OncoMerge_pathways/TCGA/PAM/'+tumor+'_somMutMC3.csv', header=0, index_col=0)
    pamTumors[tumor] = pamTumor
    tmp1 = set(pamTumor.index[pamTumor.T.mean()>=0.05]).intersection([int(i) for i in genes])
    tmp1_means = list(pamTumor.T.mean()[list(tmp1)])
    pams[tumor] = zip(tmp1, tmp1_means)

    ## CNVs
    print('\t\tCNV')
    cnvs[tumor] = {'amp':{},'del':{}}
    # Get list of significantly CNA amplified genes
    ampLoci = {}
    amp1 = pd.read_csv('C:/Users/nihan/Dropbox/OncoMerge_pathways/TCGA/GISTIC/'+tumor+'/amp_genes.conf_99.txt',index_col=0,sep='\t')
    for col1 in amp1.columns:
        if float(amp1[col1]['residual q value'])<=0.05 and not (col1[0]=='X' or col1=='Y'):
            #tmpA = mj.querymany([i.lstrip('[').rstrip(']').split('|')[0] for i in list(amp1[col1].dropna()[3:])], scopes='symbol', fields='entrezgene', species='human', as_dataframe=True, df_index=True, entrezonly=True, verbose=False)
            tmpA = [convertMe[i] for i in [i.lstrip('[').rstrip(']').split('|')[0] for i in list(amp1[col1].dropna()[3:])] if i in convertMe.index]
            if len(tmpA)>1:
                ampLoci[col1] = [int(i) for i in tmpA]

    # Get list of significantly CNA deleted genes
    print('\t\tGISTIC')
    delLoci = {}
    del1 = pd.read_csv('C:/Users/nihan/Dropbox/OncoMerge_pathways/TCGA/GISTIC/'+tumor+'/del_genes.conf_99.txt',index_col=0,sep='\t')
    for col1 in del1.columns:
        if float(del1[col1]['residual q value'])<=0.05 and not (col1[0]=='X' or col1=='Y'):
            tmpD = [convertMe[i] for i in [i.lstrip('[').rstrip(']').split('|')[0] for i in list(del1[col1].dropna()[3:])] if i in convertMe.index]
            #tmpD = mj.querymany([i.lstrip('[').rstrip(']').split('|')[0] for i in list(del1[col1].dropna()[3:])], scopes='symbol', fields='entrezgene', species='human', as_dataframe=True, df_index=True, entrezonly=True, verbose=False)
            if len(tmpD)>1:
                delLoci[col1] = [int(i) for i in tmpD]

    # Load CNV thresholded gene data
    allThresh = pd.read_csv('C:/Users/nihan/Dropbox/OncoMerge_pathways/TCGA/GISTIC/'+tumor+'/all_thresholded.by_genes.txt', header=0, index_col=1, sep='\t').drop(columns=['Gene Symbol','Cytoband'])
    allThreshs[tumor] = allThresh
    # Get amplified loci freq >=0.05
    tmpAmp = list(set([gene for loci in ampLoci for gene in ampLoci[loci] if gene in allThresh.index and float(sum(allThresh.loc[gene]>=2))/allThresh.shape[1]>=0.05]))
    tmpAmp_mean = [float(sum(allThresh.loc[gene]>=2))/allThresh.shape[1] for gene in tmpAmp]
    cnvs[tumor]['amp'] = dict(zip(tmpAmp, tmpAmp_mean))
    tmpDel = list(set([gene for loci in delLoci for gene in delLoci[loci] if gene in allThresh.index and float(sum(allThresh.loc[gene]<=-2))/allThresh.shape[1]>=0.05]))
    tmpDel_mean = [float(sum(allThresh.loc[gene]<=-2))/allThresh.shape[1] for gene in tmpDel]
    cnvs[tumor]['del'] = dict(zip(tmpDel, tmpDel_mean))

notPams = dict([[tumor, [i[0] for i in geneSets[tumor]['pq_mff'] if not i[0] in [j[0] for j in pams[tumor]]]] for tumor in tumors])
valueAdded = dict([[tumor, [i[0] for i in geneSets[tumor]['pq_mff'] if not ((i[0] in [j[0] for j in pams[tumor]]) or (i[0] in [m for m in [k for k in cnvs[tumor]['amp']]+[l for l in cnvs[tumor]['del']]]))]] for tumor in tumors])
pd.DataFrame([[i,str(len(valueAdded[i])),' '.join([str(j) for j in valueAdded[i]])] for i in valueAdded]).to_csv('valueAdded_TCGA_8_10_2021.csv')
valueAddedPlot = pd.DataFrame([[tumor,i] for tumor in tumors for i in valueAdded[tumor]])
valueAddedPlot.columns = ['Tumor','Genes Added']
valueAddedPlot.index = valueAddedPlot['Tumor']

deltaOverPAM = pd.DataFrame([[i,tumor] for tumor in tumors for i in oncomergeSummary[tumor]['pq_mff']['Delta_over_PAM'].dropna()])
deltaOverPAM.columns = ['Delta','Tumor']
deltaOverPAM.index = deltaOverPAM['Tumor']

with PdfPages('valueAddedPlots_pq_mff.pdf') as pp:
    fig, ax = plt.subplots(2,1,figsize=(8,6)) #sharex='col'
    # Plot final mutation counts
    #fig, ax = plt.subplots(figsize=(8,4))
    sns.countplot(x='Tumor', data=valueAddedPlot.loc[o1], ax=ax[0], order=o1) #, dodge=True)
    ax[0].set_xticklabels(ax[0].get_xticklabels(),rotation=90)
    ax[0].set_ylabel('Genes Added')
    # Plot final mutation frequencies by tumor type
    sns.swarmplot(x='Tumor',y='Delta', data=deltaOverPAM, ax=ax[1], order=o1, size=2, edgecolor='gray',linewidth=0.01, dodge=True) #, palette= sns.color_palette('deep', 33))
    ax[1].set_xticklabels(ax[1].get_xticklabels(),rotation=90)
    ax[1].set_ylabel('Frequency added over PAM')
    ax[1].set_xlabel('Tumor type')
    ax[1].set(ylim=(0,1))
    fig.tight_layout()
    fig.savefig(pp, format='pdf')

with PdfPages('summaryAndValueAddedPlots_pq_mff.pdf') as pp:
    fig, ax = plt.subplots(4,1,figsize=(8,10)) #, sharex='all')
    # Plot final mutation counts
    #fig, ax = plt.subplots(figsize=(8,4))
    tmp = copy.copy(summaryTable['pq_mff'])
    del tmp['Patients']
    tmp.T.reindex(tmp.T.sum().sort_values().index, axis=1).T.plot(kind='bar', stacked=True, ax=ax[0], color=[colors1[i] for i in tmp.columns if not i=='Patients'])
    ax[0].set_ylabel('Number of somatic mutations')
    o1 = tmp.T.sum().sort_values().index
    sns.swarmplot(x='tumor',y='freq', data=mutFreq, ax=ax[1], order=o1, size=2, edgecolor='gray',linewidth=0.01, dodge=True)
    ax[1].set_xticklabels(ax[1].get_xticklabels(),rotation=90)
    ax[1].set_ylabel('Somatic mutation frequency')
    ax[1].set_xlabel(None)
    ax[1].set(ylim=(0,1))
    # Plot final mutation frequencies by tumor type
    sns.swarmplot(x='Tumor',y='Delta', data=deltaOverPAM, ax=ax[2], order=o1, size=2, edgecolor='gray',linewidth=0.01, dodge=True)
    ax[2].set_xticklabels(ax[0].get_xticklabels(),rotation=90)
    ax[2].set_ylabel('Frequency added over PAM')
    ax[2].set_xlabel(None)
    #ax[2].set_xlabel('Tumor type')
    ax[2].set(ylim=(0,1))
    # Plot final mutation counts
    sns.countplot(x='Tumor', data=valueAddedPlot.loc[o1], ax=ax[3], order=o1)
    ax[3].set_xticklabels(ax[0].get_xticklabels(),rotation=90)
    ax[3].set_ylabel('Genes added')
    ax[3].set_xlabel('Tumor type')
    plt.setp(ax[3].get_xticklabels(), visible=True)
    fig.tight_layout()
    fig.savefig(pp, format='pdf')

print('Done.')

## List of top 15 mutations with increased frequency: [[oncomergeSummary[tumor]['pq_mff_mhc'][['Symbol','Delta_over_PAM']].dropna().sort_values('Delta_over_PAM').tail(15),tumor] for tumor in tumors]

## [[oncomergeSummary[tumor]['pq_mff_mhc'][['Symbol','Final_freq']].dropna().sort_values('Final_freq').query('Final_freq >= 0.5'),tumor] for tumor in tumors]

#%%
##################
## Cross cancer ##
##################

print('Cross cancer...')

crossCancer = []
for tumor in tumors:
    tmp1 = oncomergeSummary[tumor]['pq_mff'][['PAM_freq','Fusion_freq','CNA_freq','Act_freq','LoF_freq','Final_freq','Final_mutation_type']].loc[oncomergeSummary[tumor]['pq_mff']['Final_mutation_type'].dropna().index]
    tmp1['tumor'] = [tumor] * tmp1.shape[0]
    crossCancer.append(tmp1)

crossCancerDF = pd.concat(crossCancer).sort_index()
dupCrossCancerDF = crossCancerDF.loc[crossCancerDF.index.duplicated(keep=False)]
dupCrossCancerDF.to_csv('crossCancer_SomMuts_PQ_MFF.csv')
#genes1 = mj.querymany(list(set(dupCrossCancerDF.index)), scopes='entrezgene', fields='symbol', species='human', as_dataframe=True, df_index=True, entrezonly=True, verbose=False)
genes1 = [convertMe.index[convertMe==i][0] for i in list(set(dupCrossCancerDF.index)) if i in list(convertMe)]
forPlotting = {}
dupCrossCancer_hm1 = pd.DataFrame(data=0,index = sorted(set(dupCrossCancerDF['tumor'])), columns = [i for i in sorted(set(dupCrossCancerDF.index)) if sum(dupCrossCancerDF.index==i)>=5])
for dup1 in list(set(dupCrossCancerDF.index)):
    if sum(dupCrossCancerDF.index==dup1)>=5:
        cur1 = {'PAM':0,'Fusion':0,'CNAamp':0,'Act':0,'CNAdel':0,'LoF':0}
        tmp2 = dupCrossCancerDF.loc[dup1]['Final_mutation_type'].value_counts()
        for i in list(tmp2.index):
            cur1[i] = tmp2[i]
        #forPlotting[convertMe.index[convertMe==dup1][0]] = cur1
        forPlotting[dup1] = cur1
        for gn1, row1 in dupCrossCancerDF.loc[dup1].iterrows():
            dupCrossCancer_hm1.loc[row1['tumor'],dup1] = colors2[row1['Final_mutation_type']]

forPlotting = pd.DataFrame(forPlotting)
# Add gene names
#forPlotting.columns = [convertMe.index[convertMe==i][0] for i in forPlotting.columns if i in list(convertMe)]
#dupCrossCancer_hm1.columns = [convertMe.index[convertMe==i][0] for i in dupCrossCancer_hm1.columns if i in list(convertMe)]

# Make plot of cross-cancer somatic mutation genes by final mutation type
# %%
with PdfPages('crossCancer_gt5_8_10_2021.pdf') as pp:
    fig, ax = plt.subplots(1,1,figsize=(8,6)) #sharex='col'
    #fig = plt.figure(figsize=(8,3))
    #gs = GridSpec(1,7)
    # Plot final mutation counts
    #fig, ax = plt.subplots(figsize=(8,4))
    spacer1 = pd.DataFrame({'':{'PAM':0,'Fusion':0,'CNAamp':0,'Act':0,'CNAdel':0,'LoF':0}})
    try:   
        tmp1 = forPlotting.T[forPlotting.loc['LoF']>0].T
        #tmp1 = tmp1.T.drop(int(convertMe['KMT2C'])).T
        tmp1 = tmp1.T.drop(58508).T
    except:
        pass
    tmp1 = tmp1.reindex(tmp1.sum().sort_values().index, axis=1).T
    tmp2 = pd.DataFrame(forPlotting[58508].T)
    #tmp2 = pd.DataFrame(forPlotting[int(convertMe['KMT2C'])].T)
    tmp3 = forPlotting.T[forPlotting.loc['LoF']==0].T
    tmp3 = tmp3.reindex(tmp3.sum().sort_values().index, axis=1).T
    tmp4 = pd.concat([tmp1,spacer1.T,tmp2.T,spacer1.T,tmp3],axis=0)
    '''tmp4_index = []
    for i in tmp4.index:
        if i in list(convertMe):
            tmp4_index.append(convertMe.index[convertMe==i][0])
        else:
            tmp4_index.append('')
    tmp4.index = tmp4_index
    '''
    tmp5 = pd.Series(tmp4.index)
    tmp5[tmp4.index.isna()] = ''
    tmp4.index = tmp5
    tmp4.T.loc[['PAM','Fusion','CNAamp','Act','CNAdel','LoF']].T.plot(kind='bar', stacked=True, ax=ax, ylim=[0,25], color=[colors1[i] for i in ['PAM','Fusion','CNAamp','Act','CNAdel','LoF']])
    ax.set_ylabel('Tumors with mutation')
    fig.tight_layout()
    fig.savefig(pp, format='pdf')

# Heatmap
with PdfPages('crossCancer_gt5_8_10_2021_heatmap.pdf') as pp:
    fig, ax = plt.subplots(1,1,figsize=(10,5)) #sharex='col'
    #fig = plt.figure(figsize=(8,3))
    sns.set(font_scale=0.4)
    sns.heatmap(dupCrossCancer_hm1.loc[dupCrossCancer_hm1.T.sum()>0,sorted(list(tmp1.index))+list(tmp2.columns)+sorted(list(tmp3.index))], cmap=colors3, xticklabels=True, yticklabels=True)
    #sns.clustermap(dupCrossCancer_hm1.loc[dupCrossCancer_hm1.T.sum()>0], cmap=colors3, xticklabels=True, yticklabels=True)
    fig.tight_layout()
    fig.savefig(pp, format='pdf')

# Heatmap
'''with PdfPages('crossCancer_gt5_4_14_clustermap.pdf') as pp:
    #fig, ax = plt.subplots(1,1,figsize=(8,3)) #sharex='col'
    #fig = plt.figure(figsize=(8,3))
    fig = sns.clustermap(dupCrossCancer_hm1.loc[dupCrossCancer_hm1.T.sum()>0,], cmap=colors3, method='complete')
    #fig.tight_layout()
    fig.savefig(pp, format='pdf')
'''

# Both oncogene and tumor suppressor
forPlotting = {}
dupCrossCancer_hm1 = pd.DataFrame(data=0, index = sorted(set(dupCrossCancerDF['tumor'])), columns = list(set(dupCrossCancerDF.index)))
for dup1 in list(set(dupCrossCancerDF.index)):
    if sum(dupCrossCancerDF.index==dup1)>=0:
        cur1 = {'PAM':0,'Fusion':0,'CNAamp':0,'Act':0,'CNAdel':0,'LoF':0}
        tmp2 = dupCrossCancerDF.loc[dup1]['Final_mutation_type'].value_counts()
        for i in list(tmp2.index):
            cur1[i] = tmp2[i]
        #forPlotting[convertMe.index[convertMe==dup1][0]] = cur1
        forPlotting[dup1] = cur1
        for gn1, row1 in dupCrossCancerDF.loc[dup1].iterrows():
            dupCrossCancer_hm1.loc[row1['tumor'],dup1] = colors2[row1['Final_mutation_type']]

forPlotting = pd.DataFrame(forPlotting)
forPlotting = forPlotting[forPlotting.columns[((forPlotting.loc['Act']>0) | (forPlotting.loc['CNAamp']>0)) & ((forPlotting.loc['LoF']>0) | (forPlotting.loc['CNAdel']>0))]]
dupCrossCancer_hm1 = dupCrossCancer_hm1[forPlotting.columns[((forPlotting.loc['Act']>0) | (forPlotting.loc['CNAamp']>0)) & ((forPlotting.loc['LoF']>0) | (forPlotting.loc['CNAdel']>0))]]
# Add gene names
#forPlotting.columns = [convertMe.index[convertMe==i][0] for i in forPlotting.columns if i in list(convertMe)]
#dupCrossCancer_hm1.columns = [convertMe.index[convertMe==i][0] for i in dupCrossCancer_hm1.columns if i in list(convertMe)]


# Make plot of cross-cancer somatic mutation genes by final mutation type
with PdfPages('crossCancer_dualRole_8_10_2021.pdf') as pp:
    fig, ax = plt.subplots(1,1,figsize=(8,3))
    tmp1 = forPlotting
    tmp1 = tmp1.reindex(tmp1.sum().sort_values().index, axis=1).T
    tmp1.T.loc[['PAM','Fusion','CNAamp','Act','CNAdel','LoF']].T.plot(kind='bar', stacked=True, ax=ax, color=[colors1[i] for i in ['PAM','Fusion','CNAamp','Act','CNAdel','LoF']])
    ax.set_ylabel('Tumors with mutation')
    fig.tight_layout()
    fig.savefig(pp, format='pdf')


# Heatmap
with PdfPages('crossCancer_dualRole_8_10_2021_heatmap.pdf') as pp:
    fig, ax = plt.subplots(1,1,figsize=(3,3)) #sharex='col'
    sns.set(font_scale=0.4)
    sns.heatmap(dupCrossCancer_hm1.loc[dupCrossCancer_hm1.T.sum()>0,sorted(list(tmp1.index))], cmap=colors3, xticklabels=True, yticklabels=True)
    fig.tight_layout()
    fig.savefig(pp, format='pdf')


print('Done.')

#%%
###############################
## Tumor type specific plots ##
###############################

print('Tumor type specific plots...')

#tumor = 'ESCA'
with PdfPages('somaticMutations_byTumor_8_10_2021.pdf') as pp:
    for tumor in tumors:
        tmp1 = oncomergeSummary[tumor]['pq_mff'].loc[oncomergeSummary[tumor]['pq_mff']['Final_freq']>=0.1]
        tmp1.index = tmp1.index
        pam1 = tmp1['PAM_freq']
        pam1.loc[pam1.isna()] = 0
        cna1 = tmp1['CNA_freq']
        cna1.loc[cna1.isna()] = 0
        overlap1 = (pam1+cna1)-tmp1['Final_freq']
        pam1 = pam1-overlap1
        cna1 = cna1-overlap1
        pam1 = tmp1['PAM_freq']
        pam1.loc[pam1.isna()] = 0
        cna1 = tmp1['CNA_freq']
        cna1.loc[cna1.isna()] = 0
        overlap1 = (pam1+cna1)-tmp1['Final_freq']
        overlap1[overlap1<=0] = 0
        pam1 = pam1-overlap1
        cna1 = cna1-overlap1
        if tmp1.shape[0]>1:
            cnaDel1 = cna1.loc[tmp1['CNA_type']=='Del']
            cnaDel1[np.isnan(cnaDel1)] = 0
            cnaAmp1 = cna1.loc[tmp1['CNA_type']=='Amp']
            cnaAmp1[np.isnan(cnaAmp1)] = 0
        else:
            if list(tmp1['CNA_type'])[0]=='Del':
                cnaDel1 = cna1
                cnaAmp1 = pd.Series([np.nan], index=cna1.index)
            elif list(tmp1['CNA_type'])[0]=='Amp':
                cnaDel1 = pd.Series([np.nan], index=cna1.index)
                cnaAmp1 = cna1
            else:
                cnaDel1 = pd.Series([np.nan], index=cna1.index)
                cnaAmp1 = pd.Series([np.nan], index=cna1.index)
        #genes1 = mj.querymany(list(tmp1.index), scopes='entrezgene', fields='symbol', species='human', as_dataframe=True, df_index=True, entrezonly=True, verbose=False)
        df1 = pd.DataFrame([pam1, overlap1, cnaAmp1, cnaDel1], index=['PAM','Overlap','CNAamp','CNAdel']).T
        #df1.reindex(tmp1['Symbol'])

        # Make plot of cross-cancer somatic mutation genes by final mutation type
        fig, ax = plt.subplots(1,1,figsize=(8,3)) #sharex='col'
        df1.T.reindex(df1.T.sum().sort_values().index, axis=1).T.plot(kind='bar', stacked=True, ax=ax, title=tumor, ylim=[0,1], color=[colors1[i] for i in df1.columns])
        ax.set_xlabel('Somatically mutated genes')
        ax.set_ylabel('Mutation frequency')
        fig.tight_layout()
        fig.savefig(pp, format='pdf')
# %%
# ['TP53','CDKN2A','RB1','SMAD4','MYC','KRAS','NRAS','ERBB2','BRAF','KMT2C']
genes2 = [7157,1029, 5925,4089,4609,3845,4893,2064,673, 58508]
gdt = {}
for gene2 in genes2:
    tmp1 = []
    tumor1 = []
    for tumor in tumors:
        tmp2 = oncomergeSummary[tumor]['pq_mff']
        tmp2.index = tmp2.index
        if gene2 in tmp2.index and not np.isnan(tmp2.loc[gene2]['Final_freq']):
            tmp1.append(tmp2.loc[gene2])
            tumor1.append(tumor)
    gdt[gene2] = pd.DataFrame(tmp1)
    gdt[gene2].index = tumor1
# %%
with PdfPages('somaticMutations_byGene_8_10_2021.pdf') as pp:
    for gene2 in genes2:
        tmp1 = gdt[gene2]
        pam1 = tmp1['PAM_freq']
        pam1.loc[pam1.isna()] = 0
        cna1 = tmp1['CNA_freq']
        cna1.loc[cna1.isna()] = 0
        overlap1 = (pam1+cna1)-tmp1['Final_freq']
        overlap1[overlap1<=0] = 0
        pam1 = pam1-overlap1
        cna1 = cna1-overlap1
        cnaDel1 = cna1.loc[tmp1['CNA_type']=='Del']
        cnaAmp1 = cna1.loc[tmp1['CNA_type']=='Amp']
        df1 = pd.DataFrame([pam1, overlap1, cnaAmp1, cnaDel1], index=['PAM','Overlap','CNAamp','CNAdel']).T
        df1['CNAamp'][np.isnan(df1['CNAamp'])] = 0
        df1['CNAdel'][np.isnan(df1['CNAdel'])] = 0
        df1.index = gdt[gene2].index

        # Make plot of cross-cancer somatic mutation genes by final mutation type
        fig, ax = plt.subplots(1,1,figsize=(8,3)) #sharex='col'
        df1.T.reindex(df1.T.sum().sort_values().index, axis=1).T.plot(kind='bar', stacked=True, ax=ax, title=gene2, color=[colors1[i] for i in df1.columns])
        ax.set_xlabel('Tumor type')
        ax.set_ylabel('Mutation frequency')
        fig.tight_layout()
        fig.savefig(pp, format='pdf')

crossCancer = []
for tumor in tumors:
    tmp1 = oncomergeSummary[tumor]['pq_mff'][['PAM_freq','CNA_freq','Act_freq','LoF_freq','Final_freq','Final_mutation_type']].loc[oncomergeSummary[tumor]['pq_mff']['Final_mutation_type'].dropna().index]
    tmp1['tumor'] = [tumor] * tmp1.shape[0]
    crossCancer.append(tmp1)

crossCancerDF = pd.concat(crossCancer).sort_index()
dupCrossCancerDF = crossCancerDF.loc[crossCancerDF.index.duplicated(keep=False)]
dupCrossCancerDF.to_csv('crossCancer_SomMuts_PQ_MFF.csv')
#genes1 = mj.querymany(list(set(dupCrossCancerDF.index)), scopes='entrezgene', fields='symbol', species='human', as_dataframe=True, df_index=True, entrezonly=True, verbose=False)
forPlotting = {}
for dup1 in list(set(dupCrossCancerDF.index)):
    if sum(dupCrossCancerDF.index==dup1)>=5:
        cur1 = {'PAM':0,'CNAamp':0,'Act':0,'CNAdel':0,'LoF':0}
        tmp2 = dupCrossCancerDF.loc[dup1]['Final_mutation_type'].value_counts()
        for i in list(tmp2.index):
            cur1[i] = tmp2[i]
        forPlotting[dup1] = cur1

forPlotting = pd.DataFrame(forPlotting)
# ['TP53','CDKN2A','RB1','MYC','ERBB2','BRAF']
genes2 = [7157,1029, 5925,4609,2064,673]
gdt = {}
for gene2 in genes2:
    tmp1 = []
    tumor1 = []
    for tumor in tumors:
        tmp2 = oncomergeSummary[tumor]['pq_mff']
        tmp2.index = tmp2.index
        if gene2 in tmp2.index and not np.isnan(tmp2.loc[gene2]['Final_freq']):
            tmp1.append(tmp2.loc[gene2])
            tumor1.append(tumor)
    gdt[gene2] = pd.DataFrame(tmp1)
    gdt[gene2].index = tumor1

params = {'legend.fontsize': 7,
          'legend.handlelength': 2}
plt.rcParams.update(params)
# Make plot of cross-cancer somatic mutation genes by final mutation type
with PdfPages('crossCancer_gt5_genes_8_10_2021.pdf') as pp:
    fig, ax = plt.subplots(7,1,figsize=(8,11)) #sharex='col'
    #fig = plt.figure(figsize=(8,3))
    #gs = GridSpec(1,7)
    # Plot final mutation counts
    #fig, ax = plt.subplots(figsize=(8,4))
    spacer1 = pd.DataFrame({'':{'PAM':0,'CNAamp':0,'Act':0,'CNAdel':0,'LoF':0}})
    try:   
        tmp1 = forPlotting.T[forPlotting.loc['LoF']>0].T
        #tmp1 = tmp1.T.drop(int(convertMe['KMT2C'])).T
        tmp1 = tmp1.T.drop(58508).T
    except:
        pass
    tmp1 = tmp1.reindex(tmp1.sum().sort_values().index, axis=1).T
    try:   
        tmp2 = pd.DataFrame(forPlotting[58508].T)
    except:
        pass
    tmp3 = forPlotting.T[forPlotting.loc['LoF']==0].T
    #tmp3 = tmp3.drop(labels='RNA5SP251',axis=1)
    tmp3 = tmp3.reindex(tmp3.sum().sort_values().index, axis=1).T
    tmp4 = pd.concat([tmp1,spacer1.T,tmp2.T,spacer1.T,tmp3],axis=0)
    tmp4.T.loc[['PAM','CNAamp','Act','CNAdel','LoF']].T.plot(kind='bar', stacked=True, ax=ax[0], ylim=[0,25], color=[colors1[i] for i in df1.columns])
    ax[0].set_ylabel('Tumors with\nmutation')
    #fig.tight_layout()
    #fig.savefig(pp, format='pdf')

    for cur1 in range(len(genes2)):
        tmp1 = gdt[genes2[cur1]]
        pam1 = tmp1['PAM_freq']
        pam1.loc[pam1.isna()] = 0
        cna1 = tmp1['CNA_freq']
        cna1.loc[cna1.isna()] = 0
        overlap1 = (pam1+cna1)-tmp1['Final_freq']
        overlap1[overlap1<=0] = 0
        pam1 = pam1-overlap1
        cna1 = cna1-overlap1
        cnaDel1 = cna1.loc[tmp1['CNA_type']=='Del']
        cnaAmp1 = cna1.loc[tmp1['CNA_type']=='Amp']
        df1 = pd.DataFrame([pam1, overlap1, cnaAmp1, cnaDel1], index=['PAM','Overlap','CNAamp','CNAdel']).T
        df1['CNAamp'][np.isnan(df1['CNAamp'])] = 0
        df1['CNAdel'][np.isnan(df1['CNAdel'])] = 0
        df1.index = gdt[genes2[cur1]].index

        # Make plot of cross-cancer somatic mutation genes by final mutation type
        df1.T.reindex(df1.T.sum().sort_values().index, axis=1).T.plot(kind='bar', stacked=True, ax=ax[cur1+1], title=genes2[cur1], color=[colors1[i] for i in df1.columns]) #, ylim=[0,0.6])
        #ax[cur1+1].set_xlabel('Tumor type')
        ax[cur1+1].set_ylabel('Mutation frequency')
    fig.tight_layout()
    fig.savefig(pp, format='pdf')

print('Done.')


# Load up MSI Hypermutator censored samples
mhcSamples = pd.read_csv('C:/Users/nihan/Dropbox/OncoMerge_pathways/TCGA/blacklist/blacklist_29850653_29625053.csv')
allSamples = {i:pamTumors[i].shape[1] for i in pamTumors}
plt.style.use('default')
with PdfPages('MSI_Hypermutation_8_10_2021.pdf') as pp:
    # Old order
    fig, ax = plt.subplots(1, 1, sharex=True, figsize=(8,1))
    tmp1 = mhcSamples.groupby(['Reason','Tumor Type'])['Tumor Type'].count().unstack('Tumor Type').fillna(0)
    tmp1[['DLBC','THCA','KIRP','TGCT','PCPG','UVM']] = [[0]*6]*3
    old = ['THCA','KICH','KIRP','THYM','KIRC','MESO','TGCT','PCPG','UVM','PAAD','BRCA','GBM','LGG','ACC','CHOL','LIHC','READ','PRAD','DLBC','ESCA','OV','HNSC','SARC','UCS','CESC','LUSC','LUAD','BLCA','UCEC','STAD','COAD']
    (tmp1[old]/[allSamples[i] for i in old]).loc[['MSI','MSI/Hypermutator','Hypermutator']].T.plot(ylabel='Frequncy of patient tumors', kind='bar', stacked=True, color=['#00429d', '#dfdfc1', '#6592bf'], ax=ax)
    fig.savefig(pp, format='pdf')
    plt.close()

    # New order
    fig, ax = plt.subplots(1, 1, sharex=True, figsize=(8,1))
    tmp1 = mhcSamples.groupby(['Reason','Tumor Type'])['Tumor Type'].count().unstack('Tumor Type').fillna(0)
    tmp1[['DLBC','THCA','KIRP','TGCT','PCPG','UVM']] = [[0]*6]*3
    new = ['THCA','KICH','KIRP','THYM','KIRC','MESO','TGCT','PCPG','UVM','PAAD','BRCA','LGG','GBM','CHOL','ACC','LIHC','READ','PRAD','CESC','DLBC','COAD','HNSC','ESCA','OV','UCS','SARC','UCEC','STAD','LUSC','LUAD','BLCA']
    (tmp1[new]/[allSamples[i] for i in new]).loc[['MSI','MSI/Hypermutator','Hypermutator']].T.plot(ylabel='Frequncy of patient tumors', kind='bar', stacked=True, color=['#feb24c', '#ffeda0', '#f03b20'], ax=ax)
    fig.savefig(pp, format='pdf')
    plt.close()


tmp2 = copy.copy(summaryTable['pq_mff'])
del tmp2['Patients']
tmp3 = copy.copy(summaryTable['pq_mff'])
del tmp3['Patients']
tmp1 = mhcSamples.groupby(['Reason','Tumor Type'])['Tumor Type'].count().unstack('Tumor Type').fillna(0)
tmp1[['DLBC','THCA','KIRP','TGCT','PCPG','UVM']] = [[0]*6]*3
pearsonr(tmp2.loc[new].T.sum(),(tmp1[new]/[allSamples[i] for i in new]).sum())
pearsonr(tmp3.loc[new].T.sum(),(tmp1[new]/[allSamples[i] for i in new]).sum())
with PdfPages('MSI_Hypermutation_scatter_8_10_2021.pdf') as pp:
    # Old order
    fig, ax = plt.subplots(1, 2, sharex=True, figsize=(8,4))
    ax[0].scatter(list(tmp2.loc[new].T.sum()),list((tmp1[new]/[allSamples[i] for i in new]).sum()))
    # New order
    ax[1].scatter(list(tmp3.loc[new].T.sum()),list((tmp1[new]/[allSamples[i] for i in new]).sum()))
    fig.savefig(pp, format='pdf')
    plt.close()


# %%
