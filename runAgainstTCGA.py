from subprocess import *
import os

# Run for each tumor type in TCGA
for tumor in ["ACC", 'BLCA', 'CESC', 'CHOL', 'COAD', 'DLBC', "BRCA", "PCPG","READ","SARC","SKCM","TGCT","THCA","THYM","UCEC","UCS","UVM","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD", "PRAD", "STAD"]:
    print(tumor)
    # Make the output directory if it doesn't exists already
    if not (os.path.exists('/output')):
        os.makedirs('/output')
    if not (os.path.exists('/output/'+tumor)):
        os.makedirs('/output/'+tumor+"/pq_mff")
    # Run with both filters
    print('\nRunning pq_mff...')
    cmd4 = ['python Untitled-2.py',
            '-gp C:/Users/nihan/Dropbox/OncoMerge_pathways/TCGA/GISTIC/'+tumor,
            '-aaf C:/Users/nihan/Dropbox/OncoMerge_pathways/TCGA/OncoMerge_input_g2e_converter.csv',
            '-ln "Locus ID"',
            '-pam C:/Users/nihan/Dropbox/OncoMerge_pathways/TCGA/PAM/'+tumor+'_somMutMC3.csv',
            '-mscv C:/Users/nihan/Dropbox/OncoMerge_pathways/TCGA/MutSig2cv/'+tumor+'_sig2cv.csv',
            '-fus C:/Users/nihan/Dropbox/OncoMerge_pathways/TCGA/FUSIONS/'+tumor+'_fusions.csv',
            '-op /output/'+tumor+'/pq_mff',
            '-pq 0.1',
            '-mlg 10',
            '-mpf 0',
            '-tcga True',
            '-bl C:/Users/nihan/Dropbox/OncoMerge_pathways/TCGA/blacklist/blacklist_29850653_29625053.csv']
    proc4 = Popen(' '.join(cmd4), shell=True, stdin=PIPE)
    out = proc4.communicate()

