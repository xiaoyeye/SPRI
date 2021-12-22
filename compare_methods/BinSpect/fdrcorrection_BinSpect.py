import numpy as np
import pandas as pd
import statsmodels
from statsmodels.stats.multitest import fdrcorrection

gene = []
combined_pvalue = []
adjusted_pvalue = []

def compute_FDR(df):
    genes = df['genes'].values
    pvalue = df['p.value'].values

    # sort
    pvalue, genes= (list(t) for t in zip(*sorted(zip(pvalue,genes),key=lambda x: x[0])))

    rejected,qvalue = fdrcorrection(np.array(pvalue), alpha=0.05, method='negcorr', is_sorted=True) #'negcorr' is the right version
    return(genes,pvalue,qvalue)


name ='MOB11' # MOB11, MOB12, BC, BC2020
Binspect_type = 'kmeans' 
prefix = 'output/'+name+'/'

read_path = prefix + name+'_'+Binspect_type+'.csv'
save_path = prefix+ name+'_'+Binspect_type+'_fdrcorrection.csv'

df_BC = pd.read_csv(read_path)
print(df_BC)
genes,combined_pvalue, pvalue_corrected = compute_FDR(df_BC)
genes, combined_pvalue,pvalue_corrected = (list(t) for t in zip(*sorted(zip(genes, combined_pvalue,pvalue_corrected),key=lambda x: float(x[2]))))
df = pd.DataFrame({'gene':genes,'pvalue':combined_pvalue,'qvalue':pvalue_corrected})
df.to_csv(save_path)

