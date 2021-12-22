import numpy as np
import pandas as pd
import statsmodels
from statsmodels.stats.multitest import fdrcorrection
from statsmodels.stats.multitest import multipletests

gene = []
combined_pvalue = []
adjusted_pvalue = []

def compute_FDR(df):
    genes = df['Unnamed: 0'].values
    pvalue = df['p.value'].values
    

    # sort
    pvalue, genes= (list(t) for t in zip(*sorted(zip(pvalue,genes),key=lambda x: x[0])))
    print(pvalue[0:10])

    rejected,qvalue = fdrcorrection(np.array(pvalue), alpha=0.05, method='negcorr', is_sorted=True)
    return(genes,pvalue,qvalue)


name ='BC2020' # MOB11, MOB12, BC, BC2020
prefix = 'output/'

read_path = prefix + name+'_MERINGUE.csv'
save_path = prefix+ name+'_MERINGUE_fdrcorrection.csv'

df_BC = pd.read_csv(read_path)
print(df_BC)
genes,combined_pvalue, pvalue_corrected = compute_FDR(df_BC)
genes, combined_pvalue,pvalue_corrected = (list(t) for t in zip(*sorted(zip(genes, combined_pvalue,pvalue_corrected),key=lambda x: float(x[2]))))
df = pd.DataFrame({'gene':genes,'pvalue':combined_pvalue,'qvalue':pvalue_corrected})
df.to_csv(save_path)

