# MOB Spatial Transcriptomics Analysis with SPRI

# Dataset explanation

SPATIAL research obtained spatial expression data with a Spatial Transcriptomics (ST) platform recent years. The geometric representations of ST is rectangular. 

Figure placeholder

First, we read the MOB Replicate 11 data and pre-process it.

```{r, mob-qc, fig.width=8, fig.height=3}
    if dataname == 'MOB11':
        file_path = 'raw_data/Rep11_MOB_count_matrix-1.tsv' # MOB11
```

Then, we compute the total information coefficient (TIC) for each gene and rank spatial genes according to TIC values.

```{r, mob-spatially-unaware, fig.width=8, fig.height=4}
MIC_list, TIC_list, genes = compute_genes_MIC(file_path=file_path)
MIC_list, TIC_list, genes = (list(t) for t in zip(*sorted(zip(MIC_list, TIC_list, genes),key=lambda x: x[1],reverse=True)))
df = pd.DataFrame({'gene':genes,'mic':MIC_list,'tic':TIC_list})
save_path = 'results/'+dataname+'/'+dataname+'_tic_rank.csv'
df.to_csv(save_path,index=False)

```

Identify statistically significant SE genes with permutation test.(optional)

```{r, mob-diff-gexp}
p_values, scores, genes = compute_genes_pvalue(file_path=file_path,top_genes = top_genes, method='TIC',num_perm = 10000,num_jobs=40)
scores_rank, p_values_rank, genes_rank, = (list(t) for t in zip(*sorted(zip(scores,p_values,genes),key=lambda x: x[1])))
df = pd.DataFrame({'gene':genes_rank,'pvalue':p_values_rank,'TIC':scores_rank})
df.to_csv(save_path)
```

After getting the P values for each gene, we used Benjamini-Yekutieli in the python package of “statsmodels” to control FDR.

```{r, mob-diff-gexp-plot, fig.width=4, fig.height=4}
genes, combined_pvalue, pvalue_corrected, tic_list = compute_FDR(df)
genes, combined_pvalue, pvalue_corrected, tic_list = (list(t) for t in zip(*sorted(zip(genes, combined_pvalue,pvalue_corrected,tic_list),key=lambda x: float(x[2]))))
df = pd.DataFrame({'gene':genes,'pvalue':combined_pvalue,'qvalue':pvalue_corrected,'TIC':tic_list})
save_path = 'results/'+dataname+'/'+dataname+'_pvals.csv'
df.to_csv(save_path)
```
