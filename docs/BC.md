# Breast Cancer Spatial Transcriptomics Analysis with SPRI

# Dataset explanation

SPATIAL research obtained spatial expression data with a Spatial Transcriptomics (ST) platform recent years. The geometric representations of ST is rectangular. 

The slice "Breast cancer layer 2"  data has been prepared for you and is available in the folder "raw_data".

First, we read the "Breast cancer layer 2" data and pre-process it.
1) selected spots with at least ten total read counts
2) filtered out genes that are expressed in less than 10% of the array spots 

```{r, mob-qc, fig.width=8, fig.height=3}
    if dataname == 'BC':
        file_path = 'raw_data/Layer2_BC_count_matrix-1.tsv'

    pd_file = pd.read_csv(file_path, sep='\t')
    cell_x,cell_y,pd_file = filter_df(pd_file)
    cortext_svz_counts = pd_file
```

Secondlly, background correction is performed..

```{r, mob-spatially-unaware, fig.width=8, fig.height=4}
pd_BG = background_correction(cell_x,cell_y)

# You can plot the original spatial coordinates and the added background coordinates to see if the background correction is correct
import matplotlib.pyplot as plt
plt.scatter(cell_x,cell_y, edgecolor='b', facecolor='none', alpha=0.5,s=1)
plt.xlabel("x"); plt.ylabel("y")
plt.savefig('original_spot.png')
plt.close()

cell_BG_x = pd_BG['cell_x_noise'].values.tolist()
cell_BG_y = pd_BG['cell_y_noise'].values.tolist()
plt.scatter(cell_BG_x,cell_BG_y, edgecolor='b', facecolor='none', alpha=0.5,s=1)
plt.xlabel("x"); plt.ylabel("y")
plt.savefig('background_spot.png')
plt.close()

cell_x_final = cell_x.tolist()+ cell_BG_x
cell_y_final = cell_y.tolist() + cell_BG_y
plt.scatter(cell_x_final,cell_y_final, edgecolor='b', facecolor='none', alpha=0.5,s=1)
plt.xlabel("x"); plt.ylabel("y")
plt.savefig('final_spot.png')
plt.close()

```

Then, we compute the TIC (total information coefficient) for each gene and rank genes according to the TIC values.

```{r, mob-spatially-unaware, fig.width=8, fig.height=4}
MIC_list, TIC_list, genes = compute_genes_MIC(cell_x = cell_x,cell_y = cell_y,pd_file = pd_file,pd_BG = pd_BG)
MIC_list, TIC_list, genes = (list(t) for t in zip(*sorted(zip(MIC_list, TIC_list, genes),key=lambda x: x[1],reverse=True)))
df = pd.DataFrame({'gene':genes,'mic':MIC_list,'tic':TIC_list})
save_path = 'results/'+dataname+'/'+dataname+'_tic_rank.csv'
df.to_csv(save_path,index=False)
```

Identify statistically significant SE genes with permutation test.(optional)

```{r, mob-diff-gexp}
# you can choose the number of permutation times and the number of cpus.
p_values, scores, genes = compute_genes_pvalue(file_path=file_path,top_genes = genes, method='TIC',num_perm = 10000,num_jobs=40)
scores_rank, p_values_rank, genes_rank, = (list(t) for t in zip(*sorted(zip(scores,p_values,genes),key=lambda x: x[1])))
df = pd.DataFrame({'gene':genes_rank,'pvalue':p_values_rank,'TIC':scores_rank})
```

After getting the P values for each gene, we used Benjamini-Yekutieli in the python package of “statsmodels” to control FDR.

```{r, mob-diff-gexp-plot, fig.width=4, fig.height=4}
genes, combined_pvalue, pvalue_corrected, tic_list = compute_FDR(df)
genes, combined_pvalue, pvalue_corrected, tic_list = (list(t) for t in zip(*sorted(zip(genes, combined_pvalue,pvalue_corrected,tic_list),key=lambda x: float(x[2]))))
df = pd.DataFrame({'gene':genes,'pvalue':combined_pvalue,'qvalue':pvalue_corrected,'TIC':tic_list})
save_path = 'results/'+dataname+'/'+dataname+'_pvals.csv'
df.to_csv(save_path)
```
