
# 10x Visium Analysis with SPRI


# Dataset explanation

10X genomics obtained spatial expression data with a a Visium platform. The geometric representations of ST is rectangular. 
10X genomics recently launched a new platform to obtain spatial expression data using a Visium Spatial Gene Expression slide. Visium platform produces a sampling distribution that is hexagonal. The geometric representations of Visium technologt is hexagonal. 

Visium technology:

Figure placeholder

In this tutorial, LIBD data was tested.  The LIBD of slice 151507 data to run this tutorial can be found [here](http://science.sciencemag.org/content/353/6294/78).

The slice 151507 data has been prepared for you and is available in the folder "raw_data".

```{r, mob-init}
dataname = '151507'
save_prefix = 'processed_data/' + dataname 

# generate save path
if os.path.exists(save_prefix):
    print(save_prefix)
else:
    os.mkdir(save_prefix)

save_prefix = save_prefix + '/'
```

First, we read the 10x Visium data and pre-process it. Secificlly, we filtered out spots less than ten total read counts and selected genes expressed at least 10% of the spots. 

```{r, mob-qc, fig.width=8, fig.height=3}
# read data
from scanpy import read_10x_h5

#Read data : 

# load gene expression: expressiion_matrix
save_prefix = 'processed_data/' + dataname 

adata = read_10x_h5('data/LIBD/' + dataname + '/'+dataname+'_filtered_feature_bc_matrix.h5') #151673_filtered_feature_bc_matrix.h5 is expression_matrix.h5
spatial=pd.read_csv('data/LIBD/'+ dataname + '/tissue_positions_list.txt',sep=",",header=None,na_filter=False,index_col=0) 
adata.obs["x1"]=spatial[1]
adata.obs["x2"]=spatial[2]
adata.obs["x3"]=spatial[3]
adata.obs["x4"]=spatial[4]
adata.obs["x5"]=spatial[5]
adata.obs["x_array"]=adata.obs["x2"]
adata.obs["y_array"]=adata.obs["x3"]
adata.obs["x_pixel"]=adata.obs["x4"]
adata.obs["y_pixel"]=adata.obs["x5"]

adata=adata[adata.obs["x1"]==1]
adata.var_names=[i.upper() for i in list(adata.var_names)]
adata.var["genename"]=adata.var.index.astype("str") # 33538 genes

x_array=adata.obs["x_array"].tolist() 
y_array=adata.obs["y_array"].tolist()
x_pixel=adata.obs["x_pixel"].tolist() 
y_pixel=adata.obs["y_pixel"].tolist() 

# generate cell location
location_list = []
for i in range(len(x_pixel)):
    location = str(x_pixel[i]) + 'x' + str(y_pixel[i])
    location_list.append(location)

# preprocess
adata.var_names_make_unique() #(3639, 33538) -> (3639, 33538)

cell_names = adata.obs_names.tolist() # 3639
gene_names = adata.var["genename"].tolist() #19130
expressiion_matrix = adata.X.A  #(3639, 33538) -> (3639, 19130) after data preprocess
# expressiion_matrix = adata.X
print('cell_names: ',len(cell_names))
print('gene_names: ',len(gene_names))
print('x_pixel: ',len(x_pixel))
print('y_pixel: ',len(y_pixel))
print('expressiion_matrix: ',np.shape(expressiion_matrix))

# generate gene expression file
df_expression = pd.DataFrame(expressiion_matrix)
df_expression.index = location_list
df_expression.columns = gene_names
df_expression = df_expression.reset_index() #(2264, 19465)
df_expression.to_csv(save_prefix + dataname + '_count_matrix.csv')
print(df_expression)
# filter spots and genes
cell_x,cell_y,pd_file = filter_df(df_expression) # (2264, 19465)
# background correction
pd_BG = background_correction(cell_x,cell_y)
```

Then, background correction was applied. In other words, we added uniformly distributed background spot locations to convert the spots shape to hexagonal.

```{r, mob-qc, fig.width=8, fig.height=3}
# background correction
pd_BG = background_correction(cell_x,cell_y)
```

plot the original spots and the corrected background spots........

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
save_path = 'results/'+dataname+'/'+dataname+'_pvals.csv'
df.to_csv(save_path)
```

After getting the P values for each gene, we used Benjamini-Yekutieli in the python package of “statsmodels” to control FDR.

```{r, mob-diff-gexp-plot, fig.width=4, fig.height=4}
# fdr correction: q-val
genes, combined_pvalue, pvalue_corrected, tic_list = compute_FDR(df)
genes, combined_pvalue, pvalue_corrected, tic_list = (list(t) for t in zip(*sorted(zip(genes, combined_pvalue,pvalue_corrected,tic_list),key=lambda x: float(x[2]))))
df = pd.DataFrame({'gene':genes,'pvalue':combined_pvalue,'qvalue':pvalue_corrected,'TIC':tic_list})
save_path = 'results/'+dataname+'/'+dataname+'_pvals.csv'
df.to_csv(save_path)
```
