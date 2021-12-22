# SPRI: Spatial Pattern Recognition using Information based method for spatial gene expression data 

The rapid development of spatially resolved transcriptomics has made it possible to analyze spatial gene expression patterns in complex biological tissues. To identify such genes, we propose a novel and robust nonparametric information-based approach, SPRI, to recognize their spatial patterns. SPRI directly models spatial transcriptome raw count data without model assumptions, which transforms the problem of spatial expression pattern recognition into the detection of dependencies between spatial coordinate pairs with gene read count as the observed frequencies. 

Framework

![image](https://github.com/xiaoyeye/SPRI/blob/main/figure/figure1.png)


# 1. Requirements 

numpy==1.16.2

pandas==0.24.2

minepy==1.2.5

statsmodels==0.9.0


# 2. Example

2.1 Raw data is put in the folder "raw_data".

2.2 Run 'SPRI_MOB11.py' to compute TIC for MOB Replicate 11 data. Statisticlly significant SE genes with permutation test also can be obtained if need.


# 3. Download all datasets used in SPRI:

Data is avalieble at https://www.spatialresearch.org/ 

