# SPRI: Spatial Pattern Recognition using Information based method for spatial gene expression data 

The rapid development of spatially resolved transcriptomics has made it possible to analyze spatial gene expression patterns in complex biological tissues. To identify such genes, we propose a novel and robust nonparametric information-based approach, SPRI, to recognize their spatial patterns. SPRI directly models spatial transcriptome raw count data without model assumptions, which transforms the problem of spatial expression pattern recognition into the detection of dependencies between spatial coordinate pairs with gene read count as the observed frequencies. 

## Framework

![image](https://github.com/xiaoyeye/SPRI/blob/main/figure/figure1.png)


## Requirements 

numpy==1.16.2

pandas==0.24.2

minepy==1.2.5

statsmodels==0.9.0

scanpy==1.8.2


## Tutorials

1. [Breast cancer Analysis](https://github.com/xiaoyeye/SPRI/blob/main/docs/BC.md)

2. [10x Visium Analysis](https://github.com/xiaoyeye/SPRI/blob/main/docs/10x_Visium.md)



## Data:

ST data is avalieble at https://www.spatialresearch.org/.

LIBD human dorsolateral prefrontal cortex, dorsolateral prefrontal cortex 10x Visium data data is avalieble at http://research.libd.org/spatialLIBD/).

Adult Mouse Brain and Kidney data is avalieble at https://www.10xgenomics.com/resources/datasets/).
