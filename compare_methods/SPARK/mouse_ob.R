##-------------------------------------------------------------
## Breast Cancer Data Analysis
##-------------------------------------------------------------

rm(list = ls())

library(SPARK)

# counts <- read.table("raw_data/Rep11_MOB_count_matrix-1.tsv", check.names = F)
counts <- read.table("raw_data/Rep12_MOB_count_matrix-1.tsv", check.names = F)
rn <- rownames(counts)

# extract the coordinates from the rawdata
info <- cbind.data.frame(x = as.numeric(sapply(strsplit(rn, split = "x"), 
    "[", 1)), y = as.numeric(sapply(strsplit(rn, split = "x"), "[", 2)))
rownames(info) <- rn

# filter genes and cells/spots and 
spark <- CreateSPARKObject(counts = t(counts), location = info[, 1:2], 
     percentage = 0.1,min_total_counts = 10)

## total counts for each cell/spot
spark@lib_size <- apply(spark@counts, 2, sum)
## Estimating Parameter Under Null
spark <- spark.vc(spark, covariates = NULL, lib_size = spark@lib_size, 
    num_core = 10, verbose = T, fit.maxiter = 500)
## Calculating pval
spark <- spark.test(spark, check_positive = T, verbose = T)

head(spark@res_mtest[,c("combined_pvalue","adjusted_pvalue")])

# write.csv(spark@res_mtest[,c("combined_pvalue","adjusted_pvalue")], file = "output/MOB_SPARK_rep11.csv", row.names = T)
write.csv(spark@res_mtest[,c("combined_pvalue","adjusted_pvalue")], file = "output/MOB_SPARK_rep12.csv", row.names = T)
