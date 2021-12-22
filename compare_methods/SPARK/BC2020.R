##-------------------------------------------------------------
## SeqFISH Data Analysis
##-------------------------------------------------------------
rm(list = ls())

library(SPARK)

counts <- t(read.csv("processed_data/BC2020_count_ambigous.csv", row.names = 1, 
    check.names = F))
info <- read.csv("processed_data/BC2020_info_ambiguous.csv", row.names = 1)
spark <- CreateSPARKObject(counts = counts, location = info[, 1:2], 
    percentage = 0.1, min_total_counts = 10)

spark@lib_size <- apply(spark@counts, 2, sum)
spark <- spark.vc(spark, covariates = NULL, lib_size = spark@lib_size, 
    num_core = 10, verbose = T, fit.maxiter = 500)
spark <- spark.test(spark, check_positive = T, verbose = T)

# save(spark, file = "./output/seqFISH_seqFISH+_SPARK.rds")
# spark@res_mtest[,c("combined_pvalue","adjusted_pvalue")]
write.csv(spark@res_mtest[,c("combined_pvalue","adjusted_pvalue")], file = "output/BC2020_SPARK_qvalue.csv", row.names = T)