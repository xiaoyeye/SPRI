
suppressMessages(library(MERINGUE))

# BC
read_path = "raw_data/Layer2_BC_count_matrix-1.tsv"
save_path = "output/BC_MERINGUE.csv"


counts <- read.table(read_path, check.names = F)
rn <- rownames(counts)

# extract the coordinates from the rawdata
info <- cbind.data.frame(x = as.numeric(sapply(strsplit(rn, split = "x"), 
    "[", 1)), y = as.numeric(sapply(strsplit(rn, split = "x"), "[", 2)))
rownames(info) <- rn

counts <- t(counts)

pos <- info
cd <- counts

# Remove poor datasets and genes
counts <- cleanCounts(counts = cd, 
                      min.reads = 100, 
                      min.lib.size = 100, 
                      plot=TRUE,
                      verbose=TRUE)


pos <- pos[colnames(counts),]

# CPM normalize
mat <- normalizeCounts(counts = counts, 
                       log=FALSE,
                       verbose=TRUE)

w <- getSpatialNeighbors(pos, filterDist = 2.5)
plotNetwork(pos, w)

# Identify sigificantly spatially auto-correlated genes
start_time <- Sys.time()
I <- getSpatialPatterns(mat, w)
end_time <- Sys.time()
print(end_time - start_time)
write.csv(I,save_path)

results.filter <- filterSpatialPatterns(mat = mat,
                                        I = I,
                                        w = w,
                                        adjustPv = TRUE,
                                        alpha = 0.05,
                                        minPercentCells = 0.05,
                                        verbose = TRUE)

# write.csv(results.filter,"output/MOB_filter.csv")