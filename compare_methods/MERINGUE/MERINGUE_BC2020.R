
suppressMessages(library(MERINGUE))

counts <- t(read.csv("processed_data/BC2020_count_ambigous.csv", row.names = 1, 
    check.names = F))
info <- read.csv("processed_data/BC2020_info_ambiguous.csv", row.names = 1)

pos <- info
cd <- counts

# Remove poor datasets and genes
counts <- cleanCounts(counts = cd, 
                      min.reads = 1, 
                      min.lib.size = 1, 
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
save_path = "output/BC2020_MERINGUE.csv"
write.csv(I,save_path)

results.filter <- filterSpatialPatterns(mat = mat,
                                        I = I,
                                        w = w,
                                        adjustPv = TRUE,
                                        alpha = 0.05,
                                        minPercentCells = 0.05,
                                        verbose = TRUE)

# write.csv(results.filter,"output/MOB_filter.csv")