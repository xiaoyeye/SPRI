library(Giotto)

##### BC2020
path_to_matrix = 'processed_data/BC2020/BC2020_count.txt'
path_to_locations = 'processed_data/BC2020/BC2020_loc.txt'
save_path = "output/BC2020/BC2020_kmeans.csv"


# 2. use an existing matrix and data.table
expression_matrix = readExprMatrix(path_to_matrix) # fast method to read expression matrix
cell_locations = data.table::fread(path_to_locations)

print(cell_locations)

my_giotto_object = createGiottoObject(raw_exprs = expression_matrix,
                                      spatial_locs = cell_locations)

# processing
my_giotto_object <- filterGiotto(gobject = my_giotto_object,
                              expression_threshold = 1,
                              gene_det_in_min_cells = 1,
                              min_det_genes_per_cell = 1,
                              expression_values = c('raw'))
my_giotto_object <- normalizeGiotto(gobject = my_giotto_object, scalefactor = 6000)


visium_brain <- createSpatialNetwork(gobject = my_giotto_object, 
                                     method = 'kNN', k = 5, 
                                     maximum_distance_knn = 400, 
                                     name = 'spatial_network')
km_spatialgenes = binSpect(visium_brain, calc_hub = T, hub_min_int = 5,
                  spatial_network_name = 'spatial_network',bin_method = 'kmeans')


print(km_spatialgenes)
write.csv(km_spatialgenes,save_path)
