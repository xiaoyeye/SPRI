import numpy as np
import pandas as pd

import NaiveDE
import SpatialDE
import time


def get_coords(index):
    coords = pd.DataFrame(index=index)
    coords['x'] = index.str.split('x').str.get(0).map(float)
    coords['y'] = index.str.split('x').str.get(1).map(float)
    return coords

# df = pd.read_table('raw_data/Rep11_MOB_count_matrix-1.tsv', index_col=0)
df = pd.read_table('raw_data/Rep12_MOB_count_matrix-1.tsv', index_col=0)

df = df.T[df.sum(0) >= 3].T  
sample_info = get_coords(df.index)
sample_info['total_counts'] = df.sum(1)
sample_info = sample_info.query('total_counts > 10')  
df = df.loc[sample_info.index]

X = sample_info[['x', 'y']]

start_time = time.time()
dfm = NaiveDE.stabilize(df.T).T
res = NaiveDE.regress_out(sample_info, dfm.T, 'np.log(total_counts)').T
res['log_total_count'] = np.log(sample_info['total_counts'])
results = SpatialDE.run(X, res)
elapsed_time = time.time() - start_time

de_results = results[(results.qval < 0.05)].copy()
ms_results = SpatialDE.model_search(X, res, de_results)


# results.to_csv('output/Rep11_MOB_spe.csv')
# ms_results.to_csv('output/Rep11_MOB_ms_spe.csv')

results.to_csv('output/Rep12_MOB_spe.csv')
ms_results.to_csv('output/Rep12_MOB_ms_spe.csv')



