import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
from collections import Counter
import os

def filter_df(df):
    # process data
    # 1.filtered out genes that are expressed in less than 10% of the array spots 
    # 2.and selected spots with at least ten total read counts

    # gene start from 1
    df_genes = df.columns.values
    df_genes = df_genes[1:]

    # obtain the cell spots
    cellcentroids = df['index']

    # delete the corridinats
    df.drop(['index'],axis=1,inplace=True)

    # compute the read counts of every spot
    df['Col_sum'] = df.apply(lambda x: x.sum(), axis=1)
    
    # obtain the cell spots after delete the spots less than ten total read counts
    cellcentroids = cellcentroids.drop(df[df['Col_sum']<=10].index)
    cellcentroids = cellcentroids.values
    cell_x = []
    cell_y = []
    for cell_i in cellcentroids:
        cell_x.append(cell_i.split('x')[0])
        cell_y.append(cell_i.split('x')[1])

    # spots less than ten total read counts
    df = df.drop(df[df['Col_sum']<=10].index)
    
    ########## # filtered out genes that are expressed in less than 10% of the array spots
    # threshold_num = 0.1 * df.shape[0]
    # delete_list = [] 
    # for gene_i in df_genes:
    #     gene_count = df[gene_i].values
    #     index = np.nonzero(gene_count)
    #     count = gene_count[index]
    #     if len(count)<threshold_num:
    #         delete_list.append(gene_i)
    # delete_list.append('Col_sum')
    # df.drop(delete_list,axis=1,inplace=True)

    ########## only filtered out col_sum
    df.drop(['Col_sum'],axis=1,inplace=True)

    return(np.array(cell_x),np.array(cell_y),df)


def prefilter_specialgenes(adata,Gene1Pattern="ERCC",Gene2Pattern="MT-"):
    id_tmp1=np.asarray([not str(name).startswith(Gene1Pattern) for name in adata.var_names],dtype=bool)
    id_tmp2=np.asarray([not str(name).startswith(Gene2Pattern) for name in adata.var_names],dtype=bool)
    id_tmp=np.logical_and(id_tmp1,id_tmp2)
    adata._inplace_subset_var(id_tmp)

def prefilter_genes(adata,min_counts=None,max_counts=None,min_cells=10,max_cells=None):
    if min_cells is None and min_counts is None and max_cells is None and max_counts is None:
        raise ValueError('Provide one of min_counts, min_genes, max_counts or max_genes.')
    id_tmp=np.asarray([True]*adata.shape[1],dtype=bool)
    id_tmp=np.logical_and(id_tmp,sc.pp.filter_genes(adata.X,min_cells=min_cells)[0]) if min_cells is not None  else id_tmp
    id_tmp=np.logical_and(id_tmp,sc.pp.filter_genes(adata.X,max_cells=max_cells)[0]) if max_cells is not None  else id_tmp
    id_tmp=np.logical_and(id_tmp,sc.pp.filter_genes(adata.X,min_counts=min_counts)[0]) if min_counts is not None  else id_tmp
    id_tmp=np.logical_and(id_tmp,sc.pp.filter_genes(adata.X,max_counts=max_counts)[0]) if max_counts is not None  else id_tmp
    adata._inplace_subset_var(id_tmp)


def background_correction(cell_x,cell_y):

    cell_x = [int(x) for x in cell_x]
    cell_y = [int(x) for x in cell_y]
    print(len(cell_x)) # 2264
    print(len(cell_y)) # 2264

    plt.scatter(cell_x,cell_y, edgecolor='b', facecolor='none', alpha=0.5,s=0.1)
    plt.xlabel("x"); plt.ylabel("y")
    # plt.show()
    plt.savefig(save_prefix + 'original_spot.png')
    plt.close()

    cell_x_copy = cell_x
    cell_y_copy = cell_y

    xMin = min(cell_x);xMax = max(cell_x)
    yMin = min(cell_y);yMax = max(cell_y)
    print('xMin:',xMin) #8278
    print('xMax:',xMax) #23669
    print('yMin:',yMin) #6306
    print('yMax:',yMax) #20688

    df_location = pd.DataFrame({'cell_x':cell_x,'cell_y':cell_y})
    df_location.to_csv(save_prefix + dataname + '_cell_location.csv')

    cell_x_BG = []
    cell_y_BG = []

    cell_x_select = []
    cell_y_select = []

    all_point = []
    for x, y in zip(cell_x, cell_y):
        all_point.append((x,y))
    # sorted_all_point = sorted(all_point, key=lambda x: x[0])
    sorted_all_point = sorted(all_point, key=lambda x: (x[0],x[1]))

    y_thred_list = []
    prev_point = sorted_all_point[0]
    for point in sorted_all_point[1:]:
        y_thred_list.append(point[1]-prev_point[1])
        prev_point = point
    collection_words = Counter(y_thred_list)
    y_thred = max(collection_words,key=collection_words.get)+1

    former_point = sorted_all_point[0]
    valleys = []
    valleys.append(former_point[0])
    for point in sorted_all_point[1:]:
        if point[0] != former_point[0] and point[1] < former_point[1]:
            valleys.append(point[0])
        former_point = point

    print("valleys: ", valleys) 
    print(len(valleys))

    if dataname == '151674':
        valleys.append(sorted_all_point[-1][0])

    peaks=valleys[1:] + [xMax+1]
    print("peaks: ", peaks)

    turn = 0 
    for start, end in zip(valleys, peaks):
        turn += 1
        #if turn < 76:
        #    continue
        one_colum_x = []
        one_colum_y = []

        df_select = df_location[df_location['cell_x']<end]
        df_select = df_select[df_select['cell_x']>=start]
        df_select = df_select.sort_values(by='cell_x')
        df_select = df_select.reset_index()
        cell_x_select = df_select['cell_x'].values.tolist()
        cell_y_select = df_select['cell_y'].values.tolist()

        dic = {}
        point_list = []
        for k,v in zip(cell_x_select,cell_y_select):
            dic[k] = v
            point_list.append((k,v))
        
        min_x = min(cell_x_select)
        max_x = max(cell_x_select)
        min_y = min(cell_y_select) 
        max_y = max(cell_y_select) 

        sorted_point_list = sorted(point_list, key=lambda x: x[1])
        print("sorted_point_list: ", sorted_point_list)
        print("len(sorted_point_list): ", len(sorted_point_list))
        
        his_cnt = len(cell_x_BG)
        # small
        temp_x = min_x - 1
        temp_y = min_y - y_thred
        while temp_y >= yMin:
            cell_x_BG.append(temp_x)
            cell_y_BG.append(temp_y)
            temp_x -= 1
            temp_y -= y_thred


        s_cnt = len(cell_x_BG) - his_cnt
        print("small add: ", s_cnt)
        his_cnt = len(cell_x_BG)
        
        former_x = sorted_point_list[0][0]
        former_y = sorted_point_list[0][1]
        for point in sorted_point_list[1:]:
            cur_x = point[0]
            cur_y = point[1]
            if cur_y - former_y > y_thred:
                point_num = int((cur_y - former_y) / y_thred)
                print("current add point_num: ", cur_y, former_y, point_num)
                for _ in range(0, point_num):
                    cell_x_BG.append(former_x + 1)
                    cell_y_BG.append(former_y + y_thred)
                    former_x = former_x + 1
                    former_y = former_y + y_thred
            former_x = cur_x
            former_y = cur_y
        m_cnt = len(cell_x_BG) - his_cnt
        print("middle add: ", m_cnt)
        his_cnt = len(cell_x_BG)
        
        # large
        temp_x = max_x + 1
        temp_y = max_y + y_thred
        while temp_y <= yMax:
            cell_x_BG.append(temp_x)
            cell_y_BG.append(temp_y)
            temp_x += 1
            temp_y += y_thred

        l_cnt = len(cell_x_BG) - his_cnt
        print("large add: ", l_cnt)
        #break
        #if turn == 2:
        #    break

    cell_x_new = cell_x_copy + cell_x_BG
    cell_y_new = cell_y_copy + cell_y_BG

    print(len(cell_x_BG))
    print(len(cell_y_BG))

    cell_x_noise = cell_x_BG 
    cell_y_noise = cell_y_BG 
    df = pd.DataFrame({'cell_x_noise':cell_x_noise,'cell_y_noise':cell_y_noise})
    df.to_csv(save_prefix + dataname + '_background_spot.csv')

    #Plotting
    plt.scatter(cell_x_new,cell_y_new, edgecolor='b', facecolor='none', alpha=0.5,s=0.1)
    plt.xlabel("x"); plt.ylabel("y")
    plt.savefig(save_prefix + 'final_spot.png')
    plt.close()

    plt.scatter(cell_x_BG, cell_y_BG, edgecolor='b', facecolor='none', alpha=0.5,s=0.1)
    plt.xlabel("x"); plt.ylabel("y")
    plt.savefig(save_prefix + 'background_spot.png')
    plt.close()



from scanpy import read_10x_h5

#Read data : 
dataname = '151507' #151507, 151508, 151509, 151510, 151669, 151670, 151671, 151672, 151673, 151674, 151675, 151676

# load gene expression: expressiion_matrix
save_prefix = 'processed_data/' + dataname 

if os.path.exists(save_prefix):
    print(save_prefix)
else:
    os.mkdir(save_prefix)

save_prefix = save_prefix + '/'

adata = read_10x_h5('raw_data/LIBD/' + dataname + '/'+dataname+'_filtered_feature_bc_matrix.h5') #151673_filtered_feature_bc_matrix.h5 is expression_matrix.h5
spatial=pd.read_csv('raw_data/LIBD/'+ dataname + '/tissue_positions_list.txt',sep=",",header=None,na_filter=False,index_col=0) 
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

x_array=adata.obs["x_array"].tolist() #3639
y_array=adata.obs["y_array"].tolist()
x_pixel=adata.obs["x_pixel"].tolist() #3639
y_pixel=adata.obs["y_pixel"].tolist() #3639

# generate cell location
location_list = []
for i in range(len(x_pixel)):
    location = str(x_pixel[i]) + 'x' + str(y_pixel[i])
    location_list.append(location)

# preprocess
adata.var_names_make_unique() #(3639, 33538) -> (3639, 33538)
# prefilter_genes(adata,min_cells=3) # avoiding all genes are zeros: (3639, 33538) -> (3639, 19151)
# prefilter_specialgenes(adata) #util.py

cell_names = adata.obs_names.tolist() # 3639
gene_names = adata.var["genename"].tolist() #19130
expressiion_matrix = adata.X.A  #(3639, 33538) -> (3639, 19130) after data preprocess
# expressiion_matrix = adata.X
print('cell_names: ',len(cell_names))
print('gene_names: ',len(gene_names))
print('x_pixel: ',len(x_pixel))
print('y_pixel: ',len(y_pixel))
print('expressiion_matrix: ',np.shape(expressiion_matrix))

# generate 
df_expression = pd.DataFrame(expressiion_matrix)
df_expression.index = location_list
df_expression.columns = gene_names
df_expression = df_expression.reset_index() #(2264, 19465)
df_expression.to_csv(save_prefix + dataname + '_count_matrix.csv')
print(df_expression)

# data pre-process
cell_x,cell_y,pd_file = filter_df(df_expression) # (2264, 19465)
print(pd_file)

background_correction(cell_x,cell_y)
