import numpy as np
from minepy import MINE
import pandas as pd
from sklearn.externals.joblib import Parallel,delayed
import random
#import matplotlib.pyplot as plt
from statsmodels.stats.multitest import fdrcorrection
from multiprocessing import Pool
import math
import os 
import time

def print_stats(mine):
    print ("MIC", mine.mic())
    print ("MAS", mine.mas())
    print ("MEV", mine.mev())
    print ("MCN (eps=0)", mine.mcn(0))
    print ("MCN (eps=1-MIC)", mine.mcn_general())
    print ("GMIC", mine.gmic())
    print ("TIC", mine.tic())

def compute_FDR(pvalue,genes,score):
    # genes = df['gene'].values
    # combined_pvalue = df['combined_pvalue'].values

    # sort
    pvalue, genes,score = (list(t) for t in zip(*sorted(zip(pvalue,genes,score),key=lambda x: x[0])))
    print(pvalue[0:10])

    rejected,qvalue = fdrcorrection(np.array(pvalue), alpha=0.05, method='negcorr', is_sorted=True)
    return(genes,pvalue,qvalue,score)

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
    df.drop(['Col_sum','Unnamed: 0'],axis=1,inplace=True)

    return(np.array(cell_x),np.array(cell_y),df)

# def calculate_mic_gene(cell_x,cell_y,gene,gene_count,scoring):
def calculate_mic_gene(cell_x,cell_y,gene_count,scoring):
    x_count_list = []
    y_count_list = []
    # print(gene)
    # gene_count = cortext_svz_counts[gene].values
    index = np.nonzero(gene_count)
    count = gene_count[index]

    x = cell_x[index]
    y = cell_y[index]
    # count the number of x, y
    # for i in range(len(count)):
    #     for j in range(count[i]):
    #         x_count_list.append(x[i])
    #         y_count_list.append(y[i])
    for i in range(len(count)):
        x_count_list.extend([x[i]]*count[i])
        y_count_list.extend([y[i]]*count[i])
            
    assert np.sum(gene_count) == len(x_count_list) == len(y_count_list)

    mine = MINE(alpha=0.6, c=15, est="mic_approx")
    mine.compute_score(x_count_list, y_count_list)
    # print(mine.mic())

    if scoring=="MIC":
        return(mine.mic())
    elif scoring=="TIC":
        return(mine.tic())

def calculate_mic_gene_both(cell_x,cell_y,gene_count,scoring):
    x_count_list = []
    y_count_list = []
    # print(gene)
    # gene_count = cortext_svz_counts[gene].values
    index = np.nonzero(gene_count)
    count = gene_count[index]

    x = cell_x[index]
    y = cell_y[index]
    # count the number of x, y
    # for i in range(len(count)):
    #     for j in range(count[i]):
    #         x_count_list.append(x[i])
    #         y_count_list.append(y[i])
    for i in range(len(count)):
        x_count_list.extend([x[i]]*count[i])
        y_count_list.extend([y[i]]*count[i])
            
    assert np.sum(gene_count) == len(x_count_list) == len(y_count_list)

    start = time.time()
    # mine = MINE(alpha=0.6, c=15, est="mic_approx")
    mine = MINE(alpha=0.6, c=15, est="mic_e")
    mine.compute_score(x_count_list, y_count_list)
    end = time.time()
    print('time: ', end-start)
    print('len(x_count_list): ',len(x_count_list))
    print('len(y_count_list): ',len(y_count_list))
    return(mine.mic(),mine.tic())


def build_dict():
    file_path = 'results/permutation/MOB/MOB_EM_noise_qvalue.csv'
    pd_file = pd.read_csv(file_path)
    print(pd_file)
    genes = pd_file['gene'].values
    pvalue = pd_file['pvalue'].values
    qvalue = pd_file['qvalue'].values
    print(genes)
    dict_pvalue = {}
    dict_qvalue = {}
    for i in range(len(genes)):
        dict_pvalue[genes[i]]=pvalue[i]
        dict_qvalue[genes[i]]=qvalue[i]
    return(dict_pvalue,dict_qvalue)


def compute_genes_MIC(file_path,method,num_perm,num_jobs):
    # preprocess the dataset
    # pd_file = pd.read_csv(file_path, sep='\t')
    pd_file = pd.read_csv(file_path)
    print(pd_file)
    cell_x,cell_y,pd_file = filter_df(pd_file)
    print(pd_file) #[3639x33539]
    print('cell_x',len(cell_x)) #3639
    print('cell_y',len(cell_y)) #3639
    cortext_svz_counts = pd_file

    dataname = file_path.split('/')[1]
    print(dataname)

    genes = cortext_svz_counts.columns.values
    # df_new = pd.read_csv('results/permutation/BC/BC_SPARK_qvalue.csv')
    # genes = df_new['gene'].values
    print('total genes:',len(genes)) #33539
    # print(genes[0:10])

    # dict_pvalue,dict_qvalue = build_dict()

    # simulate cell location by poisson
    pd_noise = pd.read_csv('processed_data/'+dataname+'/'+dataname+'_background_spot.csv')
    cell_x_noise = pd_noise['cell_x_noise'].values
    cell_y_noise = pd_noise['cell_y_noise'].values
    print('cell_x_noise',len(cell_x_noise))

    # updata the location 
    cell_x_new = np.hstack((cell_x,cell_x_noise)) 
    cell_y_new = np.hstack((cell_y,cell_y_noise)) 
    print('cell_x_new',len(cell_x_new))

    ##### later: we can compute for ten times and obtain the highest result
    # for isid in range(1:10):
    # set.seed(isid)
    MIC_list = []
    TIC_list = []
    # pvalue_list = []
    # qvalue_list = []
    for gene in genes:
        print('gene:',gene)
        # score, permutation_scores, pvalue = permutation_test(cell_x,cell_y,gene,scoring="MIC",n_permutations=10,
        #             n_jobs=None, random_state=0)
        gene_count = cortext_svz_counts[gene].values

        Q1 = np.quantile(gene_count,0.25)
        Q2 = np.quantile(gene_count,0.5)
        Q3 = np.quantile(gene_count,0.75)
        EM = 0.5*Q2 + 0.25*(Q1+Q3)
        mean = int(EM)
        print('mean',mean)

        gene_count_noise = mean*np.ones(np.shape(cell_x_noise),dtype=np.int64)

        gene_count_new = np.hstack((gene_count + mean,gene_count_noise)) 
        gene_count_new = gene_count_new.astype(np.int64)

        mic_score,tic_score = calculate_mic_gene_both(cell_x=cell_x_new,cell_y=cell_y_new,gene_count=gene_count_new,scoring=method)
        print('mic:',mic_score)
        print('tic:',tic_score)
        MIC_list.append(mic_score)
        TIC_list.append(tic_score)
        # pvalue_list.append(dict_pvalue[gene])
        # qvalue_list.append(dict_qvalue[gene])

    # return(MIC_list,TIC_list,genes,pvalue_list,qvalue_list)
    return(MIC_list,TIC_list,genes)


# generates a null distribution by calculating the MIC on 1000 different permutations of the one gene, 
# where counts remain the same but cordinate(x,y) undergo different permutations.
def permutation_test(cell_x,cell_y,gene,gene_count,scoring="TIC",n_permutations=1000,
                    n_jobs=None, random_state=0,verbose=0):
    # print('n_permutations',n_permutations)

# def permutation_test_score(estimator, X, y, *, groups=None, cv=None,
#                            n_permutations=100, n_jobs=None, random_state=0,
#                            verbose=0, scoring=None, fit_params=None):

    # random_state = np.random.mtrand._rand
    #rng = np.random.RandomState(random_state)

    # We clone the estimator to make sure that all the folds are
    # independent, and that it is pickle-able.
    score_true = _permutation_test_score(cell_x, cell_y,gene,gene_count,scoring)
    '''
    permutation_scores = Parallel(n_jobs=n_jobs, verbose=verbose)(
        delayed(_permutation_test_score)(
             cell_x, cell_y,gene, _shuffle(gene_count,rng),scoring)
        for _ in range(n_permutations))
    '''
    num_processor = n_jobs
    pool = Pool(num_processor)

    result = []
    split_num = n_permutations // num_processor
    left_num = n_permutations%num_processor
    process_list = [split_num] * num_processor
    for i in range(left_num):
        process_list[i] += 1

    print("split_num:", split_num)
    print("left_num:", left_num)

    def per_process_shuffle(a, repeat_num, i): 
        rng = np.random.RandomState(i)
        res = []
        for _ in range(repeat_num):
            rng.shuffle(a)
            copy_a = a.copy()
            res.append(copy_a)   
        return res

    for i in range(num_processor):
        result.append(
            pool.apply_async(
                _permutation_test_score_list, args=(cell_x, cell_y,gene, np.array(per_process_shuffle(gene_count.tolist(), process_list[i], i)),scoring)
            )
        )
    permutation_scores = []
    for i  in result:
        for score in i.get():
            permutation_scores.append(score)
    
    pool.close()
    pool.join()

    permutation_scores = np.array(permutation_scores)
    save_path = 'results/permutation_score/BC/'+gene+'.npy'
    np.save(save_path, permutation_scores)
    pvalue = (np.sum(permutation_scores >= score_true) + 1.0) / (n_permutations + 1)
    return score_true, permutation_scores, pvalue


def _permutation_test_score(cell_x,cell_y,gene,gene_count,scoring):
    """Auxiliary function for permutation_test_score"""
    # Adjust length of sample weights
    # cell_x = corrdinate[0]
    # cell_y = corrdinate[1]
    # print(gene_count)
    # return calculate_mic_gene(cell_x,cell_y,gene,gene_count,scoring)
    return calculate_mic_gene(cell_x,cell_y,gene_count,scoring)

def _permutation_test_score_list(cell_x,cell_y,gene,gene_count_list,scoring):
    """Auxiliary function for permutation_test_score"""
    # Adjust length of sample weights
    # cell_x = corrdinate[0]
    # cell_y = corrdinate[1]
    # print(gene_count)
    score_list = []
    for gene_count in gene_count_list:
        # score_list.append(calculate_mic_gene(cell_x,cell_y,gene,gene_count,scoring))
        score_list.append(calculate_mic_gene(cell_x,cell_y,gene_count,scoring))
    return score_list


def _shuffle(gene_count,rng):
    # 'exact': all possible permutation is considered
    # 'approximate' : the number of drawn samples is given by 'num_rounds'
    """Return a shuffled copy of corrdinate"""
    # the same null-dist is used for all genes for comparibility???

    # indices = random_state.permutation(len(y))
    # return _safe_indexing(y, indices)

    # indices = np.random.permutation(len(cell_x))
    # cell_x = cell_x[indices]
    # cell_y = cell_y[indices]

    rng.shuffle(gene_count)

    return gene_count


if __name__ == '__main__': 
    # load raw data
    dataname = '151507' #'151507','151676'
    file_path = 'processed_data/'+ dataname + '/'+ dataname + '_count_matrix.csv'

    MIC_list,TIC_list,genes = compute_genes_MIC(file_path=file_path,method='TIC',num_perm = 10000,num_jobs=40)
    MIC_list, TIC_list, genes = (list(t) for t in zip(*sorted(zip(MIC_list, TIC_list, genes),key=lambda x: x[1],reverse=True)))
    df = pd.DataFrame({'gene':genes,'mic':MIC_list,'tic':TIC_list})

    save_prefix = 'results/'+ dataname 
    if os.path.exists(save_prefix):
        print(save_prefix)
    else:
        os.mkdir(save_prefix)

    save_path = 'results/'+ dataname + '/' + dataname + '_mic_tic.csv'
    df.to_csv(save_path)



