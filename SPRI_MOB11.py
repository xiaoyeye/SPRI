import numpy as np
import pandas as pd
from minepy import MINE
# from sklearn.externals.joblib import Parallel,delayed
import random
import math
from multiprocessing import Pool
# import matplotlib.pyplot as plt
import statsmodels
from statsmodels.stats.multitest import fdrcorrection

def print_stats(mine):
    print ("MIC", mine.mic())
    print ("MAS", mine.mas())
    print ("MEV", mine.mev())
    print ("MCN (eps=0)", mine.mcn(0))
    print ("MCN (eps=1-MIC)", mine.mcn_general())
    print ("GMIC", mine.gmic())
    print ("TIC", mine.tic())

def filter_df(df):
    # process data
    # 1.filtered out genes that are expressed in less than 10% of the array spots 
    # 2.and selected spots with at least ten total read counts
    # the same data preprocess with SPARK

    # gene start from 1
    df_genes = df.columns.values
    df_genes = df_genes[1:]

    # obtain the cell spots
    cellcentroids = df['Unnamed: 0']

    # delete the corridinats
    df.drop(['Unnamed: 0'],axis=1,inplace=True)

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
    
    # filtered out genes that are expressed in less than 10% of the array spots
    threshold_num = 0.1 * df.shape[0]
    delete_list = [] 
    for gene_i in df_genes:
        gene_count = df[gene_i].values
        index = np.nonzero(gene_count)
        count = gene_count[index]
        if len(count)<threshold_num:
            delete_list.append(gene_i)
    delete_list.append('Col_sum')

    df.drop(delete_list,axis=1,inplace=True)
    return(np.array(cell_x),np.array(cell_y),df)

def filter_df_TIC(df):
    # process data
    # only selected spots with at least ten total read counts
    # dose not filter gene

    # gene start from 1
    df_genes = df.columns.values
    df_genes = df_genes[1:]

    # obtain the cell spots
    cellcentroids = df['Unnamed: 0']

    # delete the corridinats
    df.drop(['Unnamed: 0'],axis=1,inplace=True)

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
    df.drop('Col_sum',axis=1,inplace=True)

    return(np.array(cell_x),np.array(cell_y),df)

def background_correction(cell_x,cell_y):
    cell_x = [float(x) for x in cell_x]
    cell_y = [float(x) for x in cell_y]

    # plt.scatter(cell_x,cell_y, edgecolor='b', facecolor='none', alpha=0.5 )
    # plt.xlabel("x"); plt.ylabel("y")
    # plt.show()

    cell_x_copy = cell_x
    cell_y_copy = cell_y

    cell_x = [round(x) for x in cell_x]
    cell_y = [round(x) for x in cell_y]

    xMin = min(cell_x);xMax = max(cell_x);
    yMin = min(cell_y);yMax = max(cell_y);

    cell_x, cell_y = (list(t) for t in zip(*sorted(zip(cell_x,cell_y),key=lambda x: x[0])))
    corridinats = []
    for i in range(len(cell_x)):
        corridinats.append((cell_x[i],cell_y[i]))

    cell_x_noise = []
    cell_y_noise = []
    for x in range(xMin,xMax+1,1):
        for y in range(yMin,yMax+1,1):
            if (x,y) not in corridinats:
                cell_x_noise.append(x)
                cell_y_noise.append(y)

    cell_x_new = cell_x_copy + cell_x_noise
    cell_y_new = cell_y_copy + cell_y_noise

    # plt.scatter(cell_x_new,cell_y_new, edgecolor='b', facecolor='none', alpha=0.5 )
    # plt.xlabel("x"); plt.ylabel("y")
    # plt.show()

    df = pd.DataFrame({'cell_x_noise':cell_x_noise,'cell_y_noise':cell_y_noise})
    return df

def calculate_mic_gene_both(cell_x,cell_y,gene_count):
    x_count_list = []
    y_count_list = []
    index = np.nonzero(gene_count)
    count = gene_count[index]

    x = cell_x[index]
    y = cell_y[index]
    for i in range(len(count)):
        x_count_list.extend([x[i]]*count[i])
        y_count_list.extend([y[i]]*count[i])
            
    assert np.sum(gene_count) == len(x_count_list) == len(y_count_list)

    mine = MINE(alpha=0.6, c=15, est="mic_approx")
    mine.compute_score(x_count_list, y_count_list)
    return(mine.mic(),mine.tic())

def calculate_mic_gene(cell_x,cell_y,gene_count,scoring):
    x_count_list = []
    y_count_list = []
    index = np.nonzero(gene_count)
    count = gene_count[index]

    x = cell_x[index]
    y = cell_y[index]
    for i in range(len(count)):
        x_count_list.extend([x[i]]*count[i])
        y_count_list.extend([y[i]]*count[i])
            
    assert np.sum(gene_count) == len(x_count_list) == len(y_count_list)

    mine = MINE(alpha=0.6, c=15, est="mic_approx")
    mine.compute_score(x_count_list, y_count_list)

    if scoring=="MIC":
        return(mine.mic())
    elif scoring=="TIC":
        return(mine.tic())

def compute_genes_MIC(file_path):
    # preprocess the dataset
    pd_file = pd.read_csv(file_path, sep='\t')
    cell_x,cell_y,pd_file = filter_df_TIC(pd_file)
    cortext_svz_counts = pd_file

    genes = cortext_svz_counts.columns.values
    print('total genes: ',len(genes))
    print('total spots: ',len(cell_x))

    # simulate cell location by poisson
    pd_BG = background_correction(cell_x,cell_y)
    cell_x_BG = pd_BG['cell_x_noise'].values
    cell_y_BG = pd_BG['cell_y_noise'].values
    print(pd_BG)

    # updata the location 
    cell_x_new = np.hstack((cell_x,cell_x_BG)) 
    cell_y_new = np.hstack((cell_y,cell_y_BG)) 

    MIC_list = []
    TIC_list = []
    for gene in genes:
        print('gene:',gene)
        gene_count = cortext_svz_counts[gene].values

        Q1 = np.quantile(gene_count,0.25)
        Q2 = np.quantile(gene_count,0.5)
        Q3 = np.quantile(gene_count,0.75)
        EM = 0.5*Q2 + 0.25*(Q1+Q3)
        mean = int(EM)

        gene_count_noise = mean*np.ones(np.shape(cell_x_BG),dtype=np.int64)

        gene_count_new = np.hstack((gene_count + mean,gene_count_noise)) 
        gene_count_new = gene_count_new.astype(np.int64)

        mic_score,tic_score = calculate_mic_gene_both(cell_x=cell_x_new,cell_y=cell_y_new,gene_count=gene_count_new)
        print('mic:',mic_score)
        print('tic:',tic_score)
        MIC_list.append(mic_score)
        TIC_list.append(tic_score)

    return(MIC_list,TIC_list,genes)


def compute_genes_pvalue(file_path,top_genes,method,num_perm,num_jobs):
    # preprocess the dataset
    pd_file = pd.read_csv(file_path, sep='\t')
    cell_x,cell_y,pd_file = filter_df_TIC(pd_file)
    cortext_svz_counts = pd_file

    # genes = cortext_svz_counts.columns.values
    genes = top_genes
    print('top genes:',len(genes))

    # simulate cell location by poisson
    pd_BG = background_correction(cell_x,cell_y)
    cell_x_noise = pd_BG['cell_x_noise'].values
    cell_y_noise = pd_BG['cell_y_noise'].values

    # updata the location 
    cell_x_new = np.hstack((cell_x,cell_x_noise)) 
    cell_y_new = np.hstack((cell_y,cell_y_noise)) 

    p_values = []
    scores = []
    for gene in genes:
        print('gene:',gene)
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

        score, permutation_scores, pvalue = permutation_test(cell_x_new,cell_y_new,gene,gene_count_new,scoring=method,n_permutations=num_perm,
                    n_jobs=num_jobs, random_state=1)
        print('pvalue:',pvalue)
        print('score:',score)
        p_values.append(pvalue)
        scores.append(score)
    return(p_values,scores,genes)

def permutation_test(cell_x,cell_y,gene,gene_count,scoring="TIC",n_permutations=1000,
                    n_jobs=None, random_state=0,verbose=0):
    score_true = _permutation_test_score(cell_x, cell_y,gene,gene_count,scoring)

    num_processor = n_jobs
    pool = Pool(num_processor)

    result = []
    split_num = n_permutations // num_processor
    left_num = n_permutations%num_processor
    process_list = [split_num] * num_processor
    for i in range(left_num):
        process_list[i] += 1

    # print("split_num:", split_num)
    # print("left_num:", left_num)

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
    pvalue = (np.sum(permutation_scores >= score_true) + 1.0) / (n_permutations + 1)
    return score_true, permutation_scores, pvalue

def _permutation_test_score(cell_x,cell_y,gene,gene_count,scoring):
    """Auxiliary function for permutation_test_score"""
    return calculate_mic_gene(cell_x,cell_y,gene_count,scoring)

def _permutation_test_score_list(cell_x,cell_y,gene,gene_count_list,scoring):
    """Auxiliary function for permutation_test_score"""
    score_list = []
    for gene_count in gene_count_list:
        score_list.append(calculate_mic_gene(cell_x,cell_y,gene_count,scoring))
    return score_list

def compute_FDR(df):
    genes = df['gene'].values
    pvalue = df['pvalue'].values
    tic_list = df['TIC'].values

    # sort
    pvalue, genes,tic_list= (list(t) for t in zip(*sorted(zip(pvalue,genes,tic_list),key=lambda x: x[0])))

    rejected,qvalue = fdrcorrection(np.array(pvalue), alpha=0.05, method='negcorr', is_sorted=True)
    return(genes,pvalue,qvalue,tic_list)


if __name__ == '__main__': 
    # load raw data
    dataname = 'MOB11'
    if dataname == 'MOB11':
        file_path = 'raw_data/Rep11_MOB_count_matrix-1.tsv' # MOB11

    ############################## compute TIC rank 
    MIC_list, TIC_list, genes = compute_genes_MIC(file_path=file_path)
    MIC_list, TIC_list, genes = (list(t) for t in zip(*sorted(zip(MIC_list, TIC_list, genes),key=lambda x: x[1],reverse=True)))
    df = pd.DataFrame({'gene':genes,'mic':MIC_list,'tic':TIC_list})
    save_path = 'results/'+dataname+'/'+dataname+'_tic_rank.csv'
    df.to_csv(save_path,index=False)

    ############################## if perform permutation on top 10% genes:p-val
    thred = 0.1*len(df)
    thred = int(thred)
    top_genes = genes[0:thred]

    p_values, scores, genes = compute_genes_pvalue(file_path=file_path,top_genes = top_genes, method='TIC',num_perm = 10000,num_jobs=40)
    scores_rank, p_values_rank, genes_rank, = (list(t) for t in zip(*sorted(zip(scores,p_values,genes),key=lambda x: x[1])))
    df = pd.DataFrame({'gene':genes_rank,'pvalue':p_values_rank,'TIC':scores_rank})
    
    # fdr correction: q-val
    genes, combined_pvalue, pvalue_corrected, tic_list = compute_FDR(df)
    genes, combined_pvalue, pvalue_corrected, tic_list = (list(t) for t in zip(*sorted(zip(genes, combined_pvalue,pvalue_corrected,tic_list),key=lambda x: float(x[2]))))
    df = pd.DataFrame({'gene':genes,'pvalue':combined_pvalue,'qvalue':pvalue_corrected,'TIC':tic_list})
    save_path = 'results/'+dataname+'/'+dataname+'_pvals.csv'
    df.to_csv(save_path)
