#得到每个样本下物种的cnt个数
import sys
import pandas as pd
import numpy as np
import os
from collections import Counter

env = sys.argv[1]
path_home = sys.argv[2]
path_data = os.path.join(path_home,"data")
# srr_strat = int(sys.argv[2])
# srr_end  =  int(sys.argv[3])
# nohup python 8_strian_cnt_and_abu.py test 0 1


#存放sp read的kmer,用于计算 strain 的个数
path_read_kmer_root = os.path.join(path_data,"env/data/3_sp_read_kmer_sorted",env)

# path_strian_kmer_root = os.path.join("/data/huixingqi/data/env/data/5_strain_unikmer",env,"kmer")

#存放所有srr下strain和read kmer的交集
path_kmer_root = os.path.join(path_data,"env/data/6_read_strain_inter",env+"_new")
#保存每个srr的sp的strian个数
path_cnt_root = os.path.join(path_data,"env/data/8_strain_cnt_and_abu",env,"cnt")
if not os.path.exists(path_cnt_root):
    os.makedirs(path_cnt_root)
#保存每个srr的sp的strian的丰度
path_abu_root = os.path.join(path_data,"env/data/8_strain_cnt_and_abu",env,"abu")
if not os.path.exists(path_abu_root):
    os.makedirs(path_abu_root)

list_srr = sorted(os.listdir(path_kmer_root)) #这个是strain unikmer 和 sp kmer的交集
# 所以我们要保留一个非uni的kmer

def get_k50(list_cnt):
    list_cnt.sort() 
    #先排序吧，不然可能unique_labels的顺序是kmer数第一次出现的顺序，而不是从小到大
    #如果不是从小到大，那么就不是N50了
    unique_labels,counts = np.unique(list_cnt, return_counts=True)
    sorted_indices = np.argsort(unique_labels)
    sorted_labels = unique_labels[sorted_indices]
    sorted_counts = counts[sorted_indices]
    # most_common_cluster_idx = np.argmax(counts)#出现次数最多的kmer的下标
    sum_lc_all = sum(list_cnt)#kmer的总数 
    sum_lc=0#用于计算n50
    for l,c in zip(sorted_labels,sorted_counts):
        sum_lc+=l*c
        if sum_lc>sum_lc_all*0.5:
            l50 = l
            break
    return l50


for srr in list_srr:
    # print(srr)
    
    df_srr_cnt = pd.DataFrame(columns=["sp","cnt_strain"])
    df_srr_abu = pd.DataFrame(columns=["sp","strain","abu"])
    
    # 6_read_strain_inter/SRR111
    path_srr = os.path.join(path_kmer_root,srr) 
    list_sp = sorted(os.listdir(path_srr)) 
    
    #记录每个sp的名字（不是所有sp都满足条件） 
    # list_sp_name =[] 
    for sp in list_sp:
        
        list_cnt_no_ref2 = []#一个sp一个
        path_sp = os.path.join(path_srr,sp)
        list_fa = sorted(os.listdir(path_sp))
        max_noref = 0 #记录最多的noref
        if len(list_fa)!=0:#这个物种有fa
            
            path_sp_kmer = os.path.join(path_read_kmer_root,srr,sp+".fa")
            df_sp_kmer = pd.read_csv(path_sp_kmer,header=None,sep="\t")
            df_sp_kmer.columns=["kmer","cnt_read"]
            list_read_cnt = df_sp_kmer["cnt_read"].to_list()
            # k50_read = get_k50(list_read_cnt[:])
            k50_read = get_k50(np.unique(list_read_cnt))
            cnt_sp = np.max(list_read_cnt) / k50_read
            
            float_sp = cnt_sp - int(cnt_sp)
            if float_sp >= 0.1:
                 cnt_sp =  int(cnt_sp) + 1
            else:
                cnt_sp =  int(cnt_sp)


            # list_sp_name.append(sp)
            cnt_fa = 0 #记录数据库中strain个数
            
            norefs = [] # 记录每个strain计算得到的noref个数
            abu_norefs = []
            list_max_cnt_fa =[]
            for fa in list_fa: #遍历这个物种的strain
                path_fa = os.path.join(path_sp,fa)
                if os.path.getsize(path_fa)!=0:
                    df_kmer = pd.read_csv(path_fa,header=None,sep="\t") 
                    #把strain的unikmer和sp的kmer合并起来,得到的结果
                    df_kmer.columns=["kmer","cnt"] #strain的kmer的数据
                    list_cnt = df_kmer["cnt"].tolist()
                    mean_cnt = np.mean(list_cnt)
                    max_cnt = np.max(list_cnt)
                    # strain_cnt = max_cnt/mean_cnt
                    # k50_strain = get_k50((list_cnt).copy())
                    k50_strain = get_k50(np.unique(list_cnt))
                    strain_cnt = max_cnt / k50_strain
                    
                    float_cnt = strain_cnt - int(strain_cnt)
                    if float_cnt >= 0.1:
                        noref_cnt =  int(strain_cnt)
                    else:
                        noref_cnt =  int(strain_cnt) - 1
                    if noref_cnt > 0:
                        norefs.append(noref_cnt)
                    if noref_cnt==1:
                        abu_norefs.append(max_cnt-mean_cnt)
                    cnt_fa+=1


                    values = {"sp":[sp],"strain":[fa],"abu":[mean_cnt]}
                    df_srr_abu = pd.concat([df_srr_abu,pd.DataFrame(values)],ignore_index=True,axis=0)
                    

            if len(norefs)>0: #有noref
                min_strain = max(norefs)
                max_strain = sum(norefs)
                cnt_noref = cnt_sp - cnt_fa
                if cnt_noref <min_strain:
                    cnt_noref = min_strain
                elif cnt_noref > max_strain:
                    cnt_noref = max_strain
                
                cnt_sp = cnt_noref + cnt_fa
                
                
                
            else:#没有noref
                cnt_sp = cnt_fa
                cnt_noref = 0
                
            if cnt_noref == 1:
                
                values = {"sp":[sp],"strain":["noref"],"abu":[np.mean(abu_norefs)]}
                df_srr_abu = pd.concat([df_srr_abu,pd.DataFrame(values)],ignore_index=True,axis=0)
            # cnt_fa2 = cnt_fa+cnt_noref
            values_cnt= {"sp":[sp],"cnt_strain":cnt_sp}
            df_srr_cnt = pd.concat([df_srr_cnt,pd.DataFrame(values_cnt)],ignore_index=True,axis=0)
  
    path_srr_cnt = os.path.join(path_cnt_root,srr+".txt") 
    df_srr_cnt.to_csv(path_srr_cnt,header=None,sep="\t",index=None)
    path_srr_abu = os.path.join(path_abu_root,srr+".txt")
    df_srr_abu.to_csv(path_srr_abu,header=None,sep="\t",index=None)

#统计结果并保存
#把上面按照样本生成的文件整合在一起
#一个是当作测试集，一个是训练集，一个是全部数据
path_all = os.path.join(path_data,"env/data/8_strain_cnt_and_abu",env,"all")
if not os.path.exists(path_all):
    os.makedirs(path_all)

list_srr = sorted(os.listdir(path_cnt_root))

#cnt
#训练集

df_all_train = pd.DataFrame(columns=["sp","cnt"])
for srr in list_srr[30:]: #遍历srr
    path_srr = os.path.join(path_cnt_root,srr)
    df_srr = pd.read_csv(path_srr,header=None,sep="\t")
    df_srr.columns=["sp","cnt"]
    df_all_train = pd.concat([df_all_train,df_srr],ignore_index=True,axis=0)
    
path_all_cnt_train = os.path.join(path_all,"sp_cnt_train.csv")
df_all_train.to_csv(path_all_cnt_train,header=None,index=None,sep="\t")

#测试集
df_all_test = pd.DataFrame(columns=["sp","cnt"])
for srr in list_srr[:30]: #遍历srr
    path_srr = os.path.join(path_cnt_root,srr)
    df_srr = pd.read_csv(path_srr,header=None,sep="\t")
    df_srr.columns=["sp","cnt"]
    df_all_test = pd.concat([df_all_test,df_srr],ignore_index=True,axis=0)
    
path_all_cnt_test = os.path.join(path_all,"sp_cnt_test.csv")
df_all_test.to_csv(path_all_cnt_test,header=None,index=None,sep="\t")
#all
df_all = pd.concat([df_all_test,df_all_train],ignore_index=True,axis=0)
path_all_cnt_all = os.path.join(path_all,"sp_cnt_all.csv")
# print(df_all.columns)
df_all.to_csv(path_all_cnt_all,header=None,index=None,sep="\t")


#abu
#训练集
df_all_train = pd.DataFrame(columns=["sp","abu"])
for srr in list_srr[30:]: #遍历srr 
    path_srr = os.path.join(path_abu_root,srr)
    df_srr = pd.read_csv(path_srr,header=None,sep="\t")
    df_srr.columns=["sp","strain","abu"]
    df_group = df_srr.groupby("sp")
    for key,index in df_group.groups.items():
        list_abu = list(df_srr.iloc[index]["abu"])
        # list_abu_std = [ i/sum(list_abu) for i in list_abu][:-1] 这个是没有考虑noref的
        list_abu_std = [ i/sum(list_abu) for i in list_abu] #因为现在只有strain为1的noref被保留，所以noref也可以考虑
        df_list_abu = pd.DataFrame(list_abu_std, columns=["abu"])
        df_list_abu["sp"] = key
        df_all_train = pd.concat([df_all_train,df_list_abu],ignore_index=True,axis=0)
path_all_cnt_train = os.path.join(path_all,"sp_abu_train.csv")
df_all_train.to_csv(path_all_cnt_train,header=None,index=None,sep="\t")

#测试集
df_all_test = pd.DataFrame(columns=["sp","abu"])
for srr in list_srr[:30]: #遍历srr 
    path_srr = os.path.join(path_abu_root,srr)
    df_srr = pd.read_csv(path_srr,header=None,sep="\t")
    df_srr.columns=["sp","strain","abu"]
    df_group = df_srr.groupby("sp")
    for key,index in df_group.groups.items():
        list_abu = list(df_srr.iloc[index]["abu"])
        # list_abu_std = [ i/sum(list_abu) for i in list_abu][:-1]
        list_abu_std = [ i/sum(list_abu) for i in list_abu]
        df_list_abu = pd.DataFrame(list_abu_std, columns=["abu"])
        df_list_abu["sp"] = key
        df_all_test = pd.concat([df_all_test,df_list_abu],ignore_index=True,axis=0)
path_all_cnt_test = os.path.join(path_all,"sp_abu_test.csv")
df_all_test.to_csv(path_all_cnt_test,header=None,index=None,sep="\t")

#all
df_all = pd.concat([df_all_test,df_all_train],ignore_index=True,axis=0)
path_all_cnt_all = os.path.join(path_all,"sp_abu_all.csv")
df_all.to_csv(path_all_cnt_all,header=None,index=None,sep="\t")


