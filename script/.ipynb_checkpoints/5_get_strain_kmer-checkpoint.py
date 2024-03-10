import pandas as pd
import os
from Bio import SeqIO
import pickle
import shutil
import subprocess
from collections import Counter
import sys

env = sys.argv[1]
start = int(sys.argv[2])
end = int(sys.argv[3])
batch_name = str(start)+"_"+str(end)
k_kmer = int(sys.argv[4])
path_home = sys.argv[5]
path_data = os.path.join(path_home,"data")
path_script = os.path.join(path_home,"script")
# k_kmer = 14
# k_kmer = 35
# k_kmer = env[-2:]


#存放每个物种strain的地址，执行路径
path_kmer_root = os.path.join(path_data,"env/data/4_sp_strain_path_sample",env)
list_sp = sorted(os.listdir(path_kmer_root))
print(list_sp[start:end])

# #存放临时的jf结果，输出结果路径
path_jf = os.path.join(path_data,"env/data/5_strain_unikmer",env,"sp")
# if not os.path.exists(path_jf):
#     os.makedirs(path_jf)

# #kmer相交后得到unikmer保存下来        
path_jf_unikmer = os.path.join(path_data,"env/data/5_strain_unikmer",env,"unikmer")
# if not os.path.exists(path_jf_unikmer):
#     os.makedirs(path_jf_unikmer)

# #将得到的strain的unikmer排序并且去掉空值
path_jf_unikmer_sorted = os.path.join(path_data,"env/data/5_strain_unikmer",env,"unikmer_sorted_nonull")
# if not os.path.exists(path_jf_unikmer_sorted):
#     os.makedirs(path_jf_unikmer_sorted)
    
path_jf_kmer_sp = os.path.join(path_data,"env/data/5_strain_unikmer",env,"kmer")
# if not os.path.exists(path_jf_kmer_sp):
#     os.makedirs(path_jf_kmer_sp)
#虽然不需要把strain的kmer单独保存了，但是这个文件夹还是有用的，这个是用来存储kmer的df的
    
def get_negative_strand(sequence): #得到反向互补链
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reversed_sequence = sequence[::-1]  # 反向排列序列
    negative_strand = ''.join([complement[base] for base in reversed_sequence])
    return negative_strand

def get_ne_set(list_unikmer): #得到反向互补列表
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    translation_table = str.maketrans(complement)   
    list_ne_unikmer =[ i[::-1].translate(translation_table)  for i in list_unikmer ]
    return set(list_ne_unikmer)
path_jf_batch = os.path.join(path_jf,batch_name) #存放当前物种的sp,每次遍历新物种，都会删除重建


for ls in list_sp[start:end]: # /data/huixingqi/data/env/data/4_sp_strain_path_sample/
    path_ls = os.path.join(path_kmer_root,ls)
    print(path_ls)
    # path_sp_kmer = os.path.join(path_jf_kmer_sp,ls[:-4])
    # if not os.path.exists(path_sp_kmer):
    #     os.makedirs(path_sp_kmer)
    # cmd_jf = "bash 5_2_jellyfish.sh "+path_ls+" "+env+" "+str(k_kmer)+" "+batch_name+" "+path_sp_kmer
    cmd_jf = "bash 5_2_jellyfish.sh "+path_ls+" "+env+" "+str(k_kmer)+" "+batch_name+" "+ls.replace(".txt","")+" "+path_home
    subprocess.run(cmd_jf,shell=True,cwd=path_script) #得到中间结果了
    print(cmd_jf)
    dict_kmer = dict() # 存放sp下每个strain的kmer
    dict_strain = dict() # 存放每个sp的dataframe
    
    # if not os.exists(path_jf_new):
    #     os.makedirs(path_jf_new)
    list_jf = [f for f in os.listdir(path_jf_batch) if f.endswith(".fa") ]
    #该物种所有strain的.fa文件，记录strain的kmer
    
    path_unikmer_sp = os.path.join(path_jf_unikmer,batch_name,ls.replace(".txt",""))
    # if os.path.exists(path_unikmer_sp):
    #     shutil.rmtree(path_unikmer_sp)
        
    # os.makedirs(path_unikmer_sp)
    # print(path_unikmer_sp)
    
    set_sp_kmer = set()
    for j in list_jf:
        path_jf_out = os.path.join(path_jf_batch,j) 
         # eg /data/huixingqi/software/1_tools/kraken2/kraken2/res_all/tmpdata5/kmer/sp/GCA_000577295.1.fa
        df_fa = pd.read_csv(path_jf_out,header=None,sep="\t")
        cnt_fa_name = "cnt"+j
        df_fa.columns=["unikmer",cnt_fa_name]
        strain_kmer = set(df_fa["unikmer"].to_list()) #strain的kmer,我们得到的新的.fa文件只有序列，没有>name
        set_sp_kmer = set_sp_kmer | strain_kmer
        # strain_kmer_ne = get_ne_set(strain_kmer) #已经是set了
        # dict_kmer[j]=strain_kmer | strain_kmer_ne #一个j一个记录
        dict_kmer[j]=strain_kmer #不需要获得反向互补链了，因为默认链就是最小的那个
        dict_strain[j]=df_fa
        
    count_kmer = Counter([acns_kmer for acns in dict_kmer.values() for acns_kmer in acns])
    #记录每个kmer在所有strian中出现的次数，出现次数为1的是unikmeri
    flag_merge = False              
    for key,val in dict_kmer.items():
        unique_kmers = [ element for element in val if count_kmer[element] == 1]
        list_path = os.path.join(path_unikmer_sp,key)
        list_strain_df = pd.DataFrame(unique_kmers)
        if len(unique_kmers) > 0:
            # if flag_merge == False:
            #     df_sp_kmer_merge = dict_strain[key].copy()
            #     flag_merge = True
            # else:
            #     df_sp_kmer_merge = pd.merge(df_sp_kmer_merge,dict_strain[key],on="unikmer", how='outer')
            list_strain_df.columns=["unikmer"]
            list_strain_df_cnt = pd.merge(dict_strain[key],list_strain_df,on=["unikmer"])
            list_strain_df_cnt.to_csv(list_path,header=None,index=None,sep="\t")
    # columns_to_mean = df_sp_kmer_merge.columns.drop('unikmer')  
    # df_sp_kmer_merge['cnt_mean'] = df_sp_kmer_merge[columns_to_mean].mean(axis=1, skipna=True)
    # df_sp_kmer_merge_save = df_sp_kmer_merge.loc[:,["unikmer","cnt_mean"]]
    # path_sp_df = os.path.join(path_jf_kmer_sp,ls[:-4]+".csv")
    
    # df_sp_kmer_merge_save.to_csv(path_sp_df,header=None,sep="\t",index=None)
    cmd_sort = "bash 5_2_strain_kmer_sort.sh "+ls.replace(".txt","")+" "+env+" "+batch_name+" "+path_home
    subprocess.run(cmd_sort, shell=True, cwd=path_script)
    print(cmd_sort)

#全部跑完之后再sort，这样sort的输入就比较简单，如果一个sp一次sort，那么还需要输入sp
#但是既然jellyfish需要输入物种文件，那么sort也要吧
#那最后就定下来一个物种一次sort
# cmd_sort = "bash /data/huixingqi/sim_meta/all/5_2_strain_kmer_sort.sh "+path_ls+" "+env
# subprocess.run(cmd_jf,shell=True,cwd="/data/huixingqi/sim_meta/all/")


