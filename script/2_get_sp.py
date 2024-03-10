#对于得到的out和rep进行处理
#根据物种的丰度队伍中进行截断
# python 2_get_sp.py /data/huixingqi/data/env/data/abundance/test/ /data/huixingqi/software/1_tools/kraken2/kraken2/res_all/tmpdata5/file/read_cnt.txt
# python 2_get_sp.py marine 0 10
import numpy as np
import pandas as pd
import os
import sys

env = sys.argv[1]
# start = int(sys.argv[2])
# end = int(sys.argv[3])
abu_thr = float(sys.argv[2]) #一般是0.001
path_home = sys.argv[3] # HiResim/
path_data = os.path.join(path_home,"data")
path_data_env = os.path.join(path_data,"env/data") #HiResim/data/env/data
# path_abundance_root=sys.argv[1]
# path_read_cnt = sys.argv[2]
print("start")
path_abundance_root = os.path.join(path_data_env,"1_abundance/",env)
print("path_abundance_root",path_abundance_root)
path_read_cnt =  os.path.join(path_data_env,"1_read_cnt",env+"_read_cnt.txt")
path_cut =  os.path.join(path_data_env,"2_sp_cut",env)


if os.path.exists(path_cut)==False:
    os.makedirs(path_cut)
    
def get_sp_list(path_read_cnt,path_abundance_root,path_cut,path_data = path_data):
    print("get_sp_list")
    df_read_cnt = pd.read_csv(path_read_cnt,sep=" ",header=None)
    df_read_cnt.columns=["read","sample"]
    # read_cnt_list = df_read_cnt.loc[:,"read"].to_list()#[start:end]
    df_read_cnt["name"]=df_read_cnt["sample"].apply(lambda x:x.split("/")[-1].replace(".out",""))
    dict_read_cnt = df_read_cnt.set_index("name")["read"].to_dict()
    
    
    
    path_s1 = os.path.join(path_data,"file/s1_ancs_database.tsv")
    s1_sp = pd.read_csv(path_s1,header=None,sep="\t",dtype=str)
    s1_sp.columns = ["name","s1"] 
    #得到sp和s1的对应关系，然后和每个样本进行merge,把merge的结果保存在sp_cut文件夹中
    abundance_name = sorted(os.listdir(path_abundance_root))#[start:end]
    #得到所有的rep.txt
    # cnt = 0
    for i in abundance_name:
        key = i.replace("_abundance.txt","")
        name = i.replace("_abundance.txt","_cutsp.txt")
        path_name = os.path.join(path_cut,name)
        path_abu = os.path.join(path_abundance_root,i)
        # read_cnt_flag = int(read_cnt_list[cnt]*abu_thr)
        rep = pd.read_csv(path_abu,sep="\t",header=None)
        rep.columns =["name","abu"]
        list_sp_cnt = rep["abu"].to_list()
        read_cnt_flag = int(sum(list_sp_cnt)*abu_thr) 
        # cnt+=1 
        # if read_cnt_flag>=2000:
        if read_cnt_flag>0:

            mask = (rep["abu"]<=read_cnt_flag) 
            rep_drop = rep.drop(rep[mask].index)
            df_merge  = pd.merge(rep_drop,s1_sp,on=["name"])
            df_merge_c2 = df_merge.loc[:,["name","s1"]]
            df_merge_c2.to_csv(path_name,sep="\t",index=None,header=None)
            #下面是之前的方法，根据个数判断去留
            # len_merge = len(df_merge_c2)
            # if len_merge>=10 and len_merge<300:
            # print(path_name)
                # df_merge_c2.to_csv(path_name,sep="\t",index=None,header=None)
            
        

get_sp_list(path_read_cnt,path_abundance_root,path_cut) 