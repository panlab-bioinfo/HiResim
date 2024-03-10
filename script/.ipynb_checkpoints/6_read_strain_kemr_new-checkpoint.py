#因为strain的kmer也加上了个数，对得到的交集处理一下

import pandas as pd
import os
import pickle
import sys
# nohup python 6_read_strain_kemr_new.py rhizosphere14 &> log/run_6_read_strain_kemr_new_rhizosphere14.log &
env = sys.argv[1]
path_home = sys.argv[2]
path_data = os.path.join(path_home,"data")

path_env = os.path.join(path_data,"env/data/6_read_strain_inter",env)
path_env_new = os.path.join(path_data,"env/data/6_read_strain_inter",env+"_new")

if not os.path.exists(path_env_new):
    os.makedirs(path_env_new)


list_srr = sorted(os.listdir(path_env))
for srr in list_srr:
    path_srr = os.path.join(path_env,srr)
    
    path_srr_new = os.path.join(path_env_new,srr)
    if not os.path.exists(path_srr_new):
        os.makedirs(path_srr_new)
    
    list_sp = sorted(os.listdir(path_srr))
    for sp in list_sp:
        
        path_sp = os.path.join(path_srr,sp)
        
        path_sp_new = os.path.join(path_srr_new,sp)
        if not os.path.exists(path_sp_new):
            os.makedirs(path_sp_new)
        
        list_fa = sorted(os.listdir(path_sp))
        for fa in list_fa:
            path_fa = os.path.join(path_sp,fa)
            
            path_fa_new = os.path.join(path_sp_new,fa)
            
            df_fa = pd.read_csv(path_fa,header=None,sep=" ")
            df_fa.columns=["kmer","cnt_read","cnt_strain"]
            df_fa['cnt'] = df_fa['cnt_read'] / df_fa['cnt_strain']
            df_fa_new = df_fa.loc[:,["kmer","cnt"]]
            df_fa_new.to_csv(path_fa_new,header=None,sep="\t",index=False)
            
            
        