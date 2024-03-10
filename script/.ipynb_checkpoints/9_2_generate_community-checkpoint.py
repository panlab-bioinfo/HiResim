#现在已经生成了全部的数据
#根据前面生成的数据生成模拟软件
#根据save_meanshift的结果生成模拟的数据
#  python 9_2_generate_community.py --env test --n 2
import pickle
import sys
import numpy as np
import pandas as pd
import os
from lib.generate_genus_and_sp import * 
import argparse
parser = argparse.ArgumentParser(description = "9_2 command line")
parser.add_argument("--env",help = "enviroment")
parser.add_argument("--n",type = int,default = 1,help = "n_samples")
args = parser.parse_args()
env = args.env
n_sample = args.n

abs_path = os.path.abspath(sys.argv[0])
abs_path_dir_script = os.path.dirname(abs_path) #HiResim/script
abs_path_dir_data = abs_path_dir_script.replace("script","data") #HiResim/data

path_sp_cnt = os.path.join(abs_path_dir_data,"env/data/9_sp_cnt/",env)
path_abu_root = os.path.join(path_sp_cnt,"genus_train")
path_sp_g = os.path.join(abs_path_dir_data,"dict_sp_g")

path_meanshift = os.path.join(path_sp_cnt,"meanshift")
# print(path_meanshift)

cnt_g,sp_cnt_list,lg = get_sp_cnt_new(path_meanshift)
# print(g_cnt)

#采样属
#属名称
path_key = os.path.join(path_sp_cnt,"genus","g_name_list.pkl")
with open(path_key,"rb") as f:
    key_list = pickle.load(f)
#属频率
path_prob = os.path.join(path_sp_cnt,"genus","g_abu_list.pkl")
with open(path_prob,"rb") as f:
    prob_list = pickle.load(f)
    
#属下物种名称   
path_g_sp = os.path.join(path_sp_cnt,"genus","g_sp.pkl")    
with open(path_g_sp,"rb") as f:
    dict_new_sp = pickle.load(f)
    
#属下物种丰度  
path_g_sp_abu = os.path.join(path_sp_cnt,"genus","g_sp_abu.pkl")   
with open(path_g_sp_abu,"rb") as f:
    dict_new_abu= pickle.load(f)
    
sample = choice_genus(key_list,prob_list,cnt_g)#得到的属的列表
# print(sample)
# print(len(sample),cnt_g)

#采样物种
path_sim = os.path.join(abs_path_dir_data,"env/data/10_sim/",env)
path_sp = os.path.join(path_sim,"sp_generate")
path_sp_strain_cnt_root = os.path.join(path_sim,"sp_strain_cnt")
path_sp_strain_abu_root = os.path.join(path_sim,"sp_strain_abu")

path_sp_strain_abu_sp_root = os.path.join(path_sim,"sp_strain_abu_sp")
if not os.path.exists(path_sp):
    os.makedirs(path_sp)
if not os.path.exists(path_sp_strain_cnt_root):
    os.makedirs(path_sp_strain_cnt_root)
if not os.path.exists(path_sp_strain_abu_root):
    os.makedirs(path_sp_strain_abu_root)
    
if not os.path.exists(path_sp_strain_abu_sp_root):
    os.makedirs(path_sp_strain_abu_sp_root)
path_model_root =  os.path.join(path_sp_cnt,"ctgan")

#生成strain丰度和个数

#个数
path_strain = os.path.join(abs_path_dir_data,"env/data/8_strain_cnt_and_abu/",env)
path_strain_cnt = os.path.join(path_strain,"all","sp_cnt_all.csv")
# print(path_strain_cnt)
df_strain_cnt = pd.read_csv(path_strain_cnt,header=None,sep="\t")
df_strain_cnt.columns=["sp","cnt"]
#丰度
path_strain_abu = os.path.join(path_strain,"all","sp_abu_all.csv")
# print(path_strain_abu)
df_strain_abu = pd.read_csv(path_strain_abu,header=None,sep="\t")
df_strain_abu.columns=["sp","abu"]


#获得ani
path_ani = os.path.join(abs_path_dir_data,"env/data/7_ani/",env)
path_matrix_root = os.path.join(path_ani,"matrix")
path_tree_root = os.path.join(path_ani,"tree")

path_sp_strain_ani_root = os.path.join(path_sim,"sp_strain_ani")
if not os.path.exists(path_sp_strain_ani_root):
    os.makedirs(path_sp_strain_ani_root)

path_sp_strain_tree_root = os.path.join(path_sim,"sp_strain_tree")
if not os.path.exists(path_sp_strain_tree_root):
    os.makedirs(path_sp_strain_tree_root)
#这两个的物种数不一样多，df_sp_cnt_mean没有根据strain删减
df_sp_cnt_mean_nodrop = get_sp_cnt_mean(df_strain_cnt)
df_sp_sample_ani = get_sp_ani_path(path_tree_root) 
df_sp_cnt_mean = df_sp_cnt_mean_nodrop.loc[df_sp_cnt_mean_nodrop["sp"].isin(df_sp_sample_ani["sp"].unique())]
# print(len(df_sp_cnt_mean["sp"].unique()))
# print(len(df_sp_sample_ani["sp"].unique()))
#多样本生成

for n in range(n_sample):

    list_sp = choice_sp(path_model_root,sample,dict_new_sp) #得到的物种的列表
    
    #生成物种的丰度
    dict_sp_abu = get_sp_abu(list_sp,path_abu_root,path_sp_g)
    df_sp_abu = pd.DataFrame(list(dict_sp_abu.items()),columns=["sp","abu"])
    path_sp_strain_abu_sp = os.path.join(path_sp_strain_abu_sp_root,"sp_strain_abu_"+str(n+1)+".csv")
    df_sp_abu.to_csv(path_sp_strain_abu_sp,header=None,sep="\t",index=None)
    
    #保存物种结果
    path_sp_generate = os.path.join(path_sp,"sp_generate_"+str(n+1)+".pkl")
    with open(path_sp_generate,"wb") as f:
        pickle.dump(list_sp,f)
    #cnt和abu
    list_sp_cnt,list_sp_abu = strain_cnt_abu(list_sp,df_strain_cnt,df_strain_abu)
    # print(list_sp[:2])
    # print(list_sp_cnt[:2])
    print(len(list_sp),len(list_sp_cnt),len(list_sp_abu))
    
    path_sp_strain_cnt = os.path.join(path_sp_strain_cnt_root,"sp_strain_cnt_"+str(n+1)+".pkl") 
    with open(path_sp_strain_cnt,"wb") as f:
        pickle.dump(list_sp_cnt,f)
        
    path_sp_strain_abu = os.path.join(path_sp_strain_abu_root,"sp_strain_abu_"+str(n+1)+".pkl")
    with open(path_sp_strain_abu,"wb") as f:
        pickle.dump(list_sp_abu,f)
    #ani
    path_sp_strain_ani = os.path.join(path_sp_strain_ani_root,"sim_ani_"+str(n+1))
    if not os.path.exists(path_sp_strain_ani):
        os.makedirs(path_sp_strain_ani)
        
    path_sp_strain_tree = os.path.join(path_sp_strain_tree_root,"sim_tree_"+str(n+1))
    if not os.path.exists(path_sp_strain_tree):
        os.makedirs(path_sp_strain_tree)
        
    mask = df_sp_cnt_mean["sp"].isin(list_sp)
    df_sp_cnt_mean_new = df_sp_cnt_mean.drop(df_sp_cnt_mean[mask].index)
    #这是不在模拟生成的那些物种，用这些物种来拟合模拟的sp却没有ani的情况
#     if len(list_sp)+len(df_sp_cnt_mean_new)!=len(df_sp_cnt_mean):
#         #因为list_sp是未删减的sp，而df_sp_cnt_mean中那些不符合条件的sp被去掉了
#         #所以二者之和可能不相等
        #但是删掉了，不影响
        
#         print("len",len(list_sp),len(df_sp_cnt_mean_new),len(df_sp_cnt_mean))
#         print("mask error")
    for sp in list_sp:
        index_sp = list_sp.index(sp)
        cnt_all = list_sp_cnt[index_sp]
        tree_path_list = df_sp_sample_ani[df_sp_sample_ani["sp"]==sp]["path"].to_list()
        if len(tree_path_list)==0:
            sp_alt = find_alt_sp(cnt_all,df_sp_cnt_mean_new)
            # print(sp_alt)
            # print(df_sp_cnt_mean[df_sp_cnt_mean["sp"]==sp_alt])
            tree_path_list = df_sp_sample_ani[df_sp_sample_ani["sp"]==sp_alt]["path"].to_list()
            # print(df_sp_sample_ani[df_sp_sample_ani["sp"]==sp_alt]["path"].to_list())
        # print(len(tree_path_list))
        generate_strain_ani_new(path_sp_strain_ani,path_sp_strain_tree,sp,tree_path_list,cnt_all)
        
    
    


    
    