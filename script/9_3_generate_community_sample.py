#根据某个环境以及样本，生成该样本中所有物种的strain
#输入是环境和物种列表的地址

#  python 9_2_generate_sample.py --env test 
import pickle
import sys
import numpy as np
import pandas as pd
import os
from lib.generate_genus_and_sp import * 
import argparse
parser = argparse.ArgumentParser(description = "9_2 command line")
parser.add_argument("--env",help = "enviroment")
parser.add_argument("--sample",help = "GTDB species list")
args = parser.parse_args()
env = args.env
sp_path = args.sample

abs_path = os.path.abspath(sys.argv[0])
abs_path_dir_script = os.path.dirname(abs_path) #HiResim/script
abs_path_dir_data = abs_path_dir_script.replace("script","data") #HiResim/data

sample_name = sp_path.split("/")[-1].split(".")[0]
path_sp_cnt = os.path.join(abs_path_dir_data,"env/data/9_sp_cnt/",env)
path_abu_root = os.path.join(path_sp_cnt,"genus_train")
path_sp_g = os.path.join(abs_path_dir_data,"dict_sp_g")
#这个后面要改成非绝对值地址
path_sim = os.path.join(abs_path_dir_data,"env/data/10_sim/",env+"_sample")
path_sp_strain_abu_sp_root = os.path.join(path_sim,"sp_strain_abu_sp")
if not os.path.exists(path_sp_strain_abu_sp_root):
    os.makedirs(path_sp_strain_abu_sp_root)

df_sample_sp = pd.read_csv(sp_path,header=None,sep="\t")
df_sample_sp.columns = ["sp"]
list_sp = df_sample_sp["sp"].to_list()
#生成物种的丰度
dict_sp_abu = get_sp_abu(list_sp,path_abu_root,path_sp_g)
df_sp_abu = pd.DataFrame(list(dict_sp_abu.items()),columns=["sp","abu"])
path_sp_strain_abu_sp = os.path.join(path_sp_strain_abu_sp_root,"sp_strain_abu.csv")
df_sp_abu.to_csv(path_sp_strain_abu_sp,header=None,sep="\t",index=None)
# print(df_sp_cnt_mean_nodrop[:2])

#采样物种
path_sp = os.path.join(path_sim,"sp_generate_sample")
path_sp_strain_cnt_root = os.path.join(path_sim,"sp_strain_cnt")
path_sp_strain_abu_root = os.path.join(path_sim,"sp_strain_abu")
if not os.path.exists(path_sp):
    os.makedirs(path_sp)
if not os.path.exists(path_sp_strain_cnt_root):
    os.makedirs(path_sp_strain_cnt_root)
if not os.path.exists(path_sp_strain_abu_root):
    os.makedirs(path_sp_strain_abu_root)

#strain个数
path_strain = os.path.join(abs_path_dir_data,"env/data/8_strain_cnt_and_abu/",env)
path_strain_cnt = os.path.join(path_strain,"all","sp_cnt_train.csv")
df_strain_cnt = pd.read_csv(path_strain_cnt,header=None,sep="\t")
df_strain_cnt.columns=["sp","cnt"]

df_sp_cnt_mean_nodrop = get_sp_cnt_mean(df_strain_cnt)#样本中这个sp出现了多次，我们将样本中的随机采样
# print(len(df_sp_cnt_mean_nodrop))

#strain丰度
path_strain_abu = os.path.join(path_strain,"all","sp_abu_train.csv")
df_strain_abu = pd.read_csv(path_strain_abu,header=None,sep="\t")
df_strain_abu.columns=["sp","abu"]

#获得ani
path_ani = os.path.join(abs_path_dir_data,"env/data/7_ani/",env)
path_matrix_root = os.path.join(path_ani,"matrix")
path_tree_root = os.path.join(path_ani,"tree")

#输出结果
path_sp_strain_ani_root = os.path.join(path_sim,"sp_strain_ani")
if not os.path.exists(path_sp_strain_ani_root):
    os.makedirs(path_sp_strain_ani_root)
    
path_sp_strain_tree_root = os.path.join(path_sim,"sp_strain_tree")
if not os.path.exists(path_sp_strain_tree_root):
    os.makedirs(path_sp_strain_tree_root)

df_sp_sample_ani = get_sp_ani_path(path_tree_root) #这个是返回这个物种在哪些样本中出现了
df_sp_cnt_mean = df_sp_cnt_mean_nodrop.loc[df_sp_cnt_mean_nodrop["sp"].isin(df_sp_sample_ani["sp"].unique())]
# print(len(df_sp_cnt_mean))

#生成样本

#保存物种结果
path_sp_generate = os.path.join(path_sp,"sp_generate"+"_sample.pkl")
path_sp_generate_csv = os.path.join(path_sp,"sp_generate"+"_sample.csv")
print(path_sp_generate,path_sp_generate_csv)
df_list_sp = pd.DataFrame(list_sp)
print("len",len(list_sp),len(df_list_sp))
df_list_sp.to_csv(path_sp_generate_csv,header=None,index=False,sep="\t")

with open(path_sp_generate,"wb") as f:
    pickle.dump(list_sp,f)
    
#cnt and abu
list_sp_cnt,list_sp_abu = strain_cnt_abu(list_sp,df_strain_cnt,df_strain_abu)

#cnt
path_sp_strain_cnt = os.path.join(path_sp_strain_cnt_root,"sp_strain_cnt"+"_sample.pkl") 
path_sp_strain_cnt_csv = os.path.join(path_sp_strain_cnt_root,"sp_strain_cnt"+"_sample.csv") 
df_list_sp_cnt = pd.DataFrame(list_sp_cnt)
print("len",len(list_sp_cnt),len(df_list_sp_cnt))
df_list_sp_cnt.to_csv(path_sp_strain_cnt_csv,header=None,index=False,sep="\t")
with open(path_sp_strain_cnt,"wb") as f:
    pickle.dump(list_sp_cnt,f)

#abu
path_sp_strain_abu = os.path.join(path_sp_strain_abu_root,"sp_strain_abu"+"_sample.pkl")
path_sp_strain_abu_csv = os.path.join(path_sp_strain_abu_root,"sp_strain_abu"+"_sample.csv")
df_list_sp_abu = pd.DataFrame(list_sp_abu)
print("len",len(list_sp_abu),len(df_list_sp_abu))
df_list_sp_abu.to_csv(path_sp_strain_abu_csv,header=None,index=False,sep="\t")
with open(path_sp_strain_abu,"wb") as f:
    pickle.dump(list_sp_abu,f)
    

#ani
path_sp_strain_ani = os.path.join(path_sp_strain_ani_root,"sim_ani")
if not os.path.exists(path_sp_strain_ani):
    os.makedirs(path_sp_strain_ani)

path_sp_strain_tree = os.path.join(path_sp_strain_tree_root,"sim_tree")
if not os.path.exists(path_sp_strain_tree):
    os.makedirs(path_sp_strain_tree)
        
mask = df_sp_cnt_mean["sp"].isin(list_sp)
df_sp_cnt_mean_new = df_sp_cnt_mean[mask] #之前是不在，现在是在
for sp in list_sp:
    index_sp = list_sp.index(sp)
    cnt_all = list_sp_cnt[index_sp]
    if cnt_all>1: #当cnt_all==1时，不需要生成strain,但是明天要看看是不是真的把这个也考虑到了
        
        tree_path_list = df_sp_sample_ani[df_sp_sample_ani["sp"]==sp]["path"].to_list()
        if len(tree_path_list)==0:
            sp_alt = find_alt_sp(cnt_all-1,df_sp_cnt_mean_new.copy())
            tree_path_list = df_sp_sample_ani[df_sp_sample_ani["sp"]==sp_alt]["path"].to_list()
        # if cnt_all == 1:
            # print(sp,tree_path_list)
        generate_strain_ani_new(path_sp_strain_ani,path_sp_strain_tree,sp,tree_path_list,cnt_all-1)#遍历最小生成树，生成物种间的ani
