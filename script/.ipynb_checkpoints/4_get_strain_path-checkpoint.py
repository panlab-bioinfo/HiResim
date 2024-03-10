#获得strain基因组的kmer的地址
#地址有两种，一种是strain基因组原始地址，用于跑jellfisyh
#另一种是存储地址，用于保存跑完的结果
#该脚本只生产地址，具体跑kmer看5以及5_2
import pandas as pd
import os
import sys
import shutil
#这个没必要分批次了，一方面是很快，另一方面是因为生成sp_.txt时需要考虑全部样本
#存放物种sp所有的基因组

env = sys.argv[1]
path_home = sys.argv[2] # HiResim/
path_data = os.path.join(path_home,"data")
path_tools = os.path.join(path_home,"tools")

# start = int(sys.argv[2])
# end = int(sys.argv[3])
#这个没必要分批次了，一方面是很快，另一方面是因为生成sp_.txt时需要考虑全部样本
#存放物种sp所有的基因组

strain_path = os.path.join(path_tools,"strain_path.csv")
df_strain = pd.read_csv(strain_path,sep="\t",header=None) 
df_strain.columns=["sp","acns","down","genome"]
df_strain["path"] = df_strain["genome"].apply(lambda x:os.path.join(path_tools,"strain_genome",x))

path_sp_cut = os.path.join(path_data,"env/data/2_sp_cut/",env[:-2])
#这个cut后的文件其实删减过不符合条件的样本了

#保存样本的物种strain的地址,这是所有的strian，没有经过删减过的
sample_strain_root = os.path.join(path_data,"env/data/4_sp_strain_path_all",env)
if not os.path.exists(sample_strain_root):
    os.makedirs(sample_strain_root)

#保存每个strain的kmer，按照物种存，即所有样本一起来
#因为这个只是生成了路径，没有得到unikmer。是5_得到的kemr,所有文件5开头
kemr_root = os.path.join(path_data,"env/data/5_strain_unikmer/",env,"unikmer")
if not os.path.exists(kemr_root):
    os.makedirs(kemr_root)
        
list_dir = sorted(os.listdir(path_sp_cut))
# list_dir_sub = list_dir[start:end]

#用来存放每个样本的strain,用于取并集
# list_strains = [] 

set_all = set()
for sp in list_dir:
    path_sp = os.path.join(path_sp_cut,sp)
    df_sp = pd.read_csv(path_sp,header=None,sep="\t") # cut后物种 sp s1
    df_sp.columns=["sp","s1"]
    df_merge = pd.merge(df_sp,df_strain,on=["sp"]) #df_merge之后的结果就已经把不符合条件的sp去掉了（有些sp没有strain）
    df_merge_sub = df_merge.loc[:,["sp","path"]] #merge的结果完全没问题，不需要分组再遍历了
    p_set = set(df_merge_sub["path"])
    set_all |= p_set
    sample_strain_path = os.path.join(sample_strain_root,sp.replace("_cutsp.txt",".txt"))
    
    df_merge_sub["path"] = df_merge_sub["path"].apply(lambda x:os.path.join(kemr_root,x.split("/")[-1][:15]+".fa"))
    df_merge_sub.to_csv(sample_strain_path,header=None,index=None,sep="\t")
    
set_all_df = pd.DataFrame(list(set_all))
set_all_df.columns = ["path"]
set_all_df_merge = pd.merge(df_strain,set_all_df,on=["path"])
set_all_df_merge_gp = set_all_df_merge.groupby("sp")
path_kmer_merge_sp_root = os.path.join(path_data,"env/data/4_sp_strain_path_sample",env)
#这个里面存的是所有样本的并集，所以如果之前有了，就不用生成了
if not os.path.exists(path_kmer_merge_sp_root):
    os.makedirs(path_kmer_merge_sp_root)

for group_name, group_indices in set_all_df_merge_gp.groups.items():
    group = set_all_df_merge.loc[group_indices]
    group_path = group["path"]
    path_sp_path = os.path.join(path_kmer_merge_sp_root,group_name+".txt")
    # if not os.path.exists(path_sp_path):
    group_path.to_csv(path_sp_path,header=None,index=False)

