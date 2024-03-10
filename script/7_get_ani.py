#得到了每个样本物种strain的基因组
#需要计算基因组之间的ani
#首先要确定哪些基因组需要计算ani
#我们当初得到unikemr的时候是把样本间所有的物种合在一起看的，而每个样本的物种在sp_cut中保存
#主要思路是遍历每一个样本，遍历该样本的物种，物种内stain计算ani,
#结果存入字典中，字典的key为物种名，值为dataframe,第一列为ref基因组，第二列为query基因组

# python 7_get_ani.py oral

import pandas as pd
import numpy as np
import os
import itertools
import sys
import subprocess

env = sys.argv[1]
path_home = sys.argv[2]
path_data = os.path.join(path_home,"data")
path_script = os.path.join(path_home,"script")
path_tools = os.path.join(path_home,"tools")
# start = int(sys.argv[2])
# end = int(sys.argv[3])
#这个也不能batch因为要遍历所有样本求出要query和ref的组合数，然后去掉样本间重复的

#这个是65703+14万个
# path_genome = "/data/huixingqi/sim_meta/script_sum/data/sp_path_drop_fromOld.txt"
path_genome = os.path.join(path_tools,"strain_path.csv")
path_inter_root = os.path.join(path_data,"env/data/6_read_strain_inter",env+"_new")
path_ani = os.path.join(path_data,"env/data/7_ani",env)
if not os.path.exists(path_ani):
    os.makedirs(path_ani)
    
#得到所有要比对的strain基因组的地址，用于sketch
df_genome = pd.read_csv(path_genome,header=None,sep="\t")
# df_genome.columns=["sp","path"]
df_genome.columns=["sp","acns","down","genome"]
df_genome["path"] = df_genome["genome"].apply(lambda x:os.path.join(path_tools,"strain_genome",x))
dict_acns_path = dict()
for path_g in df_genome["path"]:
    acns = path_g.split("/")[-1][:15]
    dict_acns_path[acns] = path_g 
    #获得每个strain基因组的路径
    
#看看有多少strain需要mash  
list_acns = [] #保存acns，没去重版本的
list_path_acns = [] #保存去重后acns地址
set_all =set()

#得到要比较的ref和query
#存放df_ref 和query,df有两列，一列是ref,一列是query
df_rq = pd.DataFrame() 
list_srr = sorted(os.listdir(path_inter_root))
for srr in list_srr: #SRR5936111
    path_srr_inter = os.path.join(path_inter_root,srr)
    
    for sp in sorted(os.listdir(path_srr_inter)):
        path_sp_inter = os.path.join(path_srr_inter,sp)
        list_strain = sorted(os.listdir(path_sp_inter))
        #地址
        list_acns += list_strain
        #ref和query
        ref_query = list([ list(rq) for rq in itertools.combinations(list_strain,2)]) #得到c(n,2)
        if len(ref_query)!=0:
            df_sp_rq = pd.DataFrame(ref_query,columns=["ref","query"])
            df_rq = pd.concat([df_rq,df_sp_rq])
    df_rq = df_rq.drop_duplicates()
        # for strain in sorted(os.listdir(path_sp_inter)):
        #     path_strain = os.path.join(path_sp_inter,strain)
        #     if os.stat(path_strain).size!=0:
        #         list_acns_os.append(strain)  
set_acns = set(list_acns)
list_path_acns = [dict_acns_path[acn.replace(".fa","")] for acn in set_acns]
# for acn in set_acns:
#     list_path_acns.append(dict_acns_path[acn.replace(".fa","")])
# print(len(set_acns),len(list_acns),len(list_path_acns))

df_list_acns = pd.DataFrame(list_path_acns)
path_df_strain = os.path.join(path_ani,"strain.txt")

df_list_acns.to_csv(path_df_strain,header=None,index=None)

df_rq["ref_path"] = df_rq["ref"].apply(lambda x:dict_acns_path[x.replace(".fa","")])
df_rq["query_path"] = df_rq["query"].apply(lambda x:dict_acns_path[x.replace(".fa","")])
path_df_rq = os.path.join(path_ani,"rq.csv")
df_rq.to_csv(path_df_rq,index=None,sep="\t",header=None)
#得到df_rq后好像可以batch了
#但是rq.csv和strain.txt行数不一样，不太好batch

#得到要执行的strian基因组后，先跑sketch，后跑mask

path_sketch= os.path.join(path_ani,"sketch")
if not os.path.exists(path_sketch):
    os.makedirs(path_sketch)
cmd_sketch = "bash 7_2_sketch_and_mash.sh "+env +" "+path_home
subprocess.run(cmd_sketch,shell=True,cwd=path_script) 

#7_2_sketch_and_mash.sh得到了两两ani结果存在ani.csv中

#每个样本，每个物种构建距离矩阵
# path_df_drop_root = "/data/huixingqi/software/1_tools/kraken2/kraken2/res_all/tmpdata5/sp_unikmer_drop/"
path_df_mash = os.path.join(path_ani,"ani.csv") #mash的结果
df_mash = pd.read_csv(path_df_mash,sep="\t",header=None)
df_mash.columns=["ref","query","dist","p","rate"]
#因为df_mash是前面执行完才生成的
#所以不能把下面的流程嵌合到前面
path_matrix = os.path.join(path_ani,"matrix")#保存dis的路径

#首先df_mash里面保存的是基因组地址，是acns根据dict_acns_path[acns]转换得到的
#我们就记录名字了，只需要得到矩阵

for srr in sorted(os.listdir(path_inter_root)): #SRR5936111
    path_srr_inter = os.path.join(path_inter_root,srr)
    
    path_srr_matrix = os.path.join(path_matrix,srr)
    if not os.path.exists(path_srr_matrix):
        os.makedirs(path_srr_matrix)
    
    for sp in sorted(os.listdir(path_srr_inter)):
        path_sp_inter = os.path.join(path_srr_inter,sp)
        list_strain = sorted(os.listdir(path_sp_inter))
        path_sp_matrix = os.path.join(path_srr_matrix,sp+".txt")
        #地址
        #ref和query
        ref_query = list([ list(rq) for rq in itertools.combinations(list_strain,2)]) #得到c(n,2)
        if len(ref_query)!=0:
            df_sp_rq = pd.DataFrame(ref_query,columns=["ref","query"])
            # df_rq = pd.concat([df_rq,df_sp_rq])
            df_sp_rq["ref"] = df_sp_rq["ref"].apply(lambda x:dict_acns_path[x.replace(".fa","")])
            df_sp_rq["query"] = df_sp_rq["query"].apply(lambda x:dict_acns_path[x.replace(".fa","")])
            df_merge = pd.merge(df_sp_rq,df_mash,on=["ref","query"])
            #不能只merge ref 因为可能你有 （s1 s2) 别人有[(s1 s2) (s1 s3)] 这样就多了
            distances = df_merge["dist"].tolist()
            n_acns = len(list_strain)
            distance_matrix = np.zeros((n_acns, n_acns))
            index = 0
            for i in range(n_acns-1):
                for j in range(i+1, n_acns):
                    distance_matrix[i][j] = distances[index]
                    distance_matrix[j][i] = distances[index]
                    index += 1
            np.savetxt(path_sp_matrix, distance_matrix)
        
            
