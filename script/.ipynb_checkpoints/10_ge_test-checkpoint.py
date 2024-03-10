
import pandas as pd
import numpy as np
import os
import subprocess
from lib.generate_genus_and_sp import *
import sys
from Bio import SeqIO
from Bio.Seq import Seq

path_sim = sys.argv[1] #模拟群落的文件夹
# /data/huixingqi/data/env/data/10_sim/gut21_sample

# path_abs = sys.argv[2] #data文件的父目录
# /data/huixingqi/sim_meta/mhrec

out_name = sys.argv[2] #输出文件夹

path_abs_file = os.path.abspath(sys.argv[2])
path_abs = os.path.dirname(path_abs_file)
# sim_out
sim_out = os.path.join(path_abs,out_name) #创建输出文件夹的过程放到.sh中
print(sim_out)
if not os.path.exists(sim_out):
    os.makedirs(sim_out)
# /data/huixingqi/sim_meta/mhrec/sim_out
path_strain_sim = os.path.join(sim_out,"sim_strain") #这个用来保存生成的新的 strains
if not os.path.exists(path_strain_sim):
    os.makedirs(path_strain_sim)

depth= sys.argv[3] #pbsim用的depth 1


# #下载strain的地址文件夹
# path_down_strain = os.path.join(path_abs,"data/down_strain.csv")
# df_down = pd.read_csv(path_down_strain,header=None,sep="\t")
# df_down.columns = ["sp","acns","path","genome"]
#strain的绝对地址
path_down_strain = os.path.join(path_abs,"data/file/gtdb_ref.csv")
df_down = pd.read_csv(path_down_strain,header=None,sep="\t")
df_down.columns = ["sp","genome"]

#这个既有sp又有sp的丰度
path_sp_root = os.path.join(path_sim,"sp_strain_abu_sp")
#strain的个数
path_strain_cnt_root = os.path.join(path_sim,"sp_strain_cnt")
#strian的abu
path_strain_abu_root = os.path.join(path_sim,"sp_strain_abu")
#strain的ani
path_strain_ani_root = os.path.join(path_sim,"sp_strain_tree")

# print("path_sp_root",path_sp_root)
path_sp_list = sorted([ psr for psr in os.listdir(path_sp_root) if psr.endswith(".csv")])
print(path_sp_list)
path_strain_genome = os.path.join(path_abs,"data/strain_download") #下载的基因组的路径
for n_index,ps in enumerate(path_sp_list[:1]):
    # n_index+=9
    print(n_index,ps)
    path_sp = os.path.join(path_sp_root,ps)
    df_sp = pd.read_csv(path_sp,header=None,sep="\t")
    df_sp.columns = ["sp","abu"]
    path_strain_sim_n = os.path.join(path_strain_sim,"sim"+str(n_index))
    
    #abu列表
    list_abu = df_sp["abu"].to_list()
    sum_abu = sum(list_abu)
    list_abu_std = [lb/sum_abu for lb in list_abu]
    print(sum(list_abu_std),"\n")
    print(list_abu_std)
    
    
    path_abu = os.path.join(path_strain_abu_root,ps.replace(".csv",".pkl"))
    with open(path_abu,"rb") as f:
        abu_real = pickle.load(f) #长度为这个样本的物种个数
    # path_abu = os.path.join(path_strain_abu_root,ps)
    # df_abu = pd.read_csv(path_abu,header=None,sep="\t")
    #不能用dataframe，因为dataframe的缺点是 每一个物种一行，每行的列数不一样
      

    df_merge = pd.merge(df_sp,df_down,on=["sp"])

    #这个是全部strain的统计信息
    df_strain = df_merge[["sp","genome","abu"]] #这个是所有的物种
    df_strain = df_strain.reset_index(drop=True)
    
    df_strain["abu_new"] = df_strain["abu"].apply(lambda x:float(float(x)/sum_abu))
    # print(df_strain)
    #这个包含了模拟strian的丰度，明天需要根据sp下的strain个数和ani生成新的strain
    path_strain_df = os.path.join(sim_out,"sim_strain"+str(n_index)+".txt")
    df_strain_save = df_strain[["sp","genome","abu_new"]]
    df_strain_save.to_csv(path_strain_df,sep="\t",index=False)
    
    # list_genome_down = os.listdir(path_strain_genome) #已经下好的df
    # mask_down = df_merge["genome"].isin(list_genome_down)
    # df_merge_nodow = df_merge[~mask_down]
    
    # df_no_none = df_merge_nodow.dropna()
    # print(df_no_none[:2])
    # print("start to download ",len(df_no_none)," strain genomes")

    #保存这个模拟样本要下载的strain地址
#     path_tmp_strain = os.path.join(sim_out,"strain_download_"+str(n_index)+".txt")
#     # print(path_tmp_strain)
#     if len(df_no_none)!=0:
#         df_no_none.to_csv(path_tmp_strain,header=None,sep="\t",index=False)
    
#         path_down = os.path.join(path_abs,"script/lib/down_strains.py")
#         cmd_down = "python "+path_down+" "+ path_tmp_strain+" "+path_abs
#     # print(cmd_down)
#         try:
#             result = subprocess.run(cmd_down,shell=True,cwd="./")
#             result.check_returncode()
#         except subprocess.CalledProcessError as e:
#             print("download strain genomes error")
#             sys.exit(1)
    # print(df_sp)
    # print(len(df_sp))
    
    #下载完了参考的strain，现在需要根据参考的strain和ani生成新的
    
    #每个样本中物种ani的df
    ps_name = ps.replace(".csv","")
    path_ani_this = os.path.join(path_strain_ani_root,ps_name.replace("sp_strain_abu","sim_tree"))
    #这个sim_train有时候又是sim_sample,最后统一是sim_tree，别忘记修改生成的方法哦ing
    path_abu_pkl = os.path.join(path_strain_abu_root,ps.replace(".csv",".pkl"))
    with open(path_abu_pkl,"rb") as f:
        dict_sp_strain_abu = pickle.load(f)
    # print(dict_sp_strain_abu)
    list_sp = df_strain["sp"].to_list()
    for index_s,ls in enumerate(list_sp): 
        #这个是所有物种的 ani文件，但是有的物种没有ani文件，只有一个strain的物种就没有ani文件
        genome_sp = df_strain[df_strain["sp"]==ls]["genome"].unique()[0]
        
        # genome_path = os.path.join(path_strain_genome,genome_sp)
        genome_path = genome_sp
        #打开参考基因组
        
        with open(genome_path) as f:
            for record in SeqIO.parse(f,"fasta"):
                record_seq = record.seq 
        # record_seq_list_init = list(record_seq)
        out_file_sp = os.path.join(path_strain_sim_n,ls) 
        if not os.path.exists(out_file_sp):
            os.makedirs(out_file_sp)
        # output_file_name =  os.path.join(out_file_sp,ls+"_0.fasta")       
        # with open(output_file_name, "w") as output_file:
        #     SeqIO.write(record, output_file, "fasta")
        # print(genome_sp)
        
        if len(dict_sp_strain_abu[index_s])>1:
            # print(len(dict_sp_strain_abu[index_s]))
            # print(ls)
            path_sp_ani = os.path.join(path_ani_this,ls+".csv") #这个是最小生成树
            df_sp_ani = pd.read_csv(path_sp_ani,header=None,sep="\t")
            print("path_sp_ani",path_sp_ani)
            df_sp_ani.columns = ["s1","s2","ani"]
            if (len(df_sp_ani)+1) !=len(dict_sp_strain_abu[index_s]):
                print(len(df_sp_ani),len(dict_sp_strain_abu[index_s]))
                print("match error")
                exit(1)
            # list_ani = df_sp_ani["ani"].to_list()
            # print(out_file_sp)
            generate_strain(out_file_sp,df_sp_ani,genome_path,ls) #输出的文件路径
            
            # for lani in list_ani:
            #     print(lani)
                
            
            # print(len(df_sp_ani))
        
            # print(path_sp_ani,"\n")
            
        else:
            
            # path_sp_ani = os.path.join(path_ani_this,ls+".csv")
            print("one strain",len(dict_sp_strain_abu[index_s]),path_sp_ani)
            output_file_name =  os.path.join(out_file_sp,ls+".fasta")       
            # with open(output_file_name, "w") as output_file:
            #     SeqIO.write(record, output_file, "fasta")
            circulator(record,output_file_name)
    
    
    