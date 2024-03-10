#得到截断后物种的read,存放在read_python中
import numpy as np
import pandas as pd
import os
import sys
import gzip
from Bio import SeqIO
from Bio.Seq import Seq
import subprocess

env = sys.argv[1]
start = int(sys.argv[2])
end = int(sys.argv[3])
mood = int(sys.argv[4])
k_kmer = int(sys.argv[5])
path_home = sys.argv[6] # HiResim/data/
path_data = os.path.join(path_home,"data")
read_path = sys.argv[7]
path_data_env = os.path.join(path_data,"env/data") #HiResim/data/env/data
# abs_path = os.path.abspath(sys.argv[0])
# abs_path_dir = os.path.dirname(abs_path)
abs_path_dir = os.path.join(path_home,"script")

#其实就是想将sp和readid对应起来
path_kraken_out = os.path.join(path_data_env,"1_out_new",env) #存放readid和s1的对应关系

# path_data = os.path.join("/data/huixingqi/data/env",env) #存放下载好的fastq
# path_read = os.path.join(read_path,env) #存放下载好的fastq，改成old的路径了
path_read = read_path
path_out = os.path.join(path_data_env,"3_read_python/",env) #存放每个物种的fastq文件

path_sp_cut = os.path.join(path_data_env,"2_sp_cut/",env) #存放物种和s1的对应关系

list_dir = sorted(os.listdir(path_sp_cut))
list_dir_sub = list_dir[start:end]
print(list_dir_sub)
#获得每个样本的中每个物种的read,保存在read_python中

def get_read1(list_dir_sub,abs_path_dir = abs_path_dir):
    for sp_cut in list_dir_sub:
        path_cut = os.path.join(path_sp_cut,sp_cut) #得到sp文件
        path_kraken = os.path.join(path_kraken_out,sp_cut.replace("_cutsp.txt",".out")) #kraken.out文件
        path_fastq =  os.path.join(path_out,sp_cut.replace("_cutsp.txt",""))#这个是样本srr不是物种
        
        # if not os.path.exists(path_fastq):
        #     os.makedirs(path_fastq)
        cmd_mkdir="bash 3_mkdir_read_kmer.sh "+path_fastq
        subprocess.run(cmd_mkdir,
                shell=True,
                cwd = abs_path_dir)

        sp_1 = pd.read_csv(path_cut,sep="\t",header=None,dtype=str)
        sp_1.columns=["sp","s1"]

        out1 = pd.read_csv(path_kraken,sep="\t",header=None)
        out1.columns=["read_id","s1"]

        df_dict = dict() #存放sp和s1对应关系
        for index,row in sp_1.iterrows():
            df_dict[row[1]] = row[0]

        list_sp = sp_1["s1"].tolist() #存放截断后的s1
        list_s1 = out1["s1"].tolist() #存放的每个read对应的s1

        input_1 = os.path.join(path_read,sp_cut.replace("_cutsp.txt",".fastq.gz"))

        dict_sp_1 =dict() #存放s1和read对应关系

        for i in list_sp:
            dict_sp_1[i] = []

        with gzip.open(input_1, "rt") as fastq_1:
            cnt_s1 = 0
            for record1 in SeqIO.parse(fastq_1, "fastq"):
                    if list_s1[cnt_s1] in list_sp:
                        dict_sp_1[list_s1[cnt_s1]].append(record1)
                    cnt_s1+=1

        for key in dict_sp_1.keys():
            sp_1 = df_dict[key]+".fastq"
            path_fastq_out_1 = os.path.join(path_fastq,sp_1) 
            with open(path_fastq_out_1, "w") as output_file:
                    SeqIO.write(dict_sp_1[key], output_file, "fastq")

def get_read2(list_dir_sub,abs_path_dir = abs_path_dir):
    
    for sp_cut in list_dir_sub:
        path_cut = os.path.join(path_sp_cut,sp_cut) #得到sp文件
        path_kraken = os.path.join(path_kraken_out,sp_cut.replace("_cutsp.txt",".out")) #kraken.out文件
        path_fastq =  os.path.join(path_out,sp_cut.replace("_cutsp.txt",""))
        # print(path_kraken)
        # print(path_fastq)
        # if not os.path.exists(path_fastq):
        #     os.makedirs(path_fastq)
        cmd_mkdir="bash 3_mkdir_read_kmer.sh "+path_fastq
        subprocess.run(cmd_mkdir,
                shell=True,
                cwd = abs_path_dir)
        
        # start_time = time.time()
        sp_1 = pd.read_csv(path_cut,sep="\t",header=None,dtype=str)
        sp_1.columns=["sp","s1"]

        out1 = pd.read_csv(path_kraken,sep="\t",header=None)
        out1.columns=["read_id","s1"]

        df_dict = dict() #存放sp和s1对应关系
        for index,row in sp_1.iterrows():
            df_dict[row[1]] = row[0]

        list_sp = sp_1["s1"].tolist() #存放截断后的s1
        list_s1 = out1["s1"].tolist() #存放的每个read对应的s1

        input_1 = os.path.join(path_read,sp_cut.replace("_cutsp.txt","_1.fastq.gz"))
        input_2 = os.path.join(path_read,sp_cut.replace("_cutsp.txt","_2.fastq.gz"))
        dict_sp_1 =dict() #存放s1和read对应关系
        dict_sp_2 =dict() #存放s1和read对应关系
        for i in list_sp:
            dict_sp_1[i] = []
            dict_sp_2[i] = []

        with gzip.open(input_1, "rt") as fastq_1, gzip.open(input_2, "rt") as fastq_2:
            cnt_s1 = 0
            for record1, record2 in zip(SeqIO.parse(fastq_1, "fastq"),SeqIO.parse(fastq_2, "fastq")):
                    if list_s1[cnt_s1] in list_sp:
                        dict_sp_1[list_s1[cnt_s1]].append(record1)
                        dict_sp_2[list_s1[cnt_s1]].append(record2)
                    cnt_s1+=1

        for key in dict_sp_1.keys():
            sp_1 = df_dict[key]+"_1.fastq"
            sp_2 = df_dict[key]+"_2.fastq"
            path_fastq_out_1 = os.path.join(path_fastq,sp_1) 
            path_fastq_out_2 = os.path.join(path_fastq,sp_2) 
            with open(path_fastq_out_1, "w") as output_file:
                    SeqIO.write(dict_sp_1[key], output_file, "fastq")

            with open(path_fastq_out_2, "w") as output_file:
                    SeqIO.write(dict_sp_2[key], output_file, "fastq")
# print(list_dir_sub)
if mood==1:
    get_read1(list_dir_sub)    
else:
    get_read2(list_dir_sub)

# k_kmer = 14
# k_kmer = 35
df_batch = pd.DataFrame(list_dir_sub)
#SRR3933087_cutsp.txt
path_batch_root = os.path.join(path_data_env,"3_read_python/batch",env)
# if not os.path.exists(path_batch_root):
#     os.makedirs(path_batch_root)
    
path_batch = os.path.join(path_batch_root,str(start)+"_"+str(end)+".txt")
df_batch.to_csv(path_batch,header=None,index = False)

cmd_get_sp_read_kmer = "bash 3_2_get_sp_read_kmer_batch.sh "+ env+" "+str(mood)+" "+path_batch+" "+str(k_kmer)+" "+path_data
subprocess.run(cmd_get_sp_read_kmer,
                shell=True,
                cwd = abs_path_dir)