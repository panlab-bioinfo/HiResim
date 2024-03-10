#调用pbsinm生成strain

import pandas as pd
import numpy as np
import os
import subprocess
from lib.generate_genus_and_sp import *
import sys
from Bio import SeqIO
from Bio.Seq import Seq
import pickle
# import multiprocessing
import concurrent.futures

# input_path = sys.argv[1] #输入路径 /data/huixingqi/sim_meta/HiResim/sim_oral
# input_abu = "/data/huixingqi/data/env/data/10_sim/oral_wgs21/sp_strain_abu"
# input_path = "/data/huixingqi/sim_meta/HiResim/sim_oral"
input_abu = sys.argv[1]
input_path = sys.argv[2]
depth_init = float(sys.argv[3])
# num_para = int(sys.argv[4])

# num_para = False

# if num_para == 1:
#     num_para = True

    

flag_ccs = "hifi"
# path_abs = os.path.abspath(sys.argv[0])
# path_dir = os.path.dirname(path_abs)
# print(path_dir)
# if "ONT" in model:
#     flag_ccs = "ont"

# def sim_read_pbsim3(genome,depth,method,model,sim_out,path_abs,sample_n=1,flag_ccs="hifi",l_min=100,l_max=1000000):
#     #用pbsim3模拟测序数据
#     # method有两种，qshmm 和 errhmm
#     # cmd_pbsim_head = "bash run_pbsim3.sh "+sim_out+" "+method+" "+model+" "+path_pbsim
#     pass_sum = "10"
#     cmd_pbsim_head = "bash run_pbsim3.sh "+sim_out+" "+method+" "+model #因为pbsim内置到环境中了，所以不需要指定pbsim路径了
#     # print(cmd_pbsim_head)
#     my_pre_head="sim_"
#     cnt_sim=1
#     # print(genome,depth)
#     print(len(genome),len(depth))
#     for my_g,my_d in zip(genome,depth):
#         print(my_g,my_d)
#         my_pre=my_pre_head+str(cnt_sim)
#         cnt_sim+=1
#         cmd_pbsim = cmd_pbsim_head+ " "+str(my_d)+" "+my_g+" "+my_pre+" "+pass_sum+" "+str(sample_n)+" "+flag_ccs
#         path_shell=os.path.join(path_abs,"lib")
#         # print(path_shell)
#         print(cmd_pbsim)
#         subprocess.run(cmd_pbsim,shell=True,cwd=path_shell)
def runpbsim(cmd,path_shell):
    subprocess.run(cmd, shell=True, cwd=path_shell)

def run_pbsim_parallel(commands, path_shell):
    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = [executor.submit(runpbsim, cmd, path_shell) for cmd in commands]
        try:
            concurrent.futures.wait(futures)
        except KeyboardInterrupt:
            for future in futures:
                future.cancel()
            executor.shutdown(wait=False)
            sys.exit("pbsim terminated by user")
            
            
    
def sim_read_pbsim3(genome,depth,method,model,sim_out,path_abs,sample_n=1,flag_ccs="hifi",l_min=100,l_max=1000000):
    #用pbsim3模拟测序数据
    # method有两种，qshmm 和 errhmm
    # cmd_pbsim_head = "bash run_pbsim3.sh "+sim_out+" "+method+" "+model+" "+path_pbsim
    pass_sum = "10"
    cmd_pbsim_head = "bash run_pbsim3.sh "+sim_out+" "+method+" "+model #因为pbsim内置到环境中了，所以不需要指定pbsim路径了
    # print(cmd_pbsim_head)
    my_pre_head="sim_"
    cnt_sim=1
    # print(genome,depth)
    print(len(genome),len(depth))
    # pool = multiprocessing.Pool(32)
    cmds=[]
    list_sim = []
    path_shell=os.path.join(path_abs,"lib")
    for my_g,my_d in zip(genome,depth):
        # print(my_g,my_d)
        my_pre=my_pre_head+str(cnt_sim)
        cnt_sim+=1
        cmd_pbsim = cmd_pbsim_head+ " "+str(my_d)+" "+my_g+" "+my_pre+" "+pass_sum+" "+str(sample_n)+" "+flag_ccs
        
        # print(path_shell)
        # print(cmd_pbsim)
        cmds.append(cmd_pbsim)
        list_sim.append(my_pre)
        
    
    grouped_cmds = [cmds[i:i+20] for i in range(0, len(cmds), 20)]
    # grouped_sim = [list_sim[i:i+20] for i in range(0, len(list_sim), 20)]
    
    for group in grouped_cmds:
        # print(group)
        run_pbsim_parallel(group, path_shell)



# list_sp = sorted([ip for ip in  os.listdir(input_path) if ip.endswith(".txt")]) 
list_sp = sorted([ip for ip in  os.listdir(input_abu) if ip.endswith(".pkl")]) 
# print(list_sp)
# path_df_sp = 
for index_lp, lp in enumerate(list_sp[:1]):
    path_abu = os.path.join(input_abu,lp)
    with open(path_abu,"rb") as f:
        strain_abu = pickle.load(f)
    
    print(path_abu)
    if len(list_sp)>1:
        lp = "sim_strain"+str(index_lp)+".txt"
    
    path_sp = os.path.join(input_path,lp)
    path_strain = os.path.join(input_path,"sim_strain","sim"+str(index_lp))
    print(path_strain)
    print(path_sp)
    df_strain = pd.read_csv(path_sp,sep="\t")
    # print(df_strain[:1])
    list_abu = df_strain["abu_new"].to_list()
    # print(sum(list_abu)) #sum为1
    # print(df_strain["sp"].nunique(),len(strain_abu))
    index_sp = 0
    
    list_abu_strain = []
    dict_abu_strain = dict()
    #这一段其实可以用列表解析实现
    for sp in df_strain["sp"].unique():
        sp_abu = list_abu[index_sp] # 物种的丰度
        # print(sp_abu)
        path_sp_strain = os.path.join(path_strain,sp)
        list_strain = sorted(os.listdir(path_sp_strain))
        # print(path_sp_strain)
        for index_sabu,sabu in enumerate(strain_abu[index_sp]): #stain abu
            sabu_new = sabu/sum(strain_abu[index_sp])
            list_abu_strain.append(sp_abu*sabu_new)
            dict_abu_strain[list_strain[index_sabu]] = sp_abu*sabu_new
        index_sp+=1
    # print(sum(list_abu_strain))
    
    min_abu = min(list_abu_strain)
    # depth_init = 0.1
    list_depth_strain = []
    
    # index_strain = 0
    list_strain_choice = []
    for sp in df_strain["sp"].unique(): 
        path_sp_strain = os.path.join(path_strain,sp)
        # print(path_sp_strain)
        list_strain = sorted(os.listdir(path_sp_strain))
        for lstrain in list_strain:
            path_genome = os.path.join(path_sp_strain,lstrain)
            # print(path_genome)
            depth_genome = dict_abu_strain[lstrain]/min_abu*depth_init
            list_depth_strain.append(round(depth_genome,2))
            # print(depth_genome)
            list_strain_choice+=[path_genome]
    
    
    method = "qshmm"        
    model = "/data/huixingqi/software/1_tools/pbsim3/data/QSHMM-RSII.model"
    # sim_pbsim = "/data/huixingqi/sim_meta/HiResim/sim_gut_new/pbsim"
    sim_pbsim = os.path.join(input_path,"pbsim")
    # if not os.path.exists(sim_pbsim):
    #     os.makedirs(sim_pbsim)
    path_abs = "/data/huixingqi/sim_meta/HiResim/script" 
    # pass_sum = 10
    print("---")
    path_shell = os.path.join(path_abs,"lib")
    # print(list_abu_strain)
    cmd_create_dir = "bash mk_res.sh " + sim_pbsim +" "+ str(index_lp+1)
    subprocess.run(cmd_create_dir,shell=True,cwd=path_shell)
    
    cmds_pbsim = sim_read_pbsim3(list_strain_choice,list_depth_strain,method,model,sim_pbsim,path_abs,index_lp+1,flag_ccs)
    
    path_read_sim = os.path.join(sim_pbsim,"sim_res_"+str(index_lp+1))
    print(path_read_sim)
    
    # os.rmtree()
    cmd_concat = "bash concat.sh " + sim_pbsim +" "+ str(index_lp+1)
    subprocess.run(cmd_concat,shell=True,cwd=path_shell)
    
    cmd_rm = "bash rm_tmp.sh "+path_read_sim
    print(cmd_rm)
    
    subprocess.run(cmd_rm,shell=True,cwd=path_shell)
    
    
    
    
    
    
