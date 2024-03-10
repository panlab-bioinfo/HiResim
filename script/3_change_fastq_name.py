#因为fastq的名字和out rep abu等文件对应关系不是很明确，稍微修改一下
#为了不影响old
#跑完old之后再执行这个程序
import pandas as pd
from pandas import read_csv
import os
import subprocess
import sys
env = sys.argv[1]
mood = int(sys.argv[2])
read_path = sys.argv[3]
# env = "rhizosphere"

in_file = os.path.join(read_path,env)


# mood = 2 #默认是双端哈，不默认了，反正也不能省略
# mood = 1 #是什么模式，mood==1是单端，2是双端
# start = 0



def get_report_1(in_file):
    #单端
    
    list_srr = sorted(os.listdir(in_file))
    for f in list_srr:
        srr = f.split(".")[0]
        path_f = os.path.join(in_file,f)
        path_f_new = os.path.join(in_file,srr+".fastq.gz")
        command = "mv " +path_f +" "+path_f_new
        print(command)
        subprocess.run(command,
                       shell=True,
                       cwd = in_file)

def get_report_2(in_file):
    list_srr_all = sorted(os.listdir(in_file))
    #双端的问题是有两个srr，_1和_2，所以要获得根name
    list_srr = [srr.replace("_1.fastq.gz","") for srr in list_srr_all if "_1" in srr]
    suffix = list_srr_all[0].split("_")[-1][1:]
    # print(suffix)
    # print(list_srr)
    for f in list_srr:#batch
        # srr = f.split("_")[0]
        srr = f.split(".")[0]
        f1 = f+"_1"+suffix
        f2 = f+"_2"+suffix
        path_f1 = os.path.join(in_file,f1)
        path_f2 = os.path.join(in_file,f2)
        path_f1_new = os.path.join(in_file,srr+"_1.fastq.gz")
        path_f2_new = os.path.join(in_file,srr+"_2.fastq.gz")
        cmd1 = "mv " +path_f1 +" "+path_f1_new
        cmd2 = "mv " +path_f2 +" "+path_f2_new
        # print(cmd1)
        # print(cmd2)

        subprocess.run(cmd1,
                       shell=True,
                       cwd = in_file)
        subprocess.run(cmd2,
                       shell=True,
                       cwd = in_file)
        
if mood==1:
    print("single")
    get_report_1(in_file)
else:
    print("paired")
    get_report_2(in_file)