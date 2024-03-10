#批量跑kraken
# nohup python 1_kraken.py marine 0 10 2&>log/run_1_kraken.log &
import pandas as pd
from pandas import read_csv
import os
import subprocess
import sys
#首先要判断是什么模式，有的是双端，有的是单端
#需要输入文件夹和输出文件夹,输出文件根据输入生成，不用指定了
env = sys.argv[1]
start = int(sys.argv[2])
end = int(sys.argv[3])
mood = int(sys.argv[4]) #是什么模式，mood==1是单端，2是双端
read_path = sys.argv[5]
path_db = sys.argv[6]
path_home= sys.argv[7]
path_data = os.path.join(path_home,"data")
path_res = os.path.join(path_data,"env/res")
# abs_path = os.path.abspath(sys.argv[0])
# abs_path_dir = os.path.dirname(abs_path)
abs_path_dir = os.path.join(path_home,"script")

# in_file = os.path.join("/data/huixingqi/data/env",env)
# mood = 2 #默认是双端哈，不默认了，反正也不能省略

in_file = read_path
# in_file = os.path.join(read_path,env)
def get_report_1(in_file,out_file,path_db = path_db,abs_path_dir = abs_path_dir):
    #单端
    # path_db = "/dev/shm/db/"
    list_srr = sorted(os.listdir(in_file))
    for f in list_srr:
        srr = f.split(".")[0]
        res = os.path.join(out_file,srr+".out")
        report=os.path.join(out_file,srr+"_rep.txt")
        path_f = os.path.join(in_file,f)
        command = "/data/huixingqi/software/1_tools/kraken2/kraken2/kraken2 --db "+ path_db+" --memory-mapping "+path_f+" --use-names"+" --output "+res+" --report "+report + " --threads 32 --use-mpa-style"
        print(command)
        subprocess.run(command,
                       shell=True,
                       cwd = abs_path_dir)

def get_report_2(in_file,out_file,path_db = path_db,abs_path_dir = abs_path_dir):
    # path_db = "/dev/shm/db/"
    list_srr_all = sorted(os.listdir(in_file))
    #双端的问题是有两个srr，_1和_2，所以要获得根name
    list_srr = [srr.replace("_1.fastq.gz","") for srr in list_srr_all if "_1" in srr]
    suffix = list_srr_all[0].split("_")[-1][1:]
    # print(suffix)
    # print(list_srr)
    for f in list_srr:#batch
        # srr = f.split("_")[0]
        srr = f.split(".")[0]
        # print(srr)
        res=os.path.join(out_file,srr+".out")
        report=os.path.join(out_file,srr+"_rep.txt")
        f1 = f+"_1"+suffix
        f2 = f+"_2"+suffix
        path_f1 = os.path.join(in_file,f1)
        path_f2 = os.path.join(in_file,f2)
        # print(path_f1,path_f2)
        
        #为了防止有的双端数据不全，判断一下
        a = os.path.exists(path_f1)
        b = os.path.exists(path_f2)
        if (a and b)==False:
            print("error file is not exists\n")
            print(a,b)
        # print(res)
        command = "/data/huixingqi/software/1_tools/kraken2/kraken2/kraken2 --db "+ path_db+" --memory-mapping --paired "+path_f1+" "+path_f2 \
                +" --use-names"+" --output " +res+" --report "+report +" --threads 32 --use-mpa-style"
        print(command)
        subprocess.run(command,
                       shell=True,
                       cwd = abs_path_dir)
    

out_file = os.path.join(path_res,"out",env)

# if os.path.exists(out_file)==False:
#     os.makedirs(out_file)
# print(in_file)
# print(out_file)
    
if mood==1:
     get_report_1(in_file,out_file)
else:
     get_report_2(in_file,out_file)

if mood == 1:
    srr_batch = sorted(os.listdir(in_file))[start:end]
    # print(srr_batch)
    
    
    srr_batch_new = [srr.split(".")[0] for srr in srr_batch]
    # print(srr_batch_new)
else:
    list_srr = sorted(os.listdir(in_file))
    # print(start,end)

    # # print(list_srr[0])
    # print(srr_batch)
    srr_batch = [srr.split(".")[0].replace("_1","") for srr in list_srr if "_1" in srr]
    srr_batch_new = srr_batch[start:end]
    # print(srr_batch_new)
    
    
# print(srr_batch_new)

df_batch = pd.DataFrame(srr_batch_new)
path_batch_root = os.path.join(path_res,"batch",env)
# if not os.path.exists(path_batch_root):
#     os.makedirs(path_batch_root)
path_batch = os.path.join(path_batch_root,str(start)+"_"+str(end)+".txt")
df_batch.to_csv(path_batch,header=None,index = False)

cmd_get_read_cnt = "bash 1_2_read_cnt_abu_newout.sh "+ env+" "+path_batch+" "+path_home
print(cmd_get_read_cnt)
subprocess.run(cmd_get_read_cnt,
                shell=True,
                cwd = abs_path_dir)
