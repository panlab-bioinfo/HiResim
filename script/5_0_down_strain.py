#下载当前sp的所有strain
import sys
import pandas as pd
import os
import subprocess
path_home = sys.argv[1]
env = sys.argv[2]
path_data = os.path.join(path_home,"data")
path_script = os.path.join(path_home,"script")
path_tools = os.path.join(path_home,"tools")
path_strain_genome = os.path.join(path_tools,"strain_genome")

#已经下载的sp,新的环境遇到这些物种就不需要重新下载了
#不用已经下载的sp了，不是很好用,还是用已经下载的strain
# path_used_sp = os.path.join(path_tools,"sp_down.csv")
# df_used_sp = pd.read_csv(path_used_sp)
# print(df_used_sp)
list_genome_down = os.listdir(path_strain_genome) #已经下好的df


# 这个环境下所有的物种
path_4 = os.path.join(path_data,"env/data/4_sp_strain_path_sample",env)
path_5 = os.path.join(path_data,"env/data/5_strain_unikmer",env)
list_sp = os.listdir(path_4)
df_sp = pd.DataFrame(list_sp,columns = ["sp_txt"])
df_sp["sp"] = df_sp["sp_txt"].apply(lambda x:x.replace(".txt",""))
print(len(df_sp))
print(df_sp[:1])

# 所有strian的地址
path_strain = os.path.join(path_tools,"strain_path.csv")
df_strain = pd.read_csv(path_strain,header=None,sep="\t")
df_strain.columns = ["sp","acns","path","genome"]
print(df_strain[:1])


df_merge = pd.merge(df_sp,df_strain,on = ["sp"])
print(len(df_merge))
mask_down = df_merge["genome"].isin(list_genome_down)
df_merge_nodow = df_merge[~mask_down]
df_no_none = df_merge_nodow.dropna() #这是所有要下载的文件
path_tmp_strain = os.path.join(path_5,"strain_download.txt")
df_no_none_new = df_no_none[["sp","sp_txt","acns","path","genome"]]
print(path_strain_genome)

if len(df_no_none_new)!=0:
    print(len(df_no_none_new))
    print(df_no_none_new[:1])
    print(path_tmp_strain)
    df_no_none_new.to_csv(path_tmp_strain,header=None,sep="\t",index=False)
    path_down = os.path.join(path_home,"script/lib/down_strains.py")
    cmd_down = "python "+path_down+" "+ path_tmp_strain+" "+path_home+" "+path_strain_genome
# print(cmd_down)
    try:
        result = subprocess.run(cmd_down,shell=True,cwd="./")
        result.check_returncode()
    except subprocess.CalledProcessError as e:
        print("download strain genomes error")
        sys.exit(1)


