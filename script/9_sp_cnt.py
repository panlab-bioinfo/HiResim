import os
import sys
import cv2
import pickle
import numpy as np
import pandas as pd
from sklearn.cluster import MeanShift
from sklearn.cluster import estimate_bandwidth
from ctgan import CTGAN
from ctgan import load_demo
# python 9_sp_cnt.py oral
start_train = 0
env = sys.argv[1]
path_home = sys.argv[2]
path_data = os.path.join(path_home,"data")
path_file = os.path.join(path_data,"file")
# abu_thr = float(sys.argv[2]) #一般是0.001
 
path_name_root = os.path.join(path_data,"env/data/2_sp_cut/",env[:-2])
#count
def get_num_new(data):
    clu = data.columns
    datad =[y for y in ([data[data[clu[i]]==x][clu[i+1]].nunique() \
                         for x in data[clu[i]].unique()] for i in range(len(clu)-1))] 
    return datad

path_tax = os.path.join(path_file,"gtdb-taxonomy-table_new.tsv")

taxonomy=pd.read_csv(path_tax,sep="\t")
def get_num_new_path(path_name_root,len_n = 0, taxonomy = taxonomy):

    list_name = sorted(os.listdir(path_name_root))
    flag = True
    
    data_list = []
        
    for i in list_name[len_n:]:
        path_name = os.path.join(path_name_root,i)
        name_tmp = pd.read_csv(path_name,header=None,sep="\t")
        name_tmp.columns=["Species","read_id"]
        name_tmp_sp = name_tmp.loc[:,"Species"]
        merge_tmp = pd.merge(taxonomy,name_tmp_sp,on=["Species"])
        
        data_tmp = get_num_new(merge_tmp)
        print(data_tmp)
        if len(data_tmp[0])==1:
            data_tmp[0] = [0,data_tmp[0][0]] #加个0
        data_list.append(list(data_tmp))
        #append一直出问题，直到给他加上了list(),难道return的不是list?
        if flag:
            flag=False
            data_all = data_tmp[:]
        else:
            for j in range(len(data_all)):
                data_all[j]=data_all[j]+data_tmp[j]
                #两个列表相加，其实是合并，而非求和s
    return data_all,data_list
        
data_all,data_list = get_num_new_path(path_name_root)

data_all_train,data_list_train = get_num_new_path(path_name_root,start_train)


path_sp_cnt = os.path.join(path_data,"env/data/9_sp_cnt/",env)
if not os.path.exists(path_sp_cnt):
    os.makedirs(path_sp_cnt)
path_data_cnt = os.path.join(path_sp_cnt,"sp_cnt.pkl")
path_data_train = os.path.join(path_sp_cnt,"sp_cnt_train.pkl")
with open(path_data_cnt,"wb") as f:
    pickle.dump(data_all,f)
with open(path_data_train,"wb") as f:
    pickle.dump(data_all_train,f)


#对统计得到的data_all进行聚类模拟，只是学习真实数据，但是没模拟
path_meanshift = os.path.join(path_sp_cnt,"meanshift")
if not os.path.exists(path_meanshift):
    os.makedirs(path_meanshift)
#将meanshif的结果保存起来
def save_meanshift(data_all,tail="",path_meanshift = path_meanshift):
    level = ["d","p","c","o","f","g"]

    criteria = (cv2.TERM_CRITERIA_EPS+cv2.TERM_CRITERIA_MAX_ITER,10,1) 
    #kmeans的评价指标
    K = 2
    #k代表了这个层级下面lineage的个数
    #那么我meanshif得到的聚类数目是mk （meanshift K)
    flag =True #看看是不是界，界用kmeans
    for index_l,i in enumerate(data_all):
        path_x = os.path.join(path_meanshift,"x_"+level[index_l]+tail+".pkl")
        path_label = os.path.join(path_meanshift,"label_"+level[index_l]+tail+".pkl")
        if flag:
            x = np.array(i).reshape(-1,1)
            Z = np.array(i)
            Z = np.float32(Z)
            ret,my_label,center = cv2.kmeans(Z,K,criteria,None,10,cv2.KMEANS_PP_CENTERS)
            flag=False
        else:
            x = np.array(i).reshape(-1,1)
            n_samples=min(2000,int(len(i)))
            
            bandwidth = estimate_bandwidth(x, quantile=0.3,n_jobs=32,n_samples=n_samples)
            if bandwidth == 0:
                bandwidth = 0.2
            clustering = MeanShift(bandwidth=bandwidth, bin_seeding=True,n_jobs=32).fit(x)
            my_label = clustering.labels_
        print(path_x,path_label)
        with open(path_x,"wb") as f:
            pickle.dump(x,f)
        with open(path_label,"wb") as f:
            pickle.dump(my_label,f)
            
save_meanshift(data_all)
save_meanshift(data_all_train,"_train")
print("domain",len(data_all[0]))

#我们之前保存的只是全部的样本的cnt信息，为了验证结果，我们需要保留训练集的cnt信息



#生成属丰度，属下物种，以及物种丰度的文件
print("abu")
path_abu_root = os.path.join(path_data,"env/data/1_abundance",env[:-2])
# path_read_cnt = os.path.join("/data/huixingqi/data/env/data/1_read_cnt",env+"_read_cnt.txt")

# def get_genus(path_abundance_root,path_name_root,start_id=0,taxonomy=taxonomy):
    
#     # read_cnt = pd.read_csv(path_read_cnt,sep=" ",header=None)
#     # read_cnt.columns=["read","sample"]
#     # read_cnt_list = read_cnt.loc[:,"read"].to_list()
#     # abundance_name = sorted(os.listdir(path_abundance_root))#得到所有的rep.txt
#     list_name = sorted(os.listdir(path_name_root)) #用cut后的sp就能找到物种的信息
#     # print(len(abundance_name))
#     # cnt = 0
#     dict_g = dict() #存放属丰度的字典，如果被添加前先看看是不是在字典中，如果不在就新加入，在就求和
#     dict_g_abu = dict()#用一个字典存放属下物种 字典的值为一个多层列表，每一层代表一个样本中该属的物种丰度
#     dict_g_sp = dict()#用一个字典存放属下物种 字典的值为一个多层列表，每一层代表一个样本中该属的物种丰度
#     # sp_a_g =[]
#     for i_sp in list_name[start_id:]:
#         i_abu = i_sp.replace("_cutsp","_abundance")
        
#         path_abu = os.path.join(path_abundance_root,i_abu)
#         df_abu = pd.read_csv(path_abu,sep="\t",header=None)
#         df_abu.columns =["Species","abu"]
#         path_sp = os.path.join(path_name_root,i_sp)
#         df_sp = pd.read_csv(path_sp,sep="\t",header=None)
#         df_sp.columns =["Species","cnt"]
#         list_sp = df_sp["Species"].to_list()
#         mask = df_abu["Species"].isin(list_sp)
#         rep_drop = df_abu[mask]
#         print(len(rep_drop))
#         df_merge = pd.merge(taxonomy,rep_drop,on=["Species"])
#         # df_merge = pd.merge(df_abu,df_sp,on=["Species"])
#         df_genus = df_merge.groupby("Genus")["abu"].sum()
#         #之前的方法是重新判断abu中>0.01的然后计算sp的abu，其实前面已经在生成2_sp_cut的时候就已经生成了
#         #所以只需要将之前的和现在的取个交集即可
# #         # print(path_abu)
# #         read_cnt_flag = int(read_cnt_list[cnt]*abu_thr)
# #         # print(cnt,"read_cnt_flag",read_cnt_flag)
        
# #         # if read_cnt_flag>=2000:
# #         if read_cnt_flag>=0:
# #             rep = pd.read_csv(path_abu,sep="\t",header=None)
# #             rep.columns =["Species","abu"]
# #             mask = (rep["abu"]<=read_cnt_flag) 
# #             rep_drop = rep.drop(rep[mask].index)
            
# #             df_merge = pd.merge(taxonomy,rep_drop,on=["Species"])
# #             len_merge = len(df_merge)
# #             if len_merge>=28 and len_merge<300:
# #                 print(len_merge)
# #             # print(df_merge)
# #                 df_genus = df_merge.groupby("Genus")["abu"].sum()

#             #  Domain  Phylum   Class   Order   Family   Genus   Species    abu
#             # Bacteria  Bacteroidota  Bacteroidia  Bacteroidales  Bacteroidaceae Alloprevotella  s__Alloprevotella_sp000437675  46018 

#         sum_g_read_cnt = sum(df_genus["abu"].to_list())
#         for key,val in df_genus.items():
#             if key in dict_g.keys():
                
#                 # dict_g[key] += val/read_cnt_list[cnt]
#                 dict_g[key] += val/sum_g_read_cnt
#             else:
#                 # dict_g[key] = val/read_cnt_list[cnt]
#                 dict_g[key] = val/sum_g_read_cnt

#         for g in df_merge["Genus"].unique(): 

#             if g not in list(dict_g_sp.keys()):
#                 #把这个属所有物种都放进字典，当然有的物种可能一次都没出现，后面删除即可
#                 dict_g_sp[g] = list(taxonomy[taxonomy["Genus"]==g]["Species"].unique())
#                 dict_g_abu[g] = []
#             list_g_sp = dict_g_sp[g]
#             list_g_abu = [0 for lgs in list_g_sp]
#             for index,row in df_merge[df_merge["Genus"]==g].iterrows():
#                 index_g = list_g_sp.index(row[6])
#                 list_g_abu[index_g] = row[7]
#             list_g_abu_mean = [lga/sum(list_g_abu) for lga in list_g_abu]
#             dict_g_abu[g].append(list_g_abu_mean)
#                 # if cnt==0:
#                 #     print(list_g_abu_mean)
            
#         # cnt+=1
                    
#     return dict_g,dict_g_abu,dict_g_sp

#看看属下物种的丰度信息
def get_genus(path_abundance_root,path_name_root,start_id=0,taxonomy=taxonomy):
    
    list_name = sorted(os.listdir(path_name_root)) #用cut后的sp就能找到物种的信息
    # print(len(abundance_name))
    # cnt = 0
    dict_g = dict() #存放属丰度的字典，如果被添加前先看看是不是在字典中，如果不在就新加入，在就求和
    dict_g_abu = dict()#用一个字典存放属下物种 字典的值为一个多层列表，每一层代表一个样本中该属的物种丰度
    dict_g_sp = dict()#用一个字典存放属下物种 字典的值为一个多层列表，每一层代表一个样本中该属的物种丰度
    # sp_a_g =[]
    dict_g_sp_abu = dict() #这个是物种的丰度列表，key为属，值是一个字典，字典的key为物种，值为物种的丰度
    
    for i_sp in list_name[start_id:]:
        i_abu = i_sp.replace("_cutsp","_abundance")
        
        path_abu = os.path.join(path_abundance_root,i_abu)
        df_abu = pd.read_csv(path_abu,sep="\t",header=None)
        df_abu.columns =["Species","abu"]
        path_sp = os.path.join(path_name_root,i_sp)
        df_sp = pd.read_csv(path_sp,sep="\t",header=None)
        df_sp.columns =["Species","cnt"]
        list_sp = df_sp["Species"].to_list()
        mask = df_abu["Species"].isin(list_sp)
        rep_drop = df_abu[mask]
        # print(len(rep_drop))
        df_merge = pd.merge(taxonomy,rep_drop,on=["Species"])
        # df_sp_merge = pd.merge(taxonomy,df_sp,on=["Species"]) #abu文件夹下只有物种，这个物种和属都有了
        # print(df_merge[:2])
        # print(df_sp_merge[:2]) #这个不对，这个是物种和s1编号，和abu没关系，不要这个了
        sum_sample_abu = float(sum(df_merge["abu"].to_list()))
        # print("sum_sample_abu",sum_sample_abu)
        list_g_uni = list(df_merge["Genus"].unique())
        # print("list_g_uni",list_g_uni)
        for lgu in list_g_uni:
            if lgu not in dict_g_sp_abu.keys():
                dict_g_sp_abu[lgu] = dict() #这个记录sp的丰度字典
                
            df_this_g = df_merge[df_merge["Genus"]==lgu] #这个属对应的所有物种
            for this_sp in df_this_g["Species"]:
                if this_sp not in dict_g_sp_abu[lgu].keys():
                    dict_g_sp_abu[lgu][this_sp] = []
                
                
                this_sp_abu=[float(df_this_g[df_this_g["Species"]==this_sp]["abu"].unique()[0])/sum_sample_abu]
                dict_g_sp_abu[lgu][this_sp]+=this_sp_abu
                #这个物种的abu
                
                # this_sp_abu = df_this_sp["abu"][0]
                # print(this_sp_abu)
            # print("--")
            # print("df_this_g",df_this_g)
            
        # print(len(df_sp_merge),df_sp_merge["Genus"].nunique())
        # print(df_merge)
        # df_merge = pd.merge(df_abu,df_sp,on=["Species"])
        df_genus = df_merge.groupby("Genus")["abu"].sum()
        df_genus_df = pd.DataFrame(df_genus)
        df_genus_df = df_genus_df.reset_index()
        # print(df_genus_df)
        # print(len(df_genus_df))
        #之前的方法是重新判断abu中>0.01的然后计算sp的abu，其实前面已经在生成2_sp_cut的时候就已经生成了
        #所以只需要将之前的和现在的取个交集即可

        sum_g_read_cnt = float(sum(df_genus_df["abu"].to_list()))
        # print(sum_g_read_cnt)
        for key,val in df_genus.items():
            if key in dict_g.keys():
                
                # dict_g[key] += val/read_cnt_list[cnt]
                dict_g[key] += val/sum_g_read_cnt
            else:
                # dict_g[key] = val/read_cnt_list[cnt]
                dict_g[key] = val/sum_g_read_cnt

        for g in df_merge["Genus"].unique(): 

            if g not in list(dict_g_sp.keys()):
                #把这个属所有物种都放进字典，当然有的物种可能一次都没出现，后面删除即可
                dict_g_sp[g] = list(taxonomy[taxonomy["Genus"]==g]["Species"].unique())
                dict_g_abu[g] = []
            list_g_sp = dict_g_sp[g]
            list_g_abu = [0 for lgs in list_g_sp]
            for index,row in df_merge[df_merge["Genus"]==g].iterrows():
                index_g = list_g_sp.index(row[6])
                list_g_abu[index_g] = row[7]
            list_g_abu_mean = [lga/sum(list_g_abu) for lga in list_g_abu]
            dict_g_abu[g].append(list_g_abu_mean)               
    return dict_g,dict_g_abu,dict_g_sp,dict_g_sp_abu
# print(dict_g_sp_abu)


dict_g,dict_g_abu,dict_g_sp,dict_g_sp_abu = get_genus(path_abu_root,path_name_root)
dict_g_train,dict_g_abu_train,dict_g_sp_train,dict_g_sp_abu_train = get_genus(path_abu_root,path_name_root,start_train) 

def remove_zero_columns(matrix): #去掉一次也没出现的sp?
    col_matrix = [list(col) for col in zip(*matrix) if any(col)]
    res = [list(col) for col in zip(*col_matrix)]
    removed_columns = []
    for i, col in enumerate(zip(*matrix)):
        if not any(col):
            removed_columns.append(i)
    return res,removed_columns

key_list = [] #属名字
val_list = [] #属丰度
for key,val in dict_g.items(): 
    # print(key,val)
    key_list.append(key)
    val_list.append(val)
sum_val = sum(val_list) #获得val的和
prob_list =[i/sum_val for i in val_list]#归一化
#训练集
key_list_train = [] #属名字
val_list_train = [] #属丰度
for key,val in dict_g_train.items(): 
    # print(key,val)
    key_list_train.append(key)
    val_list_train.append(val)
sum_val_train = sum(val_list_train) #获得val的和
prob_list_train =[i/sum_val_train for i in val_list_train]#归一化



dict_new_sp = dict()
dict_new_abu = dict()
for key,val in dict_g_abu.items():
    res,remove = remove_zero_columns(val)
    dict_new_sp[key] = [dict_g_sp[key][i] for i in range(len(dict_g_sp[key])) if i not in remove]
    dict_new_abu[key] = res
    
#训练集
dict_new_sp_train = dict()
dict_new_abu_train = dict()
for key,val in dict_g_abu_train.items():
    res,remove = remove_zero_columns(val)
    dict_new_sp_train[key] = [dict_g_sp_train[key][i] for i in range(len(dict_g_sp_train[key])) if i not in remove]
    dict_new_abu_train[key] = res
    
#将结果保存
path_genus = os.path.join(path_sp_cnt,"genus")
if not os.path.exists(path_genus):
    os.makedirs(path_genus)

path_genus_train = os.path.join(path_sp_cnt,"genus_train")
if not os.path.exists(path_genus_train):
    os.makedirs(path_genus_train)
    
#属名字
path_key_list = os.path.join(path_genus,"g_name_list.pkl")
with open(path_key_list,"wb") as f:
    pickle.dump(key_list,f)

#属丰度
path_val_list = os.path.join(path_genus,"g_abu_list.pkl")
with open(path_val_list,"wb") as f:
    pickle.dump(prob_list,f)
    
#属下物种
path_dict_new_sp = os.path.join(path_genus,"g_sp.pkl")
with open(path_dict_new_sp,"wb") as f:
    pickle.dump(dict_new_sp,f)
#属下物种丰度  
path_dict_new_abu = os.path.join(path_genus,"g_sp_abu.pkl")
with open(path_dict_new_abu,"wb") as f:
    pickle.dump(dict_new_abu,f)
    
#属下物种丰度  全局 
path_dict_new_abu_all = os.path.join(path_genus,"g_sp_abu_all.pkl")
with open(path_dict_new_abu_all,"wb") as f:
    pickle.dump(dict_g_sp_abu,f)
    
#训练集
path_key_list_train = os.path.join(path_genus_train,"g_name_list_train.pkl")
with open(path_key_list_train,"wb") as f:
    pickle.dump(key_list_train,f)

#属丰度
path_val_list_train = os.path.join(path_genus_train,"g_abu_list_train.pkl")
with open(path_val_list_train,"wb") as f:
    pickle.dump(prob_list_train,f)
    
#属下物种
path_dict_new_sp_train = os.path.join(path_genus_train,"g_sp_train.pkl")
with open(path_dict_new_sp_train,"wb") as f:
    pickle.dump(dict_new_sp_train,f)
#属下物种丰度  ，局部丰度
path_dict_new_abu_train = os.path.join(path_genus_train,"g_sp_abu_train.pkl")
with open(path_dict_new_abu_train,"wb") as f:
    pickle.dump(dict_new_abu_train,f)
    
#属下物种丰度  全局丰度    
path_dict_new_abu_all_train = os.path.join(path_genus_train,"g_sp_abu_all_train.pkl")
with open(path_dict_new_abu_all_train,"wb") as f:
    pickle.dump(dict_g_sp_abu_train,f)

#保存ctgan的结果
path_ctgan = os.path.join(path_sp_cnt,"ctgan")
if not os.path.exists(path_ctgan):
    os.makedirs(path_ctgan)
#ctgan 训练集
path_ctgan_train = os.path.join(path_sp_cnt,"ctgan_train")
if not os.path.exists(path_ctgan_train):
    os.makedirs(path_ctgan_train)



    
#将训练好的ctgan模型保存
ctgan = CTGAN(epochs=20)
for key,val in dict_new_abu.items():
    #val是二维列表，vl是某个样本的物种分别
    #sum(num != 0 for num in vl)是某样本不为0的物种数
    # print(key)
    #val_new应该是把总的sp个数考虑在内了
    val_new = [vl+[sum(num != 0 for num in vl)] for vl in val ]
    discrete_columns = [str(i) for i in range(len(val_new[0]))]
    df_str = pd.DataFrame(val_new)
    df_str.columns=discrete_columns
    name_save = os.path.join(path_ctgan,key+".pkl")
    ctgan.fit(df_str,discrete_columns)
    ctgan.save(name_save)

ctgan = CTGAN(epochs=20)
for key,val in dict_new_abu_train.items():
    #val是二维列表，vl是某个样本的物种分别
    #sum(num != 0 for num in vl)是某样本不为0的物种数
    # print(key)
    #val_new应该是把总的sp个数考虑在内了
    val_new = [vl+[sum(num != 0 for num in vl)] for vl in val ]
    discrete_columns = [str(i) for i in range(len(val_new[0]))]
    df_str = pd.DataFrame(val_new)
    df_str.columns=discrete_columns
    name_save = os.path.join(path_ctgan_train,key+".pkl")
    ctgan.fit(df_str,discrete_columns)
    ctgan.save(name_save)
    
    

