import pandas as pd
import numpy as np
import os
import pickle
import random
import math

from Bio import SeqIO
from Bio.Seq import Seq
# import torch
# cuda_device = torch.device("cuda:1")
# torch.cuda.set_device(1)
# current_device = torch.cuda.current_device()
# print("Currently using GPU:", current_device)

def circulator(fasta, output_file_name):
    # fasta=SeqIO.parse(input_file,"fasta")
    # print(fasta.seq)
    fasta_seq = fasta.seq
    new_seq = str(fasta_seq)+str(fasta_seq[:30000])
    # print("new_seq",new_seq)
    # print(fasta)
    fasta.seq = Seq(new_seq)
    with open(output_file_name,'w') as file:
        SeqIO.write(fasta, file, "fasta")

def get_sp_cnt_new(path_meanshift,tail=""):
    level = ["d","p","c","o","f","g"]
    list_generate = []
    for l in level:
        
        path_x = os.path.join(path_meanshift,"x_"+l+tail+".pkl")
        path_label = os.path.join(path_meanshift,"label_"+l+tail+".pkl")
       
        with open(path_x,"rb") as f:
            x = pickle.load(f)
        with open(path_label,"rb") as f:
            label = pickle.load(f)
        if l=="d":
            d1 = random.choice(x[label==0])
            d2 = random.choice(x[label==1])

            domain=[d1,d2]
            # domain_new = [dm for dm_1 in domain for dm in dm_1]
            list_generate.append(domain)
  
            K = sum(domain)
        else:
            choice_tmp = []
            #存放mk-e+1个类,其中e各类因为里面数太少,被合并为一个类,error类
            choice_e =[]
            sum_label_cnt = len(x)
            mk = len(np.unique(label))
           
            for e in range(mk): #e是label 的id,从0到mk-1个label
                x_choice = list(x[label==e].reshape(-1))
                n=(len(x_choice)/sum_label_cnt)*K
                # print("n",n)
                if n<1:
                    choice_e +=x_choice
                else:
                    choice_tmp.append(list(x_choice))
            choice_tmp.append(list(choice_e))
            n_cnt=0
            n_res = [] #保存剩下的部分作为丰度
            generate_tmp_1 =[] #保存采样的lieage
            for ct_index,ct in enumerate(choice_tmp):
                n = (len(ct)/sum_label_cnt)*K
                n_new = math.floor(n)
                #向下取整所以不会有负值
           
                n_cnt+=n_new
                n_res.append(n-n_new)
                if n_new>=1:
                    my_choice = random.sample(ct,n_new)
                    generate_tmp_1.append(list(my_choice))
            
            generate_tmp_1_new = [gt for gt_1 in generate_tmp_1 for gt in gt_1]
            
            res_cnt = K -n_cnt #看看还有多少lineage没有生产
           
            
            if res_cnt>0:
                
                label_ct = [i_ct for i_ct in range(len(choice_tmp))]
            
                pro_ct = [lc/sum(n_res) for lc in n_res]
                
                generate_label = np.random.choice(label_ct,size = res_cnt,p=pro_ct)
               
                unique_label,counts_label = np.unique(generate_label,return_counts=True)
                #choice_tmp[unique_label[index_label]]是一个聚类结果，要从该结果中抽size个数
                # generate_tmp_2 = [list(np.random.choice(choice_tmp[unique_label[index_label]], size = counts_label[index_label],replace=False)  )
                #     for index_label in range(len(unique_label))]
                generate_tmp_2 = []
                for index_label in range(len(unique_label)):
                    if counts_label[index_label] > len(choice_tmp[unique_label[index_label]]):
                        replace = True
                    else:
                        replace = False
                    # print(replace)
                    generate_tmp_2.append(list(np.random.choice(choice_tmp[unique_label[index_label]], 
                                                size=counts_label[index_label],
                                                replace=replace)))
                
                generate_tmp_2_new = [gt for gt_1 in generate_tmp_2 for gt in gt_1]
          
                generate_tmp = generate_tmp_1_new + generate_tmp_2_new
            else:
                generate_tmp = generate_tmp_1_new

            K = sum(generate_tmp)
            list_generate.append(list(generate_tmp))

        
    return  sum(list_generate[-2]),list_generate[-1],list_generate
#分别是属个数，种列表，以及全部文件
# path_model_root = "/data/huixingqi/data/env/data/9_sp_cnt/test"

#采样属
def choice_genus(key_list,prob_list,cnt_g):
    sample = np.random.choice(key_list,size =cnt_g,p=prob_list,replace=False )
    return sample

#采样物种
def choice_sp(path_model_root,sample,dict_new_sp):
    from ctgan import CTGAN
    from ctgan import load_demo
    # list_g_sp = [] #遍历样本，样本的属确定了，跑属的ctgan,得到属的物种,一个样本一组属，一个属一组物种
    #我们只生成一个样本
    list_sp = []#记录所有的属的sp,
    # print(sample)
    
    for genus in sample:
        # print(genus)
        path_model = os.path.join(path_model_root,genus+".pkl")
        # print(path_model)
        load_model = CTGAN.load(path_model)
        new_data = list(load_model.sample(1).iloc[0,:])
        new_data_sp = new_data[:-1]

        cnt_gan = new_data[-1] #gan生成的物种个数
        cnt_sum = sum(num != 0 for num in new_data[:-1])

        #判断cnt_gan和cnt_sum大小
        index_sp = [ids for ids in range(len(new_data_sp))]
        #物种下标

        while cnt_gan> len(new_data_sp):
            print("error cnt_gan")
            # print(new_data_sp)
            new_data = list(load_model.sample(1).iloc[0,:])
            new_data_sp = new_data[:-1]
            cnt_gan = new_data[-1] #gan生成的物种个数
            cnt_sum = sum(num != 0 for num in new_data[:-1])
            # index_sp = [ids for ids in range(len(new_data_sp))]

        if cnt_sum > cnt_gan:
            # print("more")
            
            prob_list = [ nds/sum(new_data_sp) for nds in  new_data_sp]
            #需要缩减到和cnt_gan一致
            new_data_choice = new_data_sp[:]
            #选择的下标

        elif cnt_sum < cnt_gan:
            # print("less")
            new_data_sp_new = new_data_sp[:]
            #需要额外补充数据
            #因为一次采样不够，多采样几次
            while cnt_sum<cnt_gan:

                new_data_add = list(load_model.sample(1).iloc[0,:])
                new_data_add_sp = new_data_add[:-1]
                new_data_sp_new = [new_data_sp_new[idx]+new_data_add_sp[idx] for idx in range(len(new_data_add_sp))]
                cnt_sum = sum(num != 0 for num in new_data_sp_new)

            # print("cnt_sum",cnt_sum)
            
            prob_list = [ nds/sum(new_data_sp_new) for nds in  new_data_sp_new]
            new_data_choice = new_data_sp_new[:]
        else:
            #二者相等的时候
            prob_list = [ nds/sum(new_data_sp) for nds in  new_data_sp]
            new_data_choice = new_data_sp[:]

        sp_list = np.random.choice(index_sp,size = int(cnt_gan),p=prob_list,replace=False)
        # print(sp_list)
        #选中的物种下标
        new_data_new = [0 for nds in new_data_sp] #这样相当于copy而不是共享
        for idx in sp_list:
            new_data_new[idx] = new_data_choice[idx]
        species = [dict_new_sp[genus][i] for i in range(len(new_data_new)) if new_data_new[i]>0]
        list_sp.append(species)
    list_sp_new = [lsn for sp_g in list_sp for lsn in sp_g]
    return list_sp_new

#得到物种的abu
def get_sp_abu(sp_list,path_abu_root,path_sp_g):
    # print()
    #输入分为三个，一个是样本中的物种，一个是属下物种丰度字典，一个是GTDB物种taxonomy信息
    path_abu_pkl = os.path.join(path_abu_root,"g_sp_abu_all.pkl")
    # huixingqi/sim_meta/sim_hre/9_sp_cnt.ipynb 得到该文件
    
    with open(path_abu_pkl,"rb") as f:
        dict_g_sp_abu = pickle.load(f)
        
    #eg /data/huixingqi/data/env/data/9_sp_cnt/gut21/genus + g_sp_abu_all.pkl
    #现在需要知道物种的属的信息，一个字典key为物种，val为属
    with open(path_sp_g,"rb") as f:
        dict_sp_g = pickle.load(f)
    list_no_genus = [] #这个是属都不存在里面
    list_no_sp = [] #这个是属存在但是物种不存在
    dict_sp_abu = dict() #存放物种采样之后的abu,key为sp，val为sp的abu
    list_sp_used = [] #使用过的sp就不要用了
    for sp in sp_list:
        
        g_sp = dict_sp_g[sp] #这个物种所在的属
        # print(sp,g_sp,len(dict_g_sp_abu.keys()))
        
        if g_sp in dict_g_sp_abu.keys():#这个物种的属在
            
            if sp in dict_g_sp_abu[g_sp].keys(): #物种存在
                dict_sp_abu[sp] = random.sample(dict_g_sp_abu[g_sp][sp],1)[0]
                list_sp_used.append(sp)
            else:#物种不存在
                list_no_sp.append(sp)
            
        else:
            list_no_genus.append(sp) #属都不存在
            
    list_sp_used_ori = list_sp_used[:] #问题在于list_sp_used是所有物种，不是这个属下的物种
    # print("dict_sp_abu1",dict_sp_abu)
    for sp in list_no_sp:
        g_sp = dict_sp_g[sp] #这个物种所在的属
        
        
        list_g_sp = list(dict_g_sp_abu[g_sp].keys())
        # if len(dict_sp_abu)==0:
        #     print("000",g_sp,sp,len(list_g_sp),len(list_no_sp))
        
        list_g_sp_noused = list(set(list_g_sp)-set(list_sp_used))
        #有两个异常情况
        # 1 这个样本中的物种包括了全部样本中的物种，这样就没有noused物种
        
        if len(list_g_sp_noused)==0:
            if len(list(set(list_g_sp)-set(list_sp_used_ori)))==0: #
                list_sp_used = [] #从已经用过的里面生成
            else:
                list_sp_used = list_sp_used_ori[:] #重新赋予new
            list_g_sp_noused = list(set(list_g_sp)-set(list_sp_used))
        # print(len(list_g_sp_noused),len(list_g_sp),len(list_sp_used_ori),g_sp,list_sp_used_ori)   
        sp_alt = random.sample(list_g_sp_noused,1)[0]
        list_sp_used.append(sp_alt)
        # print("sp,sp_alt",sp,sp_alt,len(list_g_sp_noused))
        dict_sp_abu[sp] = random.sample(dict_g_sp_abu[g_sp][sp_alt],1)[0]
    # print("dict_sp_abu2",dict_sp_abu)
    
    
    if len(dict_sp_abu)==0:
        # print("000")
        # print("1211",list_no_genus)
        #一个物种都没找到
        for sp in list_no_genus: 
            g_sp = random.sample(list(dict_g_sp_abu.keys()),1)[0]
            # print(dict_g_sp_abu[g_sp])
            sp_alt = random.sample(list(dict_g_sp_abu[g_sp].keys()),1)[0]
            dict_sp_abu[sp] = random.sample(dict_g_sp_abu[g_sp][sp_alt],1)[0]
        
    else:
        list_mean = np.mean(list(dict_sp_abu.values()))
        for sp in list_no_genus:
            # print(list(dict_sp_abu.values()),"\n")
            # dict_sp_abu[sp] = np.mean(list(dict_sp_abu.values()))
            dict_sp_abu[sp] = list_mean
    # print()
    return dict_sp_abu
            
        
    

#生成strain个数和丰度
def strain_cnt_abu(sample,df_strain_cnt,df_strain_abu):
    cnt_sp_tmp = 0 #物种个数
    list_sp_cnt = []
    list_sp_abu = []
    list_g_sp_sample =[]
    for sp in sample:
        # print("sp",sp)
        list_sp_strian_cnt = df_strain_cnt[df_strain_cnt["sp"]==sp]["cnt"].to_list()
        if len(list_sp_strian_cnt)!=0:#把不满足条件的sp去掉了
            # print("not 0")
            cnt_sp_tmp+=1
            sp_cnt = np.random.choice(list_sp_strian_cnt)
            # print(sp_cnt)
            list_sp_cnt.append(sp_cnt)
            list_g_sp_sample.append(sp)
            df_abu_sp = df_strain_abu[df_strain_abu["sp"]==sp]["abu"].to_list()
            sp_abu = list(np.random.choice(df_abu_sp,size=sp_cnt))
            list_sp_abu.append(list(sp_abu))
        else:
            # print("sp",sp)
            # print("0")
            list_sp_cnt.append(1)
            list_sp_abu.append([1])
            #如果这个物种有strain，那么就可以抽丰度，大不了抽的都一样
    return list_sp_cnt,list_sp_abu

#得到每个sp的平均strain个数
def get_sp_cnt_mean(df_strain_cnt):
    df_sp_cnt_mean = pd.DataFrame(columns=["sp","mean_cnt"])
    for sp in df_strain_cnt["sp"].unique():
        list_sp = list(df_strain_cnt[df_strain_cnt["sp"]==sp]["cnt"])
        mean_sp_cnt = np.mean(list_sp)#就是要用mean
        # mean_sp_cnt = random.sample(list_sp)
        value = {"sp":[sp],"mean_cnt":[mean_sp_cnt]}
        df_sp_cnt_mean = pd.concat([df_sp_cnt_mean,pd.DataFrame(value)],ignore_index=True)
    return df_sp_cnt_mean

#找到最适合某个没有ani的物种的相近物种的函数
def find_alt_sp(cnt_all,df_sp_cnt_mean):
    min_dis = 1000 #记录mean_cnt和cnt_all的差距，我们要找一个最小的差距
    flag_key = True
    # list_index = list(range(len(df_sp_cnt_mean)))
    list_index = list(df_sp_cnt_mean.index)
    random.shuffle(list_index)
    for key in list_index:
        if flag_key:
            best_key = key
            flag_key = False
        val1 = df_sp_cnt_mean.loc[key][1]
        dis = abs(val1-cnt_all)
        # print(val1)
        # print(val[0],val[1])
        if dis< min_dis:
            min_dis = dis
            best_key = key
    return df_sp_cnt_mean.loc[best_key][0]
#得到每个物种ani文件的地址
def get_sp_ani_path(path_tree_root):
    df_sp_sample_ani = pd.DataFrame(columns=["sp","srr","path"]) #记录每个物种在srr的出现情况
    list_srr_tree = sorted(os.listdir(path_tree_root))
    for srr in list_srr_tree:
        path_srr = os.path.join(path_tree_root,srr)
        list_sp = sorted(os.listdir(path_srr))
        for sp in list_sp:
            # print(sp)
            path_sp = os.path.join(path_srr,sp)
            if os.path.getsize(path_sp)!=0:
                sp_name = sp.replace(".txt","")
                values = {"sp":[sp_name],"srr":[srr],"path":[path_sp]}
                df_sp_sample_ani = pd.concat([df_sp_sample_ani,pd.DataFrame(values)],axis=0,ignore_index=True)
    return df_sp_sample_ani

#把生成的ani文件保存
def generate_strain_ani_new(path_strain,path_tree,sp,tree_path_list,cnt_all):
    #sim是第几次模拟，用于保存模拟结果
    #tree_path_list是当前输入的tree的地址列表，用于遍历并生成某个物种的所有ani
    #cnt_all是总共需要的strain个数
    cnt_now = 0    
    list_dis_all = []
    tree_path_list_ori = tree_path_list.copy()
    dict_tree = dict() #保存生成树
    record_strain = 0
    df_ani = pd.DataFrame() #这个是模拟群落的ani df
    list_set_index = list(range(cnt_all+1))
    list_set_index_ori = list_set_index[:]
    flag_ani = False
    list_index_used = [] #使用过的index
    while (cnt_now<cnt_all):
        # print("---cnt_now",cnt_now,"cnt_all---",cnt_all,"||||")
        if len(tree_path_list)==0:
            tree_path_list = tree_path_list_ori.copy()
        tree_path = np.random.choice(tree_path_list)
        tree_path_list = list(set(tree_path_list) - set([tree_path])) 
        # print(tree_path,"||||")
        #df_tree存放了一个样本的树
        df_tree = pd.read_csv(tree_path,header=None,sep="\t")
        df_tree.columns = ["s1","s2","dis"] #存放的是strain1和strain2之间的距离
        list_s1 = df_tree["s1"].to_list()
        list_s2 = df_tree["s2"].to_list()
        set_s = set(list_s1) | set(list_s2)
        list_s_all = list(set_s)
        
        #对于一个新的到的df_tree 直接给他换index,问题就是index超了怎么办，这样本来0-19，结果index出现了22
        #最后只需要所有超的index进行下降即可
        
        if flag_ani: #已经有df了
            # strain_s = np.random.choice(list(range(len(df_tree)+1))) #这个是这棵树中的节点
            strain_s = np.random.choice(list_s_all) #这个是这棵树中的节点
            choice_s = np.random.choice(list_index_used)#这个选一个已经有的节点
            
            # dict_csv={strain_s:choice_s}
            list_set_no_s = list(set(list_s_all) - set([strain_s])) #去掉拼接节点
            if len(list_set_index) < len(list_set_no_s): #这个df超了，说明到这个df就结束了
                cnt_sub = len(list_set_no_s)-len(list_set_index)
                list_set_index_new = list_set_index+[list_set_index[-1]+cs+1 for cs in range(cnt_sub)]
                list_set_index = list_set_index_new[:] #扩展一下index
            
            if len(list_set_index) < len(list_set_no_s):
                print("len(list_set_index) < len(list_set_no_s)",len(list_set_index) ,len(list_set_no_s))
            dict_csv = {lis: list_set_index[lis_index] for lis_index,lis in enumerate(list_set_no_s)}
            list_set_index = list_set_index[(len(list_set_no_s)):]
            list_index_used = list(set(list_set_index_ori)-set(list_set_index))

            dict_csv[strain_s] = choice_s
            #dict_csv表示了名字对应关系，现在开始改名
            try:
                df_tree["s1_new"] = df_tree["s1"].apply(lambda x:dict_csv[x])
                df_tree["s2_new"] = df_tree["s2"].apply(lambda x:dict_csv[x])
            except:
                print("dict_csv",dict_csv)
                print("list_set_no_s",list_set_no_s,strain_s)
                print(tree_path)
                exit(1)
        else:
            # strain_s = np.random.choice(list(range(len(df_tree)+1))) #这个是这棵树中的节点
            strain_s = np.random.choice(list_s_all) #这个是这棵树中的节点
            #还是要选一个strian_s ，作为这棵树的开始
            list_set_no_s = list_s_all[:] #没有拼接节点
            if len(list_set_index) < len(list_set_no_s): #这个df超了，说明到这个df就结束了
                cnt_sub = len(list_set_no_s)-len(list_set_index)
                list_set_index_new = list_set_index+[list_set_index[-1]+cs+1 for cs in range(cnt_sub)]
                list_set_index = list_set_index_new[:] #扩展一下index
            if len(list_set_index) < len(list_set_no_s):
                print("len(list_set_index) < len(list_set_no_s)",len(list_set_index) ,len(list_set_no_s))
            dict_csv = {lis: list_set_index[lis_index] for lis_index,lis in enumerate(list_set_no_s)}
            list_set_index = list_set_index[(len(list_set_no_s)):]
            list_index_used = list(set(list_set_index_ori)-set(list_set_index))
            
            #dict_csv表示了名字对应关系，现在开始改名
            try:
                df_tree["s1_new"] = df_tree["s1"].apply(lambda x:dict_csv[x])
                df_tree["s2_new"] = df_tree["s2"].apply(lambda x:dict_csv[x])
                #这里可能会报错，导致s2_new没有
                
            except:
                print("dict_csv",dict_csv)
                print("list_set_no_s",list_set_no_s,strain_s)
                print(tree_path)
                print(tree_path)
                exit(1)
            
        cnt_add = cnt_now+len(df_tree)
        flag_ani = True
        if cnt_add <= cnt_all:# 这里面每一个都能对应到list_index_set
            cnt_now = cnt_add
            list_dis_all+=df_tree["dis"].to_list()
            # 这里面每一个都能对应到list_index_set
            # dict_csv = {lis: list_set_index[lis] for lis in range(len(df_tree))}
            # list_set_index = list_set_index[len(df_tree):]
            # df_csv_s = pd.DataFrame(list(dict_csv.items()),columns=["s1","s2"])
            # df_csv_d = df_tree["dis"]
            # df_csv = pd.concat([df_csv_s,df_csv_d],axis=1)
            try:
                df_csv = df_tree[["s1_new","s2_new","dis"]]
            except:
                print("df_tree",df_tree)
                exit(1)
                
            df_ani = pd.concat([df_ani,df_csv])
            df_ani = df_ani.reset_index(drop=True)
            
        elif cnt_add >cnt_all: #这里选一部分就够用了
            list_s1 = df_tree["s1"].to_list()
            list_s2 = df_tree["s2"].to_list()
            set_s = set(list_s1) | set(list_s2)
            list_s = list(set_s)
            #变成set是为了取并集，变成list是为了后续处理方便
            list_s_new = list_s[:]
            # list_s_pre = [] #上一个节点的邻居
            flag_pre = False # 一开始没有邻居
            flag_strain_s = False
            while(cnt_now<cnt_all): #一点一点选，直到结束
                if flag_strain_s:
                    strain_s = np.random.choice(list_s_new) #先选一个节点
                flag_strain_s = True
                #第一个strain_s一定是前面的拼接节点，不然对应不起来了
                # if flag_pre:
                # list_s_nos = list(set(list_s_new) - set([strain_s])) #选了一个之后面的就不能要了
                
                #strain_s是我们随机选择的节点，从这个节点开始广度优先遍历最小生成树
                list_left = df_tree[df_tree["s1"]==strain_s]["s2"].to_list()
                list_right = df_tree[df_tree["s2"]==strain_s]["s1"].to_list()

                list_left_dis = df_tree[df_tree["s1"]==strain_s]["dis"].to_list()
                list_right_dis = df_tree[df_tree["s2"]==strain_s]["dis"].to_list()
                list_all = [] #把左相连和右相连的节点放到一起
                list_dis = [] #距离

                if list_left!=None:
                    list_all+=list_left
                    list_dis+=list_left_dis
                if list_right!=None:
                    list_all+=list_right
                    list_dis+=list_right_dis
                if len(list_all)!=0:
                    # print("list_dis",list_dis)
                    random_index = [index_list_all for index_list_all in range(len(list_all))]
                    random.shuffle(random_index) #随机打乱节点顺序

                    list_all_new = [list_all[ri] for ri in random_index]
                    list_dis_new = [list_dis[ri] for ri in random_index] # 这个是当前节点的全部邻居
                else:
                    # print()
                    print("dict_csv",dict_csv)
                    print(df_tree,strain_s)
                list_s_nos =  list_all_new[:]   
                cnt_s = cnt_now+len(list_dis_new) #可能和list_s_pre有关
                if cnt_s <= cnt_all: # 这个节点不够
                    # print("cnt_s <= cnt_all")
                    cnt_now = cnt_s
                    list_dis_all+=list_dis_new
                    # list_left_s1 = list_left+[strain_s]*len(list_right)
                    # list_left_s2 = [strain_s]*len(list_left)+list_right
                    # if len(df_ani)>0: #已经有df了
                    #     choice_s = np.random.choice(list_s_new)#这个选一个已经有的节点
                    #     dict_csv = {strain_s:choice_s}
                    #     dict_csv.undate({lis: list_set_index[lis] for lis_index,lis in enumerate(list_all_new)})
                    # list_set_index = list_set_index[len(df_tree):]
                    # df_csv_s = pd.DataFrame(list(dict_csv.items()),columns=["s1","s2"])
                    # df_csv_d = df_tree["dis"]
                    # df_csv = pd.concat([df_csv_s,df_csv_d],axis=1)
                    # df_ani = pd.concat([df_ani,df_csv])
                    # df_ani = df_ani.reset_index(drop=True)
                    df_tree_s = df_tree[(df_tree['s1'] == strain_s) | (df_tree['s2'] == strain_s)]
                    #这个节点对应的dataframe
                    df_csv = df_tree_s[["s1_new","s2_new","dis"]]
                    df_ani = pd.concat([df_ani,df_csv])
                    df_ani = df_ani.reset_index(drop=True)
                    
                else: # 这个节点够了
                    len_c = cnt_all - cnt_now
                    # print("len_c",len_c,"cnt_now",cnt_now,"cnt_s",cnt_s,"cnt_all",cnt_all)
                    # print(len(list_dis_new),len_c)
                    df_tree_s_all = df_tree[(df_tree['s1'] == strain_s) | (df_tree['s2'] == strain_s)]
                    df_tree_s = df_tree_s_all[:len_c]
                    # list_dis_all+=list_dis_new[:len_c]
                    list_dis_all+=df_tree_s["dis"].to_list()
                    # dict_csv = {lis: list_set_index[lis] for lis in range(len(len_c))}
                    # list_set_index = list_set_index[len(df_tree):]
                    # df_csv_s = pd.DataFrame(list(dict_csv.items()),columns=["s1","s2"])
                    # df_csv_d = df_tree["dis"]
                    df_csv = df_tree_s[["s1_new","s2_new","dis"]]
                    df_ani = pd.concat([df_ani,df_csv])
                    df_ani = df_ani.reset_index(drop=True)
                    cnt_now = cnt_all
                df_tree_l = df_tree[df_tree["s1"]!=strain_s]
                df_tree = df_tree_l[df_tree_l["s2"]!=strain_s]
                list_s1 = df_tree["s1"].to_list()
                list_s2 = df_tree["s2"].to_list()
                set_s = set(list_s1) | set(list_s2) #并集，新df的所有节点
                #首先set_s中是 strain_s的所有非孤立点邻居
                # list_s_pre 中是上一个节点的所有邻居
                
                # list_s = list(set_s)
                # list_s_new_all = list_s_new + list_s_pre #这里面有孤岛，即只和strain_s相邻的，这种没有别的邻居了
                if flag_pre:
                    list_s_new_all = list_s_nos[:]+list_s_pre[:] #我的邻居加上一个的邻居
                else:
                    list_s_new_all = list_s_nos[:] #第一次的时候只有我的邻居
                # list_s_new = list(set_s & set(list_s_new_all)) # 取交集之后是邻居中还能用的那些
                list_s_new = list(set_s & set(list_s_new_all)) # 取交集之后是邻居中还能用的那些
                list_s_pre = list_s_new[:] #这里面其实包括了 下次要选的strain_s，但是下次要选的strain_s肯定不在set_s中
                #所以不影响
                flag_pre = True
    # print(cnt_all,len(list_dis_all),cnt_now)
    #创建结果路径
    path_strain_ani = os.path.join(path_strain,sp+'.csv')
    path_strain_tree = os.path.join(path_tree,sp+'.csv')
    # print(path_strain_ani)
    df_list_dis_all = pd.DataFrame(list_dis_all,columns=["ani"])
    df_list_dis_all.to_csv(path_strain_ani,header=None,index=None)
    df_ani.to_csv(path_strain_tree,header=None,index=None,sep="\t")
    if len(df_ani)!=cnt_all :
        print(len(df_ani),cnt_all)
        print("error: ani not match")
    
    

#把增加变异写成函数
def get_strain_genome(seq_id_init,record,record_seq_list,ani,name_j,strain_s,path_strian):
    len_genome = len(record_seq_list)
    len_base = round(len(record_seq_list)*(1-ani))
    # print("len_base",len(record_seq_list),ani,len_base)
    #现在只考虑替换和插入删除
    len_1 = random.randint(round(len_base*0.2),round(len_base*0.8)) #snp的碱基数，不能太小，也不能太大
    len_2 = len_base - len_1 #indel的碱基数，其实indel的碱基数和-indel_count并不一一对应，
    #好像是符合某个幂律分布eg P(x) = 0.5* x^(-2)
    #但是由于目前没办法确定具体len_2，就先这样看着吧
    base_list = ["A","T","C","G"]
    index_variant = random.sample(range(len(record_seq_list)),len_base)
    insert_base=[] #存放插入的碱基
    insert_index=[] #存放插入的位置
    for i in range(len(index_variant)):
        if i <len_1:
            #替换
            base_new = random.choice([base for base in base_list if base!=record_seq_list[index_variant[i]]])
            record_seq_list[index_variant[i]] = base_new # index_variant 存放的是 需要改变的下标
        else:
            flag_indel = random.choice([0,1])
            if flag_indel==0: #0的时候删除
                record_seq_list[index_variant[i]]="N" #先赋值N，到时候统一删除
                len_genome-=1
            else: #1的时候插入
                insert_index.append(index_variant[i])
                len_genome+=1
    insert_index.sort(reverse=True)
    # print(len_1,len_2,len(insert_index))
    # print(len(record_seq_list)+2*len(insert_index)-len_2)
    for i in insert_index:
        base = random.choice(base_list)
        # print(base)
        record_seq_list.insert(i,base)
    record_seq_list_new = [base for base in record_seq_list if base!="N"]
    separator = ''
    record_seq_list_new_str = separator.join(record_seq_list_new)
    record_seq_list_seq = Seq(record_seq_list_new_str)
    record.seq = record_seq_list_seq
    # path_strian  = "/data/huixingqi/software/1_tools/insilicoseq/sim_meta/strain/"
    output_file_name = os.path.join(path_strian,name_j)
    # id_record = record.id
    record.id = seq_id_init+"_"+str(strain_s)
    print("----",strain_s,record.id,output_file_name)
    # with open(output_file_name, "w") as output_file:
    #     SeqIO.write(record, output_file, "fasta")
    circulator(record,output_file_name)
    # print("len(record_seq_list_new),len_genome",len(record_seq_list_new),len_genome)
    return record_seq_list_new,len_genome

#根据输入的最小生成树文件生成基因组
def generate_strain(path_out,df_tree,genome_path,base_name):
    #输入有三个，一个是输出基因组的位置，一个是最小生成树，一个是参考基因组
    # print(1)
    with open(genome_path) as f:
        for ref_seq in SeqIO.parse(f,"fasta"):
            record_seq = ref_seq.seq
    print("ref_seq.id",ref_seq.id)
    len_strain = len(df_tree) #除了ref之外有这么多个strain
    record_seq_list = list(record_seq)
    len_base = len(record_seq_list)
    # base_name = os.path.basename(path_out)
    
    list_s1 = df_tree["s1"].to_list()
    list_s2 = df_tree["s2"].to_list()
    set_s = set(list_s1) | set(list_s2)
    list_s = list(set_s)
    list_s_new = list_s[:]
    
    cnt_strain = 0 #记录strain的个数
    
    flag_pre = False # 一开始没有邻居
    len_df_tree = len(df_tree)
    dict_genome = dict() # 保存基因组的
    seq_id_init = ref_seq.id[:]
    while(cnt_strain<len_df_tree):
        strain_s = np.random.choice(list_s_new)
        if strain_s not in dict_genome.keys():
            print("init seq")
            dict_genome[strain_s] = record_seq_list[:]
            name_out = base_name+"_"+str(strain_s)+".fasta"
            output_file_name = os.path.join(path_out,name_out)
            
            # id_seq = ref_seq.id
            
            ref_seq.id = seq_id_init+"_"+str(strain_s)
            print("----",strain_s,ref_seq.id,output_file_name)
            # with open(output_file_name, "w") as output_file:
            #     SeqIO.write(ref_seq, output_file, "fasta")
            circulator(ref_seq,output_file_name)
        
        else:
            record_seq_list = dict_genome[strain_s]
        
                
        #strain_s是我们随机选择的节点，从这个节点开始广度优先遍历最小生成树
        list_left = df_tree[df_tree["s1"]==strain_s]["s2"].to_list()
        list_right = df_tree[df_tree["s2"]==strain_s]["s1"].to_list()

        list_left_dis = df_tree[df_tree["s1"]==strain_s]["ani"].to_list()
        list_right_dis = df_tree[df_tree["s2"]==strain_s]["ani"].to_list()
        list_all = [] #把左相连和右相连的节点放到一起
        list_dis = [] #距离

        if list_left!=None:
            list_all+=list_left
            list_dis+=list_left_dis
        if list_right!=None:
            list_all+=list_right
            list_dis+=list_right_dis
        if len(list_all)!=0:
            # print("list_dis",list_dis)
            random_index = [index_list_all for index_list_all in range(len(list_all))]
            random.shuffle(random_index) #随机打乱节点顺序

            list_all_new = [list_all[ri] for ri in random_index]
            list_dis_new = [list_dis[ri] for ri in random_index] # 这个是当前节点的全部邻居

        list_s_nos = list_all_new[:]
        for lan,ldn in zip(list_all_new,list_dis_new):
            cnt_strain+=1
            print("strain_s,lan,ldn,cnt_strain",strain_s,lan,ldn,cnt_strain)
            name_j = base_name+"_"+str(lan)+".fasta"
            print(base_name,name_j,path_out)
            my_seq,my_len = get_strain_genome(seq_id_init,ref_seq,record_seq_list[:],ldn,name_j,lan,path_out)
            print(lan,my_len)
            dict_genome[lan] = my_seq
        df_tree_l = df_tree[df_tree["s1"]!=strain_s]
        df_tree = df_tree_l[df_tree_l["s2"]!=strain_s]
        list_s1 = df_tree["s1"].to_list()
        list_s2 = df_tree["s2"].to_list()
        set_s = set(list_s1) | set(list_s2) #并集，新df的所有节点
        if flag_pre:
            list_s_new_all = list_s_nos[:]+list_s_pre[:] #我的邻居加上一个的邻居
        else:
            list_s_new_all = list_s_nos[:] #第一次的时候只有我的邻居
        # list_s_new = list(set_s & set(list_s_new_all)) # 取交集之后是邻居中还能用的那些
        list_s_new = list(set_s & set(list_s_new_all)) # 取交集之后是邻居中还能用的那些
        list_s_pre = list_s_new[:] #这里面其实包括了 下次要选的strain_s，但是下次要选的strain_s肯定不在set_s中
        print("list_s_new,list_s_pre",list_s_new,list_s_pre)
        #所以不影响
        flag_pre = True
    
        
    
    
#遍历一个节点的所有相邻节点
def search_a_node(df_tree,strain_s,acns_name,len_genome_sp,record,record_seq_list,path_strain,cnt_now,cnt_all):
#先找到和节点strain_s相连的节点，有可能是(0,3) 也有可能是(3,11)
    # print("len_genome_sp",len_genome_sp)
    list_left = df_tree[df_tree["s1"]==strain_s]["s2"].to_list()
    list_right = df_tree[df_tree["s2"]==strain_s]["s1"].to_list()

    list_left_dis = df_tree[df_tree["s1"]==strain_s]["dis"].to_list()
    list_right_dis = df_tree[df_tree["s2"]==strain_s]["dis"].to_list()

    list_all = [] #把左相连和右相连的节点放到一起
    list_dis = [] #距离
    
    dis_all = 0
    
    if list_left!=None:
        list_all+=list_left
        list_dis+=list_left_dis
    if list_right!=None:
        list_all+=list_right
        list_dis+=list_right_dis
    
    if len(list_all)!=0:
        # print("list_dis",list_dis)
        random_index = [index_list_all for index_list_all in range(len(list_all))]
        random.shuffle(random_index) #随机打乱节点顺序

        list_all_new = [list_all[ri] for ri in random_index]
        list_dis_new = [list_dis[ri] for ri in random_index]
        my_seq_list = []
        my_len_list = []
        for ld in list_dis_new:
            name_j = acns_name.replace("_genomic.fna","_"+str(cnt_now)+"_genomic.fna")
            len_base = round(len_genome_sp*(1-ld))
            my_seq,my_len = get_strain_genome(record,record_seq_list,len_base,name_j,path_strain)
            my_seq_list.append((my_seq))
            my_len_list.append(my_len)
            cnt_now+=1
            if cnt_now>=cnt_all:
                return None,None,cnt_now,my_seq_list,my_len_list,True #满了之后就不用在继续执行了
        df_tree_l = df_tree[df_tree["s1"]!=strain_s]
        df_tree_new = df_tree_l[df_tree_l["s2"]!=strain_s]
    else:
        list_all_new =[]
        df_tree_new = df_tree
        my_seq_list =[]
        my_len_list = []
        # print("看看返回值",list_all_new,list_all)
    return list_all_new,df_tree_new,cnt_now,my_seq_list,my_len_list,False

#根据输入的最小生成树文件生成ani
def generate_strain_ani(sim,sp,tree_path_list,cnt_all):
    #sim是第几次模拟，用于保存模拟结果
    #tree_path_list是当前输入的tree的地址列表，用于遍历并生成某个物种的所有ani
    #cnt_all是总共需要的strain个数
    path_genome_root = "/data/huixingqi/sim_meta/data/ncbi_genome/ncbi_new_all"
    acns_path = "/data/huixingqi/sim_meta/data/taxonomy_rank/name.txt"
    acns_path_df = pd.read_csv(acns_path,header=None,sep="\t")
    acns_path_df.columns=["acns","sp"]
    acns_name = acns_path_df[acns_path_df["sp"]==sp]["acns"].unique()[0][3:]+"_genomic.fna"
    genome_path = os.path.join(path_genome_root,acns_name) # 参考基因组的地址

    path_len_genome = "/data/huixingqi/software/1_tools/insilicoseq/sim_meta/len_genome_new.txt"
    #len_genome.txt只有1224行，应该是一个样本的，我们把他扩展成所有物种的
    # for acns in /data/huixingqi/sim_meta/data/ncbi_genome/ncbi_new_all/*.fna; do
    # > len=$(seqkit stats $acns |awk 'NR==2{gsub(",","");print$5}')
    # > echo -e "$acns\t$len" >> /data/huixingqi/software/1_tools/insilicoseq/sim_meta/len_genome_new.txt
    # > done
    
    len_genome = pd.read_csv(path_len_genome,header=None,sep="\t")
    len_genome.columns=["name","len"]
    len_genome_sp_init = len_genome[len_genome["name"]==acns_name]["len"].unique()[0]
    # len_genome_sp_init = 100


    #创建结果路径
    path_strain = os.path.join("/data/huixingqi/sim_meta/strain_genome/",sim,sp)
    #一个环境中物种的strain 新生成基因组的保存地址
    if os.path.exists(path_strain)!=True:
        os.makedirs(path_strain,)

    #打开参考基因组
    # print(genome_path)
    with open(genome_path) as f:
        for record in SeqIO.parse(f,"fasta"):
            record_seq = record.seq 
            #因为更新后的基因组只有一个序列，所以record_seq只会被赋值一次
    cnt_now = 0
    name_j = acns_name.replace("_genomic.fna","_"+str(cnt_now)+"_genomic.fna")
    output_file_name = output_file_name = os.path.join(path_strain,name_j)
    with open(output_file_name, "w") as output_file:
                SeqIO.write(record, output_file, "fasta")
            #因为更新后的基因组只有一个序列，所以record_seq只会被赋值一次
    

    flag_ref = True
    # print("flag_ref",flag_ref)
    #第一次的时候，参考的基因组是ref，后面参考的基因组就是队列里面新的生成的基因组了
        

    while (cnt_now<cnt_all):
        cnt_now = 1
        # print("---cnt_now",cnt_now,"cnt_all---",cnt_all,"||||")
        tree_path = np.random.choice(tree_path_list)
        tree_path_list = list(set(tree_path_list) - set([tree_path])) 
        # print(tree_path,"||||")
        #df_tree存放了一个样本的树
        df_tree = pd.read_csv(tree_path,header=None,sep="\t")
        df_tree.columns = ["s1","s2","dis"]#存放的是strain1和strain2之间的距离
        list_s1 = df_tree["s1"].to_list()
        list_s2 = df_tree["s2"].to_list()
        set_s = set(list_s1) | set(list_s2)
        list_s = list(set_s)
        #变成set是为了取并集，变成list是为了后续处理方便
        strain_s = np.random.choice(list_s)
        #strain_s是我们随机选择的节点，从这个节点开始广度优先遍历最小生成树

        record_seq_list_init = [list(record_seq)]
        # record_seq_list_init = [list((record_seq)[:100])]
        
        list_all = [strain_s]
        my_seq_list_all = record_seq_list_init[:]
        # my_len_all =[len_genome_sp_init]
        my_len_all =[len(my_seq_list_all)]
        
        # my_seq_list_all = ["s_ref"]
        df_tree_new = df_tree.copy()
        # print(list_all)
        # cnt_break = 0
        while(len(list_all)!=0):
            # cnt_break +=1
            list_all_new_next = [] #保存下一级的节点
            my_seq_list_next = []
            my_len_next = []
            for strain_s_now,record_seq_list, len_genome_sp in zip(list_all,my_seq_list_all,my_len_all):
                # cnt_break +=1
                print("----strain_s_now-------",strain_s_now,list_all,len_genome_sp,len(record_seq_list))
                list_all_new,df_tree_new,cnt_now,my_seq_list,my_len,flag_break = search_a_node(df_tree_new,strain_s_now,acns_name,len_genome_sp,record,record_seq_list,path_strain,cnt_now,cnt_all)
                if flag_break:
                    break
                else:
                    # print("没有break")
                    if len(list_all_new)!=0:
                        list_all_new_next.append(list(list_all_new))
                        my_seq_list_next.append(list(my_seq_list))
                        my_len_next.append(list(my_len))
                        # print(list_all_new)
                        # print(list_all_new_next)
                        # print("cnt_now",cnt_now)
            # print(list_all_new_next)
            list_all = [la for lan in list_all_new_next for la in lan]
            my_seq_list_all = [la for lan in my_seq_list_next for la in lan]
            # print("为啥",my_len_next)
            my_len_all = [la for lan in my_len_next for la in lan]
            # my_len_all = my_len_next
            # my_seq_list_all = my_seq_list_next
            # print(list_all,my_seq_list_all,my_len_all)
        # print("cnt_now",cnt_now,"\n")
            # print("最后",list_all)
            # if cnt_break>50:
            #     break