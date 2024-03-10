#根据得到的ani生成最小生成树
# python 7_2_prime.py oral
import sys
import numpy as np
import os
import pandas as pd
import shutil

env = sys.argv[1]
data_home=sys.argv[2]

#先定义最小生成树的方法
class Graph:
    def __init__(self, vertices):
        self.V = vertices
        self.graph = [[0 for _ in range(vertices)] for _ in range(vertices)]

    def min_key(self, key, mst_set):
        min_val = sys.maxsize
        min_index = None
        for v in range(self.V):
            if key[v] < min_val and not mst_set[v]:
                min_val = key[v]
                min_index = v
        return min_index

    def prim(self):
        key = [sys.maxsize] * self.V
        parent = [None] * self.V
        key[0] = 0
        mst_set = [False] * self.V
        parent[0] = -1

        for _ in range(self.V):
            u = self.min_key(key, mst_set)
            mst_set[u] = True
            for v in range(self.V):
                if self.graph[u][v] >= 0 and not mst_set[v] and key[v] > self.graph[u][v]:
                    key[v] = self.graph[u][v]
                    parent[v] = u

        return self.get_mst_edges(parent)

    def get_mst_edges(self, parent):
        mst = np.zeros((self.V, self.V), dtype=float)
        for i in range(1, self.V):
            mst[parent[i]][i] = self.graph[i][parent[i]]
            mst[i][parent[i]] = self.graph[i][parent[i]]
        return mst


path_ani = os.path.join(data_home,"data/env/data/7_ani",env)
path_matrix_root = os.path.join(path_ani,"matrix")

# if not os.path.exists(path_matrix_root):
#     os.makedirs(path_matrix_root)

path_tree_root = os.path.join(path_ani,"tree")
if not os.path.exists(path_tree_root):
    os.makedirs(path_tree_root)

list_srr = sorted(os.listdir(path_matrix_root))
for srr in list_srr:
    path_srr = os.path.join(path_matrix_root,srr)
    path_srr_tree = os.path.join(path_tree_root,srr)
    if os.path.exists(path_srr_tree):
        shutil.rmtree(path_srr_tree)
    os.makedirs(path_srr_tree)
    for sp in sorted(os.listdir(path_srr)):
        path_sp =  os.path.join(path_srr,sp)
        path_sp_tree =  os.path.join(path_srr_tree,sp)
        matrix_sp = np.loadtxt(path_sp)
        n = len(matrix_sp)
        g = Graph(n)
        g.graph = matrix_sp
        mst = g.prim()
        list_tree = [] #用来记录节点和距离
        # cnt=0
        for i in range(n-1):
            for j in range(i+1,n):
                dis_ij = mst[i][j]
                if dis_ij!=0:
                    list_tree.append([i,j,1-dis_ij])
                    # cnt+=1
        df_tree = pd.DataFrame(list_tree)
        # print(cnt,len(df_tree),path_sp_tree)
        df_tree.to_csv(path_sp_tree,sep="\t",index=None,header=None)
