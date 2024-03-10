#因为有些基因组有N以及多contig
#整合成一个基因组
#这个已经执行完了，得到了结果，所以不需要再跑了
import random
import gzip
import pandas as pd
import os
from Bio import SeqIO
from Bio.Seq import Seq

def acns_path_fun(acn_a,path_genome = "/data/huixingqi/sim_meta/data/ncbi_genome/gtdb_genomes_reps_r207"):
    name_acn = acn_a.split("_")
    # print(name_acn)
    path_a = os.path.join(path_genome,name_acn[1],name_acn[2][:3],name_acn[2][3:6],name_acn[2][6:9])
    genome = os.path.join(path_a,os.listdir(path_a)[0])
    return genome

#随机替换N
def replace_n_with_random(sequence):
    replaced_sequence = ""
    nucleotides = ['A', 'T', 'C', 'G']
    for base in sequence:
        if base not in nucleotides:
            replaced_sequence += random.choice(nucleotides)
        else:
            replaced_sequence += base
    return replaced_sequence

def process_fasta_file(input_file, output_file_name):
    with gzip.open(input_file, "rt") as fasta_gz_file:
        # has_n = any("N" in record.seq for record in SeqIO.parse(fasta_gz_file, "fasta"))
        # print("has_n",has_n)
        records = []
        # if has_n:
        cnt=0
            for record in SeqIO.parse(fasta_gz_file, "fasta"):
            if "N" in record.seq:
                replaced_sequence = replace_n_with_random(str(record.seq))
                new_record = record
                new_record.seq = Seq(replaced_sequence)
                records.append(new_record)
                # print("N",cnt)
            else:
                records.append(record)
                # print("no N",cnt)
            cnt+=1

        # print("records",len(records))

        with open(output_file_name, "w") as output_file:
            SeqIO.write(records, output_file, "fasta")


acns_path = "/data/huixingqi/sim_meta/data/taxonomy_rank/name.txt"
acns_path_df = pd.read_csv(acns_path,header=None,sep="\t")
acns_path_df.columns=["acns","sp"]
path_genome_list=[]
for i in acns_path_df["acns"].tolist():
    path_genome_list.append(acns_path_fun(i))
    
for i in path_genome_list:
    name_out = i.split("/")[-1].replace(".gz","")
    path_out = "/data/huixingqi/sim_meta/data/ncbi_genome/ncbi_new_genome/"+name_out
    # print(i,"\n",path_out)
    process_fasta_file(i,path_out)