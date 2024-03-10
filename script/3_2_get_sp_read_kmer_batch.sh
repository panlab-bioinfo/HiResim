#3_get_sp_read.py已经得到了物种的read文件，保存在read_python中
#现在把read_python中的测序数据，跑jellyfish得到read的kmer，不是uni
#3_2_get_sp_read_kmer.sh 的批量版本
# bash 3_2_get_sp_read_kmer_batch.sh marine 2 /data/huixingqi/data/env/data/3_read_python/batch/marine/0_10.txt
env=$1
mood=$2
batch=$3
k=$4
#HiResim/data
path_data=$5

path_kmer=$path_data/env/data/3_sp_read_kmer/$env$k
path_kmer_sort=$path_data/env/data/3_sp_read_kmer_sorted/$env$k
path_read=$path_data/env/data/3_read_python/$env
# if [ ! -d $path_kmer ]; then
#     mkdir $path_kmer
# fi

# if [ ! -d $path_kmer_sort ]; then
#     mkdir $path_kmer_sort
# fi

cd $path_data/env/data/3_sp_read_kmer/tmp_fastq
#保存拼接的_1和_2

if [ $mood -eq 1 ]; then
    while read -r srr_sp; do
    # for srr in $path_read/SRR*; do
        base_srr=$(basename $srr_sp _cutsp.txt)
        if [ ! -d $path_kmer/$base_srr ]; then
            mkdir $path_kmer/$base_srr
        fi
        if [ ! -d $path_kmer_sort/$base_srr ]; then
            mkdir $path_kmer_sort/$base_srr
        fi
        srr=$path_read/$base_srr
        for sp in ${srr}/*.fastq; do
            base_sp=$(basename $sp .fastq)
            jellyfish count -m $k -s 100M -t 128 -C $sp -o ${base_srr}_${base_sp}.jf
            jellyfish dump -t -c ${base_srr}_${base_sp}.jf > $path_kmer/$base_srr/${base_sp}.fa
            sort -k1 $path_kmer/$base_srr/${base_sp}.fa  > $path_kmer_sort/$base_srr/${base_sp}.fa
        done
    done < $batch
else
    while read -r srr_sp; do
    # for srr in $path_read/SRR*; do
        # echo $srr
        base_srr=$(basename $srr_sp _cutsp.txt)
        if [ ! -d $path_kmer/$base_srr ]; then
            mkdir $path_kmer/$base_srr
        fi
        if [ ! -d $path_kmer_sort/$base_srr ]; then
            mkdir $path_kmer_sort/$base_srr
        fi
        srr=$path_read/$base_srr
        echo $srr
        for sp in ${srr}/*_1.fastq; do
            echo $sp
            base_sp=$(basename $sp _1.fastq)
            sp2=${sp/"_1.fastq"/"_2.fastq"}
            # echo ${sp2}
            cat $sp $sp2 > ${base_srr}_${base_sp}.fastq
            #之前是>cat.fastq，但是batch后会冲突，所以名字要不一样了
            #不能单纯的用物种名称，因为不同样本物种可能会重
            jellyfish count -m $k -s 100M -t 128 -C ${base_srr}_${base_sp}.fastq -o ${base_srr}_${base_sp}.jf
            
            jellyfish dump -t -c ${base_srr}_${base_sp}.jf > $path_kmer/$base_srr/${base_sp}.fa
            
            sort -k1 $path_kmer/$base_srr/${base_sp}.fa  > $path_kmer_sort/$base_srr/${base_sp}.fa
        done
    done< $batch
fi



#因为read_strain_inter在3中只生成路径，交集是6生成的，所以前缀还是6,这个已经在6中生成了