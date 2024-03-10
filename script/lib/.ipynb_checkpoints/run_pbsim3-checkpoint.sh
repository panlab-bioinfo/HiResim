path_out=$1
method=$2
model=$3
depth=$4
genome=$5
prefix=$6
pass_num=$7
sample_n=$8
flag_css=$9

sim_res=sim_res_${sample_n}
sim_concat=sim_concat_${sample_n}
# sim_strain=sim_strain${sample_n}
# if [ ! -d $path_out/$sim_res ];then
#     echo "mkdir $path_out/$sim_res"
#     mkdir -p $path_out/$sim_res
# fi

# if [ ! -d $path_out/$sim_strain ];then
#     echo "mkdir $path_out/$sim_strain"
#     mkdir -p $path_out/$sim_strain
# fi

# if [ ! -d $path_out/$sim_concat ];then
#     echo "mkdir $path_out/$sim_concat"
#     mkdir -p $path_out/$sim_concat
# fi

# if [ ! -f $path_out/$sim_concat/sim.fastq ];then
#     echo "mkdir $path_out/$sim_concat/sim.fastq "
#     > $path_out/$sim_concat/sim.fastq
# fi

# echo $genome
genome_name=$(basename $genome)
# echo $genome_name
# genome_new=$path_out/$sim_strain/$genome_name
genome_new=$genome
# echo "gzip -dc ${genome}.gz > $genome_new"
# gzip -dc ${genome}.gz > $genome_new


path_hifi=$( realpath pbsim_hifi.sh)
path_ont=$( realpath pbsim_ont.sh)

cd $path_out/$sim_res
# source activate env_mcss
# echo "/data/huixingqi/conda/envs/env_mcss/bin/pbsim --strategy wgs --method $method --$method $model --depth $depth --genome $genome_new --prefix $prefix"

# /data/huixingqi/conda/envs/env_mcss/bin/pbsim --strategy wgs --method $method --$method $model --depth $depth --genome $genome_new --prefix $prefix --pass-num $pass_num

# /data/huixingqi/conda/envs/env_mcss/bin/samtools view -bS ${prefix}_0001.sam -o ${prefix}_0001.subreads.bam
# /data/huixingqi/conda/envs/env_mcss/bin/ccs ${prefix}_0001.subreads.bam ${prefix}_0001.ccs.fastq.gz
# gunzip ${prefix}_0001.ccs.fastq.gz
# cat ${prefix}_0001.ccs.fastq >> $path_out/$sim_concat/sim.fastq

# > $path_out/$sim_concat/pbsim.log
# flag_css="hifi"
if [ $flag_css == "hifi" ]; then
    pbsim --strategy wgs --method $method --$method $model --depth $depth --genome $genome_new --prefix $prefix --pass-num $pass_num
    samtools view -bS ${prefix}_0001.sam -o ${prefix}_0001.subreads.bam
    ccs ${prefix}_0001.subreads.bam ${prefix}_0001.ccs.fastq.gz
    gunzip ${prefix}_0001.ccs.fastq.gz
    rm ${prefix}_0001.maf
    rm ${prefix}_0001.ref
    rm ${prefix}_0001.sam
    rm ${prefix}_0001.subreads.bam
    rm ${prefix}_0001.ccs.ccs_report.txt
    rm ${prefix}_0001.ccs.zmw_metrics.json.gz
    
    # python change_ccs.py ${prefix}_0001.ccs.fastq
    # cat ${prefix}_0001.ccs.fastq >> $path_out/$sim_concat/sim.fastq
    # bash $path_hifi $prefix $path_out $method $model $depth $genome_new $pass_num $sim_res $sim_concat & >> $path_out/$sim_concat/pbsim.log
else
    pbsim --strategy wgs --method $method --$method $model --depth $depth --genome $genome_new --prefix $prefix
    # cat ${prefix}_0001.fastq >> $path_out/$sim_concat/sim.fastq
    # bash $path_ont $prefix $path_out $method $model $depth $genome_new $pass_num $sim_res $sim_concat & >> $path_out/$sim_concat/pbsim.log
fi

# awk '{ if (NR % 4 == 2) { gsub(/[^ATCG]/, substr("ATCG", int(rand()*4) + 1, 1)); } } 1' sim.fastq > sim_old.fastq
# awk 'BEGIN {i=1} (NR%4 == 1) {print "@sim." i" "i"/"1; i++} NR%4 != 1 {print $0}' sim_old.fastq > sim.fastq
# rm sim_old.fastq