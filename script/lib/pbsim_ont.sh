prefix=$1
path_out=$2
method=$3
model=$4
depth=$5
genome_new=$6
pass_num=$7
sim_res=$8
sim_concat=$9

cd $path_out/$sim_res

pbsim --strategy wgs --method $method --$method $model --depth $depth --genome $genome_new --prefix $prefix
cat ${prefix}_0001.fastq >> $path_out/$sim_concat/sim.fastq