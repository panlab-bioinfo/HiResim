path_out=$1
sample_n=$2
sim_res=sim_res_${sample_n}
sim_concat=sim_concat_${sample_n}

>$path_out/$sim_concat/sim.fastq

cd $path_out/$sim_res

for fq in $path_out/$sim_res/*.fastq; do
    cat $fq >> $path_out/$sim_concat/sim.fastq
done

awk '{ if (NR % 4 == 2) { gsub(/[^ATCG]/, substr("ATCG", int(rand()*4) + 1, 1)); } } 1' $path_out/$sim_concat/sim.fastq > $path_out/$sim_concat/sim_old.fastq
awk 'BEGIN {i=1} (NR%4 == 1) {print "@sim." i" "i"/"1; i++} NR%4 != 1 {print $0}' $path_out/$sim_concat/sim_old.fastq > $path_out/$sim_concat/sim.fastq
rm $path_out/$sim_concat/sim_old.fastq