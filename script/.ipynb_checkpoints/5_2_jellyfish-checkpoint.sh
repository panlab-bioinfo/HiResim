file=$1
env=$2
k=$3
batch=$4
sp=$5
path_home=$6
path_data=$path_home/data
# path_kmer=$5
path_unikmer=$path_data/env/data/5_strain_unikmer/$env



if [ -d $path_unikmer/sp/$batch ]; then
    rm -rf $path_unikmer/sp/$batch
fi

mkdir $path_unikmer/sp/$batch

if [ ! -d $path_unikmer/unikmer/$batch ]; then
    mkdir $path_unikmer/unikmer/$batch
fi

if [ ! -d $path_unikmer/unikmer/$batch/$sp ]; then
    mkdir $path_unikmer/unikmer/$batch/$sp
fi
# mkdir $path_unikmer/unikmer/$batch/$sp



cd $path_unikmer/sp/$batch #这个文件夹用来暂存每个物种的strian的kmer

while read line ; do
    name=$(basename "$line" _genomic.fna.gz)
    name_sub=${name:0:15}
        jellyfish count -m $k -s 100M -t 128 -C $line -o kmer.jf
        # /data/software/bin/jellyfish count -m $k -s 100M -t 128 -C <(zcat $line) -o ${name_sub}.jf
        jellyfish dump -t -c kmer.jf > ${name_sub}.fa
        # cp ${name_sub}.fa $path_kmer
done <$file  
#读取的文件是 5_2_jellyfish.sh /data/huixingqi/data/env/data/4_sp_strain_path_sample/test/s__Anaerostipes_hadrus.txt test