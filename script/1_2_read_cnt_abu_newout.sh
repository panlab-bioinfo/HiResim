#读取srr.out,得到srr的read数
#bash 1_2_get_read_cnt.sh test batch
#把get_newout整合到这里面来了
#把get_read_abu也整合进来吧
env=$1
batch=$2
path_home=$3
path_data=$path_home/data
path_out=$path_data/env/res/out/$env
#这个其实是输入文件，只是输入的.out的结果
# file_path="/data/huixingqi/data/env/data/1_read_cnt/${1}_read_cnt.txt"
# echo $file_path
# echo $path_out
# >$file_path
#get_newout的文件
out_file=$path_data/env/data/1_out_new/$env

# out_file=$path_out_new
# if [ ! -d $out_file ]; then
#     mkdir $out_file
# fi

# for i in $path_out/*.out; do
#     # echo $i
#     wc -l $i >>$file_path

# done

#abu
path_abu="$path_data/env/data/1_abundance/"$env
# if [ ! -d $path_abu ] ; then
#     mkdir $path_abu
# fi

while read -r srr; do
    name=$path_out/${srr}.out
    rep=$path_out/${srr}_rep.txt
    # echo $name
    # wc -l $name >>$file_path
    # echo "wc -l $name >>$file_path"
    awk -F '[\t ]' '{print $2"\t"$3}' $name>$out_file/${srr}.out
    # echo "awk -F '[\t ]' '{print $2"\t"$3}' $name>$out_file/${srr}.out"
    grep "s__" $rep |awk -F '[| \t]' '{print$(NF-2)"_"$(NF-1)"\t"$NF}'  >$path_abu/${srr}_abundance.txt 
    # echo "grep "s__" $rep |awk -F '[| \t]' '{print$7"_"$8"\t"$9}'   >$path_abu/${srr}_abundance.txt"
    
done < $batch

