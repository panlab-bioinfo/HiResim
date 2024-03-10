size=$1
env=$2
mood=$3
#HiResim/data
path_home=$4
echo $path_home
path_data=$path_home/data
read_path=$5
path_db=$6
# HiResim/data
path_res=$path_data/env/res
path_out=$path_res/out
path_batch=$path_res/batch
path_batch_env=$path_batch/$env
path_out_env=$path_out/$env
out_file=$path_data/env/data/1_out_new/$env
path_abu=$path_data/env/data/1_abundance/$env


if [ ! -d $path_abu ] ; then
    mkdir $path_abu
fi


if [ ! -d $path_res ]; then
    mkdir $path_res
fi

if [ ! -d $path_out ]; then
    mkdir $path_out
fi

if [ ! -d $path_batch ]; then
    mkdir $path_batch
fi

if [ ! -d $path_batch_env ]; then
    mkdir $path_batch_env
fi

if [ ! -d $path_out_env ]; then
    mkdir $path_out_env
fi

if [ ! -d $out_file ]; then
    mkdir $out_file
fi


run_background_commands(){
    local_arr=("$@")
    bg_pids=()
    i_local=0
    for localcmd in "${local_arr[@]}"; do
        echo "$localcmd" 
        # echo "$localcmd" >${array_log[$i_local]}
        $localcmd &>>${array_log[$i_local]} &
        # $localcmd &
        
        bg_pids+=($!)
        ((i_local++))
    done
    normal_exit=false
    wait
    normal_exit=true
}

normal_exit=true
kill_background_commands() {
    for pid in "${bg_pids[@]}"; do
        kill "$pid"
    done
}
echo "size" $size
my_array=()
for((i=0;i<=$size;i+=10)); do
    echo $i
    my_array+=($i)
    # ((i++))
done

echo my_array "${my_array[@]}"

len_array_all=${#my_array[@]}
echo len_array_all $len_array_all
if [ $len_array_all -eq 1 ]; then
    # echo "111"
    my_array+=($size)
fi

if ((my_array[-1]<$size)); then
    my_array[-1]=$size
fi

echo my_array "${my_array[@]}"
len_array=$((${#my_array[@]}-1))
echo "len_array" $len_array

arry_cmd=()
array_log=()
for ((i=0;i<len_array;i+=1)); do
    # echo $i" "$((i+1))
    echo ${my_array[$i]}" "${my_array[$((i+1))]}
    start=${my_array[$i]}
    end=${my_array[$((i+1))]}
    echo "1 kraken" >/data/huixingqi/sim_meta/HiResim/1.txt
    arry_cmd+=("python 1_kraken.py $env $start $end $mood $read_path $path_db $path_home")
    echo "python 1_kraken.py $env $start $end $mood" >"log/1_kraken${env}_${start}_${end}_${mood}.log"
    array_log+=("log/1_kraken${env}_${start}_${end}_${mood}.log")
    # nohup 
done

# echo "${arry_cmd[@]}"

run_background_commands "${arry_cmd[@]}"
trap 'if [ "$normal_exit" = false ]; then kill_background_commands; fi' EXIT



file_path="$path_data/env/data/1_read_cnt/${env}_read_cnt.txt"
>$file_path
for srr in $path_out_env/*.out; do
    wc -l $srr >>$file_path
done