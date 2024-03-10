size=$1
env=$2
mood=$3
kmer=$4
#HiResim/data
path_home=$5
read_path=$6
path_data=$path_home/data
path_batch=$path_data/env/data/3_read_python/batch/$env
path_kmer=$path_data/env/data/3_sp_read_kmer/$env$kmer
path_kmer_sort=$path_data/env/data/3_sp_read_kmer_sorted/$env$kmer
path_sp=$path_data/env/data/3_read_python/$env

echo "path_home" $path_home
echo $path_sp
if [ ! -d $path_sp ]; then
    mkdir $path_sp
    echo ok
fi
echo $path_batch
if [ ! -d $path_batch ]; then
    mkdir $path_batch
    echo ok
    
fi
echo $path_kmer
if [ ! -d $path_kmer ]; then
    mkdir $path_kmer
    echo ok
    
fi
echo $path_kmer_sort
if [ ! -d $path_kmer_sort ]; then
    mkdir $path_kmer_sort
    echo ok
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

my_array=()
for((i=0;i<=$size;i+=10)); do
    
    my_array+=($i)
    # ((i++))
done

len_array_all=${#my_array[@]}
echo len_array_all $len_array_all
if [ $len_array_all -eq 1 ]; then
    # echo "111"
    my_array+=($size)
fi

if ((my_array[-1]<$size)); then
    my_array[-1]=$size
fi
# echo "${my_array[@]}"
len_array=$((${#my_array[@]}-1))
echo $len_array

arry_cmd=()
array_log=()
for ((i=0;i<len_array;i+=1)); do
    # echo $i" "$((i+1))
    echo ${my_array[$i]}" "${my_array[$((i+1))]}
    start=${my_array[$i]}
    end=${my_array[$((i+1))]}
    arry_cmd+=("python 3_get_sp_read.py $env $start $end $mood $kmer $path_home $read_path")
    echo "python 3_get_sp_read.py $env $start $end $mood $kmer $path_home $read_path" >"log/3_get_sp_read${env}_${start}_${end}_${mood}.log"
    array_log+=("log/3_get_sp_read${env}_${start}_${end}_${mood}.log")
    # nohup 
done

# echo "${arry_cmd[@]}"

run_background_commands "${arry_cmd[@]}"
trap 'if [ "$normal_exit" = false ]; then kill_background_commands; fi' EXIT

echo "finished"




