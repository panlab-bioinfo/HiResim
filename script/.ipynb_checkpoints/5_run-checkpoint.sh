size=$1
env=$2
#HiResim
path_home=$3
kmer=$4
path_data=$path_home/data
path_unikmer=$path_data/env/data/5_strain_unikmer/$env

if [ ! -d $path_unikmer ]; then
    mkdir $path_unikmer
fi

path_jf=$path_unikmer/sp
if [ ! -d $path_jf ]; then
    mkdir $path_jf
fi
path_jf_unikmer=$path_unikmer/unikmer
if [ ! -d $path_jf_unikmer ]; then
    mkdir $path_jf_unikmer
fi

path_jf_unikmer_sorted=$path_unikmer/unikmer_sorted_nonull
if [ ! -d $path_jf_unikmer_sorted ]; then
    mkdir $path_jf_unikmer_sorted
fi

# path_jf_kmer_sp=$path_unikmer/kmer
# if [ ! -d $path_jf_kmer_sp ]; then
#     mkdir $path_jf_kmer_sp
# fi

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
for((i=0;i<=$size;i+=80)); do
    
    my_array+=($i)
    # ((i++))
done

len_array_all=${#my_array[@]}
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
    arry_cmd+=("python 5_get_strain_kmer.py $env $start $end $kmer $path_home")
    echo "python 5_get_strain_kmer.py $env $start $end $kmer $path_home" >"log/5_get_strain_kmer_${env}_${start}_${end}.log"
    array_log+=("log/5_get_strain_kmer_${env}_${start}_${end}.log")
    # nohup 
done

echo "${arry_cmd[@]}"

#下载基因组
# echo "python 5_0_down_strain.py $path_home $env"
# python  5_0_down_strain.py $path_home $env

# run_background_commands "${arry_cmd[@]}"
# trap 'if [ "$normal_exit" = false ]; then kill_background_commands; fi' EXIT

