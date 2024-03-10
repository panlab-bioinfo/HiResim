size=$1
env=$2
path_home=$3
path_data=$path_home/data



> $path_data/env/data/6_read_strain_inter/read_batch_path/${env}_path.txt
for srr in $path_data/env/data/3_sp_read_kmer_sorted/$env/*; do
echo $srr>>$path_data/env/data/6_read_strain_inter/read_batch_path/${env}_path.txt
done



run_background_commands(){
    local_arr=("$@")
    bg_pids=()
    i_local=0
    for localcmd in "${local_arr[@]}"; do
        echo "$localcmd"
        # $localcmd &>>${array_log[$i_local]} &
        $localcmd &>>${array_log[$i_local]} &
        # $localcmd
        # ((i_local++))
        bg_pids+=($!)
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
if [ $len_array_all -eq 1 ]; then
    # echo "111"
    my_array+=($size)
fi

if ((my_array[-1]<$size)); then
    my_array[-1]=$size
fi
# echo "${my_array[@]}"
len_array=$((${#my_array[@]}-1))
echo len_array $len_array
echo my_array "${my_array[@]}"

arry_cmd=()
array_log=()
for ((i=0;i<len_array;i+=1)); do
    echo $i" "$((i+1))
    echo ${my_array[$i]}" "${my_array[$((i+1))]}
    start=$((${my_array[$i]}+1))
    end=${my_array[$((i+1))]}
    arry_cmd+=("bash 6_read_strain_join.sh $env $start $end $path_home")
    echo "bash 6_read_strain_join.sh $env $start $end $path_home" >"log/6_read_strain_join_${env}_${start}_${end}.log"
    array_log+=("log/6_read_strain_join_${env}_${start}_${end}.log")
    # nohup 
done

# echo "${arry_cmd[@]}"

run_background_commands "${arry_cmd[@]}"
trap 'if [ "$normal_exit" = false ]; then kill_background_commands; fi' EXIT

