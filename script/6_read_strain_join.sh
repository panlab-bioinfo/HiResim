#把read的kmer和strain的unikemr相交
env=$1
srr_start=$2
srr_end=$3
path_home=$4
path_data=$path_home/data
echo "$srr_start $srr_end"
echo "$path_data/env/data/6_read_strain_inter/read_batch_path/${env}_path.txt"
echo $(sed -n "$srr_start,${srr_end}p" $path_data/env/data/6_read_strain_inter/read_batch_path/${env}_path.txt)

path_unikmer=$path_data/env/data/5_strain_unikmer/$env/unikmer_sorted_nonull 


path_df_stats_root=$path_data/env/data/6_read_strain_inter/$env



if [ ! -d $path_df_stats_root ]; then
    mkdir $path_df_stats_root
fi

while read -r srr; do
    basename_srr=$( basename $srr)
    if [ ! -d $path_df_stats_root/$basename_srr ]; then
        mkdir -p $path_df_stats_root/$basename_srr
    fi
    #创建保存交集的路径,每个样本一个文件夹，因为read是按照样本生成的
    for sp in $srr/*.fa; do
        basename_sp=$( basename $sp .fa)
        echo $path_unikmer/${basename_sp}
        if [ -d $path_unikmer/${basename_sp} ]; then
            # echo "yes"
            #样本中每个物种都有自己的文件夹，保存物种下strian的unikmer和readkmer的交集
            if [ ! -d $path_df_stats_root/$basename_srr/$basename_sp ]; then
                mkdir $path_df_stats_root/$basename_srr/$basename_sp
            fi
        
            for fa in $path_unikmer/${basename_sp}/*; do
                basename_fa=$( basename $fa)
                join "$sp" "$fa" > $path_df_stats_root/$basename_srr/$basename_sp/${basename_fa}
                if [ ! -s $path_df_stats_root/$basename_srr/$basename_sp/${basename_fa} ]; then
                    rm $path_df_stats_root/$basename_srr/$basename_sp/${basename_fa}
                fi
            done
        fi
    done

done < <(sed -n "$srr_start,${srr_end}p" $path_data/env/data/6_read_strain_inter/read_batch_path/${env}_path.txt)

