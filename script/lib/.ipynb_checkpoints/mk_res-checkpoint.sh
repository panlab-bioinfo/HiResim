path_out=$1
sample_n=$2

sim_res=sim_res_${sample_n}
sim_concat=sim_concat_${sample_n}

if [ ! -d $path_out/$sim_res ];then
    echo "mkdir $path_out/$sim_res"
    mkdir -p $path_out/$sim_res
fi

if [ ! -d $path_out/$sim_concat ];then
    echo "mkdir $path_out/$sim_concat"
    mkdir -p $path_out/$sim_concat
fi