#"bash /data/huixingqi/sim_meta/all/7_2_sketch.sh "+env
env=$1
path_home=$2
path_data=$path_home/data
path_sketch=$path_data/env/data/7_ani/$env/sketch

while read -r strain ; do
    # echo $strain
    base_name=$( basename $strain)
    # echo $base_name
    acns=${base_name:0:15}
    # echo $acns
    # echo "mash sketch $strain > $path_sketch/${acns}.fna.msh "
    mash sketch $strain -o $path_sketch/${acns}.fna.msh
done <$path_data/env/data/7_ani/$env/strain.txt

path_ani=$path_data/env/data/7_ani/$env/ani.csv
> $path_ani
while read -r acns_1 acns_2 ref query ; do
    base_ref_sketch=$path_sketch/${acns_1:0:15}.fna.msh
    base_query_sketch=$path_sketch/${acns_2:0:15}.fna.msh
    mash dist $base_ref_sketch $base_query_sketch >> $path_ani
done < $path_data/env/data/7_ani/$env/rq.csv