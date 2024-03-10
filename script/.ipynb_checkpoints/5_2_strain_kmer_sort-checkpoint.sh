sp=$1
env=$2
batch=$3
path_home=$4
path_data=$path_home/data
path_unikmer=$path_data/env/data/5_strain_unikmer/$env
path_sp=$path_unikmer/unikmer/$batch
path_sort=$path_unikmer/unikmer_sorted_nonull


# if [ ! -d $path_sort ]; then
#     mkdir $path_sort
# fi

if [ ! -d $path_sort/$sp ]; then
    mkdir $path_sort/$sp
fi

for fa in $path_sp/$sp/GC*.fa; do
    if [ -s "$fa" ]; then
        basename_fa=$( basename $fa)
        sort -k1 $fa > $path_sort/$sp/$basename_fa
    fi
done  
# for sp in $path_sp/*; do
#     basename_sp=$( basename $sp )
    
#     if [ ! -d $path_sort/$basename_sp ]; then
#         mkdir $path_sort/$basename_sp
#     fi
#     for fa in $sp/GC*.fa; do
#         if [ -s "$fa" ]; then
#             basename_fa=$( basename $fa)
#             sort -k1 $fa > $path_sort/$basename_sp/$basename_fa
#         fi
#     done
# done


