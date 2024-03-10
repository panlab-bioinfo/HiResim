env=$1
len_sample=$2
mood=$3
abu=$4
n_sample=$5
read_path=$6
path_db=$7
kmer=21
#kraken
path_home=$(pwd)
echo $path_home
path_data=$path_home/data
path_4=$path_data/env/data/4_sp_strain_path_sample/$env
cd script

path_6_batch=$path_data/env/data/6_read_strain_inter/read_batch_path
if [ ! -d $path_6_batch ]; then
    mkdir $path_6_batch
fi 

# echo "bash 1_run.sh $len_sample $env $mood $path_home $read_path $path_db &> log/run_1_kraken_${env}.log"
# bash 1_run.sh $len_sample $env $mood $path_home $read_path $path_db &> log/run_1_kraken_${env}.log

# echo "2 get sp"
# echo "python 2_get_sp.py $env $abu $path_home &> log/run_2_get_sp_${env}.log"
# python 2_get_sp.py $env $abu $path_home &> log/run_2_get_sp_${env}.log

# python 3_change_fastq_name.py $env $mood &> log/run_3_change_fastq_name_${env}.log
# echo "3 sp read"
# echo "bash 3_run.sh $len_sample $env $mood $kmer $path_home $read_path &> log/run_3_get_sp_read_${env}.log"
# bash 3_run.sh $len_sample $env $mood $kmer $path_home $read_path &> log/run_3_get_sp_read_${env}.log
# echo "4 strain path"
# echo "python 4_get_strain_path.py $env $path_home &> log/run_4_get_strain_path_${env}.log"
# python 4_get_strain_path.py $env$kmer $path_home &> log/run_4_get_strain_path_${env}.log

# echo "5 strain kmer"
# len_sp=$(ls $path_4$kmer |wc -l)
# echo "bash 5_run.sh $len_sp $env$kmer $path_home $kmer &> log/run_5_${env}.log"
# bash 5_run.sh $len_sp $env$kmer $path_home $kmer &> log/run_5_${env}.log
# echo "6 determine strain"
# echo "bash 6_run.sh $len_sample $env$kmer $path_home &> log/run_6_${env}.log"
# bash 6_run.sh $len_sample $env$kmer $path_home &> log/run_6_${env}${kmer}.log
# echo "python 6_read_strain_kemr_new.py $env$kmer $path_home &> log/run_6_read_strain_kemr_new_${env}${kmer}.log"
# python 6_read_strain_kemr_new.py $env$kmer $path_home &> log/run_6_read_strain_kemr_new_${env}${kmer}.log
# echo "7 strain ani "
# echo "python 7_get_ani.py $env$kmer $path_home &> log/run_7_get_ani_${env}.log"
# python 7_get_ani.py $env$kmer $path_home &> log/run_7_get_ani_${env}.log
# echo "python 7_2_prime.py $env$kmer $path_home &> log/run_7_2_prime_${env}.log"
# python 7_2_prime.py $env$kmer $path_home &> log/run_7_2_prime_${env}.log
# echo "8 strain cnt and abu "
# echo "python 8_strian_cnt_and_abu2.py $env$kmer $path_home &> log/run_8_strian_cnt_and_abu_${env}.log"
# python 8_strian_cnt_and_abu2.py $env$kmer $path_home &> log/run_8_strian_cnt_and_abu_${env}.log
echo "9 sp cnt"
echo "python 9_sp_cnt.py $env$kmer $path_home &> log/run_9_sp_cnt_${env}.log"
python 9_sp_cnt.py $env$kmer $path_home &> log/run_9_sp_cnt_${env}.log
# python 9_2_generate_community_train.py --env $env --n $n_sample &> log/run_9_2_generate_community_train_${env}.log