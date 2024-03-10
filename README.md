# HiResim

## installation

1. clone repository
```
git clone https://github.com/panlab-bioinfo/HiResim
```

2.  enter HiResim path
```
cd HiResim
```
3. create python environment
```
# The tool is currently developed based on amd64 POSIX and has undergone testing only on the Linux operating system.
# We are dedicated to resolving any issues encountered by users on different operating systems and have plans to introduce support for additional systems in future updates.
conda env create -f env_HiResim.yaml
pip install ctgan
```
    
4. activate environment
```
conda activate env_HiResim / source activate env_HiResim
```

## usage

1. generate community 

```
python script/9_2_generate_community.py --env gut21 --n 1
   
```

2. generate strains 

```
python script/10_generate_strain.py PATH_Hiresim/data/env/data/10_sim/gut21 sim_gut
   
```

3. generate reads 

```
python script/10_pbsim.py PATH_Hiresim/data/env/data/10_sim/gut21/sp_strain_abu /data/huixingqi/sim_meta/HiResim/sim_gut 1
   
```

### citation
if you use HiResim, please cite:





