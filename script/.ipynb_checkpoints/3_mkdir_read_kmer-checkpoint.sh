#每次跑3的时候总会跑不全，怀疑是os.makedirs的问题，现在吧创建文件夹这一步用bash实现
my_dir=$1
if [ ! -d $my_dir ]; then
    mkdir $my_dir
fi