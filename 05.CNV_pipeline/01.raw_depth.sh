for i in `ls ../00.cram|grep cram$|cut -f 1 -d .`
do
echo -ne "/public/agis/chengshifeng_group/xianwenfei/software/anaconda2/envs/gcc7/bin/mosdepth --no-per-base -f /vol3/agis/chengshifeng_group/shiyan/03-pea_WGS/ZW6_ref/01-index/ZW6_ref.fa  -t 4 -b ../unique.gene.bed -Q 0 ./output/$i ../00.cram/${i}.deduped.cram " >$i.depth.sh
done
