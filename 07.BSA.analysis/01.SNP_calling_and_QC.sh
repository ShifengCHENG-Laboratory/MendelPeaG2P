for i in `cat  chr.list`
            do
            echo "bcftools mpileup -f /vol3/agis/chengshifeng_group/shiyan/03-pea_WGS/ZW6_ref/01-index/ZW6_ref.fa --redo-BAQ --min-BQ 30 --per-sample-mF --annotate FORMAT/AD,FORMAT/DP --regions ${i} -Ou --bam-list  cram.list |bcftools call -mv |bgzip -c >F250.${i}.vcf.gz
bcftools index F250.${i}.vcf.gz " >F250.${i}.call.sh
            done
for i in `ls|grep gz$|cut -f 1-2 -d .`
do
echo "bcftools filter -g 3 -G10 -e '%QUAL<10 ||(RPB<0.1 && %QUAL<15) || (AC<2 && %QUAL<15)||MQ<30||MQSB <=0.1' ${i}.vcf.gz |bgzip -c >${i}.filter.vcf.gz
bcftools index ${i}.filter.vcf.gz
done

