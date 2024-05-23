#step1
for i in `ls ../00.cram|grep cram$|cut -f 1 -d .`
do
	echo "perl step1.pl ../duplication.list ../01.raw.depth/output/${i}.regions.bed.gz ../01.raw.depth/output/${i}.mosdepth.summary.txt > ./output/${i}.depth" > ${i}.step1.sh
done
##step2
for i in `ls ../00.cram|grep cram$|cut -f 1 -d .`
do
        echo "perl step2.pl ../02.step1/output/${i}.depth > ./output/${i}.accession" > ${i}.step2.sh
done
##step3
for i in `ls ../00.cram|grep cram$|cut -f 1 -d .`
do
        echo "perl step3.pl ../03.step2/output/${i}.accession /vol3/agis/chengshifeng_group/shiyan/03-pea_WGS/ZW6_ref/ZW6_chromosome_genes.gff3 > ./output/${i}.cnv.vcf" > ${i}.step3.sh
done

