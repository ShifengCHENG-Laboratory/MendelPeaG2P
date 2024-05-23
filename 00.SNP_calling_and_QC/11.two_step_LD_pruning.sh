for i in `ls ../../06.remove_2_sample/*maf_retain.remove2sample.vcf.gz`
do
    chr=`basename $i | cut -f 1 -d "."`
    echo "
    plink --vcf $i --allow-extra-chr --allow-no-sex --make-bed --double-id  --out ../$chr.merge

    plink --bfile ../$chr.merge --indep-pairwise 10kb 1 0.8 --out ../$chr.step1 --allow-extra-chr
    plink --bfile ../$chr.merge --extract ../$chr.step1.prune.in --out ../$chr.prunedata.step1 --make-bed --allow-extra-chr --threads 10

    plink --bfile ../$chr.prunedata.step1 --indep-pairwise 50 1 0.8 --out ../$chr.step2 --allow-extra-chr
    plink --bfile ../$chr.prunedata.step1 --extract ../$chr.step2.prune.in --out ../$chr.prunedata.step2 --make-bed --allow-extra-chr --threads 10
    " > $chr.thin.sh

done
