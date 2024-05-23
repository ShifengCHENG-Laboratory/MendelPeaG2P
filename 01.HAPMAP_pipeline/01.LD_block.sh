for a in chr1 chr2 chr3 chr4 chr5 chr6 chr7
do
echo "
/public/agis/chengshifeng_group/fengcong/WGRS/software/plink/plink --vcf ../00.phased_vcf/${a}.snp.raw.HARD.Missing-unphasing.ID.allele2_retain.hard_retain.InbreedingCoeff_retain.missing_retain.maf_retain.ann.phased.vcf.gz --geno 0.2 --blocks-min-maf 0.05 --blocks no-pheno-req --out ${a}.ldblock --blocks-max-kb 1000 --allow-extra-chr --threads 10 --memory 200000">${a}.sh
done
