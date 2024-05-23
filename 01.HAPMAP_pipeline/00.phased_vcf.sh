for a in chr1 chr2 chr3 chr4 chr5 chr6 chr7
do
echo "
java -Xmx200g -jar /public/agis/chengshifeng_group/fengcong/WGRS/software/beagle/beagle.21Apr21.304.jar gt=/vol3/agis/chengshifeng_group/shiyan/03-pea_WGS/003.SNP_recalling_773_V2/03-ann/01-SNP/01-MAF/${a}/${a}.snp.raw.HARD.Missing-unphasing.ID.allele2_retain.hard_retain.InbreedingCoeff_retain.missing_retain.maf_retain.ann.vcf.gz out=./${a}.snp.raw.HARD.Missing-unphasing.ID.allele2_retain.hard_retain.InbreedingCoeff_retain.missing_retain.maf_retain.ann.phased nthreads=10
/public/agis/chengshifeng_group/fengcong/WGRS/software/bcftools1.9/bin/bcftools index ${a}.snp.raw.HARD.Missing-unphasing.ID.allele2_retain.hard_retain.InbreedingCoeff_retain.missing_retain.maf_retain.ann.phased.vcf.gz
">${a}.phased_popfilter.sh
done
