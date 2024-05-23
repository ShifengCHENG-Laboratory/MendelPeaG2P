for a in chr1 chr2 chr3 chr4 chr5 chr6 chr7
do
echo "/public/home/fengcong/anaconda2/envs/py3/bin/python /vol3/agis/chengshifeng_group/huangzejian/zz.tools/blocks_connection/blocks_connection_v4_mutil.py -B /vol3/agis/chengshifeng_group/shiyan/03-pea_WGS/003.SNP_recalling_773_V2/09-hapmap/01.LDblock/$a.ldblock.blocks.det -v /vol3/agis/chengshifeng_group/shiyan/03-pea_WGS/003.SNP_recalling_773_V2/09-hapmap/00.phased_vcf/${a}.snp.raw.HARD.Missing-unphasing.ID.allele2_retain.hard_retain.InbreedingCoeff_retain.missing_retain.maf_retain.ann.phased.vcf.gz -b 1000 -q 25 -c 0.98 -o /vol3/agis/chengshifeng_group/shiyan/03-pea_WGS/003.SNP_recalling_773_V2/09-hapmap/02.LDblock_identify/output">${a}.sh
done
