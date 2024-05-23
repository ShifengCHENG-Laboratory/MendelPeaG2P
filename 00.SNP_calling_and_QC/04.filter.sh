########################SNP filter#######################
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7
do
    datapath="/vol3/agis/chengshifeng_group/shiyan/03-pea_WGS/003.SNP_recalling_773_V2/02-filter/00.allele2_filter"
    pythonpath="/public/home/fengcong/anaconda2/envs/py3/bin/"
    echo "
    #get hardfilter
    mkdir -p ../00.allele2_filter
    $pythonpath/python3 /public/agis/chengshifeng_group/fengcong/WGRS/software/Fc-code/snp_qc_script/filter_allele2.py /vol3/agis/chengshifeng_group/shiyan/03-pea_WGS/003.SNP_recalling_773_V2/00-raw_vcf/01-HARD_ID/$chr.SNP.raw.HARD.Missing-unphasing.ID.vcf.gz ../00.allele2_filter/$chr.snp.raw.HARD.Missing-unphasing.ID

    mkdir -p ../01.hardfilter;
    $pythonpath/python3 /public/agis/chengshifeng_group/fengcong/WGRS/software/Fc-code/snp_qc_script/filter_hard.py $datapath/$chr.snp.raw.HARD.Missing-unphasing.ID.allele2_retain.vcf.gz ../01.hardfilter/$chr.snp.raw.HARD.Missing-unphasing.ID.allele2_retain
    /public/agis/chengshifeng_group/fengcong/WGRS/software/bcftools1.9/bin/bcftools index ../01.hardfilter/$chr.snp.raw.HARD.Missing-unphasing.ID.allele2_retain.hard_filter.vcf.gz
    /public/agis/chengshifeng_group/fengcong/WGRS/software/bcftools1.9/bin/bcftools index ../01.hardfilter/$chr.snp.raw.HARD.Missing-unphasing.ID.allele2_retain.hard_retain.vcf.gz

    #get IC info
    mkdir -p ../02.calc_IC_info;
    $pythonpath/python3 /public/agis/chengshifeng_group/fengcong/WGRS/software/Fc-code/calc_inbreeding_coefficient_v3.py ../01.hardfilter/$chr.snp.raw.HARD.Missing-unphasing.ID.allele2_retain.hard_retain.vcf.gz > ../02.calc_IC_info/$chr.snp.InbreedingCoeff.info.txt
    Fmedian=\`tail -n2  ../02.calc_IC_info/$chr.snp.InbreedingCoeff.info.txt | grep \"#Fmedian\" | cut -f 2 -d \"=\"\`;

    #filter IC
    mkdir -p ../03.IC_filter;
    $pythonpath/python3 /public/agis/chengshifeng_group/fengcong/WGRS/software/Fc-code/filter_according_inbreeding_coefficient_v2.py ../01.hardfilter/$chr.snp.raw.HARD.Missing-unphasing.ID.allele2_retain.hard_retain.vcf.gz \$Fmedian ../03.IC_filter/$chr.snp.raw.HARD.Missing-unphasing.ID.allele2_retain.hard_retain
    /public/agis/chengshifeng_group/fengcong/WGRS/software/bcftools1.9/bin/bcftools index ../03.IC_filter/$chr.snp.raw.HARD.Missing-unphasing.ID.allele2_retain.hard_retain.InbreedingCoeff_retain.vcf.gz
    /public/agis/chengshifeng_group/fengcong/WGRS/software/bcftools1.9/bin/bcftools index ../03.IC_filter/$chr.snp.raw.HARD.Missing-unphasing.ID.allele2_retain.hard_retain.InbreedingCoeff_filter.vcf.gz

    #filter missing
    mkdir -p ../04.missing_filter;
    $pythonpath/python3 /public/agis/chengshifeng_group/fengcong/WGRS/software/Fc-code/missing_filter.py ../03.IC_filter/$chr.snp.raw.HARD.Missing-unphasing.ID.allele2_retain.hard_retain.InbreedingCoeff_retain.vcf.gz ../04.missing_filter/$chr.snp.raw.HARD.Missing-unphasing.ID.allele2_retain.hard_retain.InbreedingCoeff_retain
    /public/agis/chengshifeng_group/fengcong/WGRS/software/bcftools1.9/bin/bcftools index ../04.missing_filter/$chr.snp.raw.HARD.Missing-unphasing.ID.allele2_retain.hard_retain.InbreedingCoeff_retain.missing_retain.vcf.gz
    /public/agis/chengshifeng_group/fengcong/WGRS/software/bcftools1.9/bin/bcftools index ../04.missing_filter/$chr.snp.raw.HARD.Missing-unphasing.ID.allele2_retain.hard_retain.InbreedingCoeff_retain.missing_filter.vcf.gz

    #filter maf
    mkdir -p ../05.maf_filter;
    $pythonpath/python3 /public/agis/chengshifeng_group/fengcong/WGRS/software/Fc-code/maf_filter.py ../04.missing_filter/$chr.snp.raw.HARD.Missing-unphasing.ID.allele2_retain.hard_retain.InbreedingCoeff_retain.missing_retain.vcf.gz ../05.maf_filter/$chr.snp.raw.HARD.Missing-unphasing.ID.allele2_retain.hard_retain.InbreedingCoeff_retain.missing_retain
    /public/agis/chengshifeng_group/fengcong/WGRS/software/bcftools1.9/bin/bcftools index ../05.maf_filter/$chr.snp.raw.HARD.Missing-unphasing.ID.allele2_retain.hard_retain.InbreedingCoeff_retain.missing_retain.maf_retain.vcf.gz
    /public/agis/chengshifeng_group/fengcong/WGRS/software/bcftools1.9/bin/bcftools index ../05.maf_filter/$chr.snp.raw.HARD.Missing-unphasing.ID.allele2_retain.hard_retain.InbreedingCoeff_retain.missing_retain.maf_filter.vcf.gz




    " > $chr.filter.sh


done


#############INDEL filter###########################
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7
do
mkdir -p ../02.HARD_filter
mkdir -p ../03.ICF
    echo "

    Fmedian=\`tail -n2  ../03.ICF/${chr}.InbreedingCoeff.info.txt | grep \"#Fmedian\" | cut -f 2 -d \"=\"\`;

    python3 filter_according_inbreeding_coefficient_v2.py ../02.HARD_filter/${chr}.INDEL.raw.HARD.Missing-unphasing.ID.allele2_retain.hard_retain.vcf.gz \$Fmedian ../03.ICF/${chr}.INDEL.raw.HARD.Missing-unphasing.ID.ann.allele2_retain.hard_retain
    /public/agis/chengshifeng_group/fengcong/WGRS/software/bcftools1.9/bin/bcftools index ../03.ICF/${a}*filter.vcf.gz
    /public/agis/chengshifeng_group/fengcong/WGRS/software/bcftools1.9/bin/bcftools index ../03.ICF/${a}*retain.vcf.gz
    " > ${chr}_filter.sh

done

