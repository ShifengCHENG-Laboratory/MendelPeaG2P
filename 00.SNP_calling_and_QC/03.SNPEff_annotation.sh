java -Xmx20g -jar /public/agis/chengshifeng_group/fengcong/WGRS/software/snpEff/snpEff.jar \
            -v Pisum_sativum_ZW6 \
            /vol3/agis/chengshifeng_group/shiyan/03-pea_WGS/002.SNP_recalling/zz.all.raw.HARD.Missing-unphasing.ID.vcf.gz \
            | bgzip -c > zz.all.raw.HARD.Missing-unphasing.ID.ann.vcf.gz
            cstat=${PIPESTATUS[@]};stats=${stats}":""${cstat}" && echo QSstats_${stats_num}:${cstat} && let stats_num+=1  ##--/--##

            /public/agis/chengshifeng_group/fengcong/WGRS/software/bcftools1.9/bin/bcftools index zz.all.raw.HARD.Missing-unphasing.ID.ann.vcf.gz

