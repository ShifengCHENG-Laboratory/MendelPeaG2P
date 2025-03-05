java -Xmx20g -jar /software/snpEff/snpEff.jar \
            -v Pisum_sativum_ZW6 \
            /zz.all.raw.HARD.Missing-unphasing.ID.vcf.gz \
            | bgzip -c > zz.all.raw.HARD.Missing-unphasing.ID.ann.vcf.gz
            cstat=${PIPESTATUS[@]};stats=${stats}":""${cstat}" && echo QSstats_${stats_num}:${cstat} && let stats_num+=1  ##--/--##

            /bin/bcftools index zz.all.raw.HARD.Missing-unphasing.ID.ann.vcf.gz

