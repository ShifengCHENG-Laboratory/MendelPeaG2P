/public/agis/chengshifeng_group/fengcong/WGRS/software/plink/plink --vcf /vol3/agis/chengshifeng_group/shiyan/3-pea/01.population/01.high-SNP/SNP.4dtv.vcf.gz \
            --indep-pairwise 100 50 0.2 \
            --out snp.clean \
            --allow-extra-chr \
            --make-bed

            /public/agis/chengshifeng_group/fengcong/WGRS/software/plink/plink --bfile snp.clean --extract snp.clean.prune.in --out prunData --recode 12 --allow-extra-chr


            /public/agis/chengshifeng_group/fengcong/WGRS/software/admixture_linux-1.3.0/admixture --cv prunData.ped 2 >> log.txt

            for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; do (/public/agis/chengshifeng_group/fengcong/WGRS/software/admixture_linux-1.3.0/admixture --cv prunData.ped $K | tee log${K}.out)&  done

            wait

