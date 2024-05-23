import sys
import os
import argparse
import logging

haplotype_matrix_file = \
    "/vol3/agis/chengshifeng_group/huangzejian/06.wheat_pop_analysis/02.LDblock_identify/01.plot_data/03.functional/09.all_chromosome_blocks_excel_haplotype_matrix/sorted.all.blocks.haplotype.matrix.txt"
bcftools = "/public/agis/chengshifeng_group/fengcong/WGRS/software/bcftools1.9/bin/bcftools"
data_path="/vol3/agis/chengshifeng_group/fengcong/wheat_pop_analysis/40.89mapping_and_merge/14.final_SNP_filter2/03.IC_filter"
order = "/vol3/agis/chengshifeng_group/fengcong/zz.my_jupyter_notebook/haplotype_GWAS/zz.sample_order.list2"
pca_kmean="/vol3/agis/chengshifeng_group/fengcong/zz.my_jupyter_notebook/haplotype_GWAS/pca_kmean.r"
Rscript="/public/home/fengcong/anaconda2/envs/R/bin/Rscript"
header_file="/vol3/agis/chengshifeng_group/fengcong/zz.my_jupyter_notebook/haplotype_GWAS/zz.header.vcf2"
#usage : python __file__ order_file_path haplotype_matrix_file > out.vcf

logging.basicConfig(level = logging.DEBUG,format = '[%(asctime)s] - %(name)s - %(levelname)s - %(message)s',stream = sys.stderr)
logger = logging.getLogger("haplotype GWAS")
def check_ret(ret,step):
    if  ret:
        logger.error("%s -- retrun code:%d\n"%(step,ret))
        exit(-1)

if __name__ == "__main__":
    order_file = sys.argv[1]
    haplotype_matrix_file = sys.argv[2]
    order_list = []
    with open(order_file) as f:
        for line in f:
            order_list.append(line.strip())
    ouf = sys.stdout
    # with open(header_file,"r") as f:
    #     for line in f:
    #         ouf.write(line)


    #read block file
    block_region = [] #[chr, start, end]
    with open(haplotype_matrix_file,"r") as f:
        line = f.readline()
        for line in f:
            blockid = line.strip().split()[0].split("_")
            blockid[1] = int(blockid[1])
            blockid[2] = int(blockid[2])
            block_region.append(blockid[0:3])

    for block in block_region:
        ret = os.system("""
        %s view -r %s:%d-%d  %s/%s.SNP.Missing-unphasing.ID.ann.finalSID.allele2_retain.hard_retain.InbreedingCoeff_retain.vcf.gz | bgzip -c > %s.%d.%d.vcf.gz
        
        """%(bcftools, block[0], block[1], block[2], data_path , block[0], block[0], block[1], block[2])
        )
        check_ret(ret,"%s %d %d"%(block[0], block[1], block[2]))

        prefix = "./%s.%d.%d.vcf.gz"%(block[0], block[1], block[2])
        #pca and kmean
        ret = os.system("""
        %s %s %s 1>/dev/null 2>/dev/null
        """%(Rscript, pca_kmean, prefix)
        )
        check_ret(ret,"%s %d %d pca_kmean"%(block[0], block[1], block[2]))

        cluster_d = {}
        with open(prefix+".cluster","r") as f:
            line = f.readline()
            for line in f:
                ls = line.strip().split(",")
                cluster_d[ls[0]] = ls[1]

        center_d = {}
        with open(prefix+".centers","r") as f:
            line = f.readline()
            for line in f:
                # print(line)
                ls = line.strip().split(",")
                center_d[ls[0]] = float(ls[1])

        gt0 = min(center_d, key=center_d.get)
        gt1 = max(center_d, key=center_d.get)
        all_cluster = set(list(center_d.keys()))
        gth = all_cluster - {gt0,gt1}

        gt_d = {gt0:"0/0",gt1:"1/1",list(gth)[0]:"0/1"}

        out_list = [block[0], block[1], "%s_%d_%d"%(block[0],block[1],block[2]), "C", "T" , "." , "PASS",".", "GT"]
        for sample in order_list:
            if sample in cluster_d:
                out_list.append(gt_d[cluster_d[sample]])
            else:
                out_list.append("./.")
        ouf.write("\t".join(str(_) for _ in out_list)+"\n")

        os.remove("%s"%(prefix))
        os.remove("%s.cluster"%(prefix))
        os.remove("%s.centers"%(prefix))
        os.remove("%s.gds"%(prefix))
        
