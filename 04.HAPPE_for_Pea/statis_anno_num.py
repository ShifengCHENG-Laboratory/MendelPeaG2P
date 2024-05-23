# coding=utf-8
import SNPANN_upgrade_3
import sys
import gzip
import re

various_snp_count_d = {
    "3_prime_UTR_variant":0,
    "start_gain":0,
    "5_prime_UTR_variant":0,
    "start_gain":0,
    "downstream_gene_variant":0,
    "intergenic_region":0,
    "intron_variant":0,
    "splice_acceptor_variant":0,
    "splice_donor_variant":0,
    "splice_region_variant":0,
    "splice_acceptor_variant&splice_donor_variant":0,
    "splice_acceptor_variant":0,
    "stop_retained_variant":0,
    "start_lost":0,
    "stop_lost":0,
    "missense_variant":0,
    "start_lost":0,
    "stop_gained":0,
    "stop_gained":0,
    "stop_lost":0,
    "stop_retained_variant":0,
    "synonymous_variant":0,
    "upstream_gene_variant":0
}

if __name__ == "__main__":
    vcf_file = gzip.open(sys.argv[1],
                         'rt') if sys.argv[1].endswith('.gz') else open(
                             sys.argv[1], 'r')
    line = vcf_file.readline()
    anno_dic = {} #{anno1:12,anno2:18.....}
    while line:
        if re.match('#', line):
          line = vcf_file.readline()
     else:
        INFO_field = line.strip().split('\t')[7]
       anno = (SNPANN_upgrade.SNP_ANN_of_Gene_Structure(INFO_field))
      if anno in anno_dic:
                anno_dic[anno] += 1
            else:
                anno_dic[anno] = 1
            line = vcf_file.readline()
    for gene_str in anno_dic:
        sys.stdout.write(gene_str + '\t' + str(anno_dic[gene_str]) + '\n')
    #while line:
    #    if re.match('#', line):
    #        line = vcf_file.readline()
    #    else:
    #        INFO_field = line.strip().split('\t')[7]
    #        anno_tuple = (SNPANN_upgrade_3.SNP_ANN_of_Gene_Structure(INFO_field))
    #        various_snp_count_d[anno_tuple[1]] += 1
    #        line = vcf_file.readline()

    #for key in various_snp_count_d:
    #    sys.stdout.write(key+"\t"+str(various_snp_count_d[key])+"\n")

    vcf_file.close()

