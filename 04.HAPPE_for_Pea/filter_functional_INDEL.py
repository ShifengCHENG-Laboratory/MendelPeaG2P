# filter functional SNPs
# usage: python3 filter_functional_SNP.py <input_file> | bgzip -c > <output_file>

import sys
import gzip
from typing_extensions import Annotated
import  INDELANN_V4
functional_or_not_d = {"bidirectional_gene_fusion":1,
"gene_fusion":1,
"exon_variant":1,
"intragenic_variant":1,
"splice_region_variant":1,
"intron_variant":0,
"upstream_variant":0,
"frameshift_variant":1,
"downstream_gene_variant":0,
"intergenic_region":0,
"regulatory_region_variant":0,
"5_prime_UTR_variant":0,
"3_prime_UTR_variant":0,
}

if __name__ == "__main__":
    input_file = sys.argv[1]
    inf = gzip.open(input_file, 'rt') if input_file.endswith('.gz') else open(input_file, 'r')
    outf = sys.stdout
    for line in inf:
        if line.startswith('#'):
            outf.write(line)
        else:
            ls = line.strip().split()
            featureid,annotation,c_xx,p_xx = INDELANN_V4.INDEL_ANN_of_Gene_Structure(ls[7])
            if annotation  not in functional_or_not_d:
                sys.stderr.write(line)
                sys.stderr.write(annotation)
#                sys.stderr.write("function annotaion error\n")
#                exit(-1)
            if functional_or_not_d[annotation]:
                outf.write(line)
            

    inf.close()
    outf.close()
