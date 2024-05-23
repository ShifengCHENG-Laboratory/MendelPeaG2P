#2021-4-8 14:18:23
#usage:python __file__ in.vcf/in.vcf.gz outprefix > stat.txt
#output: prefix.missing_filter.vcf.gz ; prefix.missing_retain.vcf.gz

import sys,re
import gzip
# import numpy as np
import bgzip 
from bgzip import BGZipWriter

def search_pattern(s,pattern):
    '''
    return string
    '''
    m  = re.search(pattern,s,flags=0)
    if m != None:
        # print(m.groups())
        return m.groups()[0]
    else:
        return "NA"

class TXTBGZipWriter(BGZipWriter):
    def write(self, text):
        super().write(text.encode("utf-8"))

if __name__ == "__main__":
    inf = gzip.open(sys.argv[1],"rt") if sys.argv[1].endswith(".gz") else open(sys.argv,"r")
    outprefix = sys.argv[2]

    ouf1 = open("%s.maf_filter.vcf.gz"%(outprefix),"wb")
    ouf2 = open("%s.maf_retain.vcf.gz"%(outprefix),"wb")
    w1 = TXTBGZipWriter(ouf1)
    w2 = TXTBGZipWriter(ouf2)

    #read header line
    line = inf.readline()
    while line.startswith("#"):
        w1.write(line)
        w2.write(line)
        line = inf.readline()

    # total_sample_count = len(line.split()[9:])
    # print(total_sample_count)

    non_reference_allele_frq_pattern = re.compile(r'AF=(-?\d+(\.\d+)?[eE]?[+-]?(\d+)?)')


    retain_count = 0
    filter_count = 0

    while line:
        # sample_count = 0
        # het_sample_count = 0
        # het_sample_frq = 0
        # InbreedingCoeff_calc = 0

        ls = line.strip().split()
        non_reference_allele_frq = search_pattern(ls[7],non_reference_allele_frq_pattern)

        # sample_count = int(int(allele_count)/2)
        # print(sample_count/float(total_sample_count))
        maf = min(float(non_reference_allele_frq),1-float(non_reference_allele_frq))

        if maf >=  0.01:
            w2.write(line )
            retain_count += 1
        else:
            w1.write(line)
            filter_count += 1

        

        

        line = inf.readline()

    #output Fmedian
    # sys.stdout.write("#using SNP sites that InbreedingCoeff_calc>0 & maf >= 0.05 to estimate Fmedian:%d SNPs\n"%(len(InbreedingCoeff_calc_list_for_Fmedian)))
    sys.stdout.write("#filter: %d\n#retain: %d\n"%(filter_count,retain_count))


    inf.close()

    w1.close()
    w2.close()