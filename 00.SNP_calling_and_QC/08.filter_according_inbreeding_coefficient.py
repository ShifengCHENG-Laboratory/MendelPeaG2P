#!/public/home/fengcong/anaconda2/envs/py3/bin/python
#2021-2-24 17:17:45
#update:2021-3-5 14:55:24
#much much faster than version 1, and ouput can be index using bcftools index
#there will be ouput two files , which one is filter out ,another is retain
#usage:python __file__ in.vcf/in.vcf.gz  Fmedian  out.prefix

import sys,re
import gzip
from bgzip import BGZipWriter

def search_pattern(s,pattern):
    '''
    return string
    '''
    m  = re.search(pattern,s,flags=0)
    if m != None:
        return m.groups()[0]
    else:
        return "NA"
        
class TXTBGZipWriter(BGZipWriter):
    def write(self, text):
        super().write(text.encode("utf-8"))

if __name__ == "__main__":
    inf = gzip.open(sys.argv[1],"rt") if sys.argv[1].endswith(".gz") else open(sys.argv,"r")
    

    Fmedian = float(sys.argv[2])
    out_prefix = sys.argv[3]
    ouf1 = open("%s.InbreedingCoeff_filter.vcf.gz"%(out_prefix),"wb")
    ouf2 = open("%s.InbreedingCoeff_retain.vcf.gz"%(out_prefix),"wb")
    w1 = TXTBGZipWriter(ouf1)
    w2 = TXTBGZipWriter(ouf2)

    #read header line
    line = inf.readline()
    while line.startswith("#"):
        w1.write(line)
        w2.write(line)
        line = inf.readline()
    
    incoeff_gatk_pattern = re.compile(r'InbreedingCoeff=(-?\d+(\.\d+)?[eE]?[+-]?(\d+)?)')
    excehet_pattern = re.compile(r'ExcessHet=(-?\d+(\.\d+)?[eE]?[+-]?(\d+)?)')
    non_reference_allele_count_pattern = re.compile(r'AC=(-?\d+(\.\d+)?[eE]?[+-]?(\d+)?)')
    non_reference_allele_frq_pattern = re.compile(r'AF=(-?\d+(\.\d+)?[eE]?[+-]?(\d+)?)')
    allele_count_pattern = re.compile(r'AN=(-?\d+(\.\d+)?[eE]?[+-]?(\d+)?)')

    filter_count = 0
    retain_count = 0
    while line:
        sample_count = 0
        het_sample_count = 0
        het_sample_frq = 0
        InbreedingCoeff_calc = 0

        ls = line.strip().split()
        non_reference_allele_frq = search_pattern(ls[7],non_reference_allele_frq_pattern)

        for sample in ls[9:]:
            gt = sample.split(":")[0]
            if gt[0] == gt[-1] and gt[0] == ".":
                pass
            elif gt[0] == gt[-1] and gt[0] != ".":
                sample_count +=1
            else: #het
                het_sample_count += 1
                sample_count +=1
        het_sample_frq = (het_sample_count*1.0)/sample_count


        Hexp = (2 * (float(non_reference_allele_frq)) * (1-float(non_reference_allele_frq)) )
        if het_sample_frq > 10*(1-Fmedian)*Hexp:
            filter_count +=1
            w1.write(line)
            # ouf1.write(line)
        else:
            retain_count +=1
            w2.write(line)
            # ouf2.write(line)



        line = inf.readline()

    sys.stdout.write("#filter: %d\n#retain: %d\n"%(filter_count,retain_count))

    w1.close()
    w2.close()
    inf.close()