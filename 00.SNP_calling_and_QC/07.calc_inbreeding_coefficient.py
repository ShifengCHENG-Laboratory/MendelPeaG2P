#2021-2-24 17:17:45
#usage:python __file__ in.vcf/in.vcf.gz > out.txt
#update: add maf calc
import sys,re
import gzip
import numpy as np

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

if __name__ == "__main__":
    inf = gzip.open(sys.argv[1],"rt") if sys.argv[1].endswith(".gz") else open(sys.argv[1],"r")
    #read header line
    line = inf.readline()
    while line.startswith("#"):
        line = inf.readline()

    #list for calc Fmedian
    InbreedingCoeff_calc_list_for_Fmedian = [] 

    #write header
    sys.stdout.write("sample_count\tallele_count\tnon_reference_allele_count\tnon_reference_allele_frq\tmaf\thet_sample_count\thet_sample_rate(nomissing)\tInbreedingCoeff_calc\tInbreedingCoeff_GATK\tExcessHet_GATK\tMQ\tFS\tSOR\tMQRankSum\tReadPosRankSum\tQD\n")

    incoeff_gatk_pattern = re.compile(r'InbreedingCoeff=(-?\d+(\.\d+)?[eE]?[+-]?(\d+)?)')
    excehet_pattern = re.compile(r'ExcessHet=(-?\d+(\.\d+)?[eE]?[+-]?(\d+)?)')
    non_reference_allele_count_pattern = re.compile(r'AC=(-?\d+(\.\d+)?[eE]?[+-]?(\d+)?)')
    non_reference_allele_frq_pattern = re.compile(r'AF=(-?\d+(\.\d+)?[eE]?[+-]?(\d+)?)')
    allele_count_pattern = re.compile(r'AN=(-?\d+(\.\d+)?[eE]?[+-]?(\d+)?)')

    mq_pattern = re.compile(r'MQ=(-?\d+(\.\d+)?[eE]?[+-]?(\d+)?)')
    fs_pattern = re.compile(r'FS=(-?\d+(\.\d+)?[eE]?[+-]?(\d+)?)')
    sor_pattern = re.compile(r'SOR=(-?\d+(\.\d+)?[eE]?[+-]?(\d+)?)')
    mqranksum_pattern = re.compile(r'MQRankSum=(-?\d+(\.\d+)?[eE]?[+-]?(\d+)?)')
    readsposranksum_pattern = re.compile(r'ReadPosRankSum=(-?\d+(\.\d+)?[eE]?[+-]?(\d+)?)')
    qd_pattern =re.compile(r'QD=(-?\d+(\.\d+)?[eE]?[+-]?(\d+)?)')


    while line:
        sample_count = 0
        het_sample_count = 0
        het_sample_frq = 0
        InbreedingCoeff_calc = 0

        ls = line.strip().split()
        incoeff_gatk = search_pattern(ls[7],incoeff_gatk_pattern)
        excehet_gatk = search_pattern(ls[7],excehet_pattern)
        non_reference_allele_count = search_pattern(ls[7],non_reference_allele_count_pattern)
        non_reference_allele_frq = search_pattern(ls[7],non_reference_allele_frq_pattern)
        allele_count = search_pattern(ls[7],allele_count_pattern)

        mq = search_pattern(ls[7],mq_pattern)
        fs = search_pattern(ls[7],fs_pattern)
        sor = search_pattern(ls[7],sor_pattern)
        mqranksum = search_pattern(ls[7],mqranksum_pattern)
        readposranksum = search_pattern(ls[7],readsposranksum_pattern)
        qd = search_pattern(ls[7],qd_pattern)

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

        # F = 1 - Het_frq / 2pq
        InbreedingCoeff_calc = 0
        try :
            InbreedingCoeff_calc = 1 - (het_sample_frq / (2 * (float(non_reference_allele_frq)) * (1-float(non_reference_allele_frq)) ))
        except :
            InbreedingCoeff_calc = -0.01

        #maf = minor allele frequency
        maf = min(float(non_reference_allele_frq),1-float(non_reference_allele_frq))

        #F>0 and maf >= 0.05 will be use to calc Fmedian
        if InbreedingCoeff_calc >0 and maf >= 0.05:
            InbreedingCoeff_calc_list_for_Fmedian.append(InbreedingCoeff_calc)
        
        
        sys.stdout.write("%d\t%s\t%s\t%s\t%.4f\t%d\t%.2f\t%.2f\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(sample_count,allele_count,non_reference_allele_count,non_reference_allele_frq,maf,
            het_sample_count,het_sample_frq,InbreedingCoeff_calc,incoeff_gatk,excehet_gatk,
            mq,fs,sor,mqranksum,readposranksum,qd))



        line = inf.readline()

    #output Fmedian
    sys.stdout.write("#using SNP sites that InbreedingCoeff_calc>0 & maf >= 0.05 to estimate Fmedian:%d SNPs\n"%(len(InbreedingCoeff_calc_list_for_Fmedian)))
    sys.stdout.write("#Fmedian=%.3f\n"%(np.median(InbreedingCoeff_calc_list_for_Fmedian)))


    inf.close()