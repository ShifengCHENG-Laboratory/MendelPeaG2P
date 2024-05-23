import sys
import gzip
from bgzip import BGZipWriter
# python __file__ in.vcf.gz out_prefix 

inf = gzip.open(sys.argv[1],"rt")


        
class TXTBGZipWriter(BGZipWriter):
    def write(self, text):
        super().write(text.encode("utf-8"))
out_prefix = sys.argv[2]
ouf1 = open("%s.allele2_filter.vcf.gz"%(out_prefix),"wb")
ouf2 = open("%s.allele2_retain.vcf.gz"%(out_prefix),"wb")
w1 = TXTBGZipWriter(ouf1)
w2 = TXTBGZipWriter(ouf2)
filter_count,retain_count = 0,0
#read header line
line = inf.readline()
while line.startswith("#"):
    w1.write(line)
    w2.write(line)
    line = inf.readline()

## dict
multi_allele_dict = {
    2:{"N":0 , "*":0},
    3:{"N":0 , "*":0},
    4:{"N":0 , "*":0},
    5:{"N":0 , "*":0},
    6:{"N":0 , "*":0},
    7:{"N":0 , "*":0},
    8:{"N":0 , "*":0},
    9:{"N":0 , "*":0},
    10:{"N":0 , "*":0}
    }

while line :
    ls = line.split()
    if ls[4].find(",") != -1:
        w1.write(line)
        filter_count += 1
        allele_ls = ls[4].split(",")
        if "*" in allele_ls:
            multi_allele_dict[len(allele_ls)]["*"] += 1
        else:
            multi_allele_dict[len(allele_ls)]["N"] += 1
    else:
        w2.write(line)
        retain_count += 1

    line = inf.readline()


sys.stderr.write("#filter: %d\n#retain: %d\n"%(filter_count,retain_count))
for i in range(2,11):
    sys.stderr.write("multi_allele_%d\ttotal:%d\tN:%d\t*:%d\n"%(i, sum(multi_allele_dict[i].values()) ,multi_allele_dict[i]["N"],multi_allele_dict[i]["*"]))


inf.close()
w1.close()
w2.close()
