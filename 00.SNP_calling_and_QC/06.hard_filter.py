import sys
import gzip
from bgzip import BGZipWriter
# python __file__ in.vcf.gz out_prefix 

inf = gzip.open(sys.argv[1],"rt")


        
class TXTBGZipWriter(BGZipWriter):
    def write(self, text):
        super().write(text.encode("utf-8"))
out_prefix = sys.argv[2]
ouf1 = open("%s.hard_filter.vcf.gz"%(out_prefix),"wb")
ouf2 = open("%s.hard_retain.vcf.gz"%(out_prefix),"wb")
w1 = TXTBGZipWriter(ouf1)
w2 = TXTBGZipWriter(ouf2)
filter_count,retain_count = 0,0
#read header line
line = inf.readline()
while line.startswith("#"):
    w1.write(line)
    w2.write(line)
    line = inf.readline()

while line :
    ls = line.split()
    if ls[6] != "PASS":
        w1.write(line)
        filter_count += 1
    else:
        w2.write(line)
        retain_count += 1

    line = inf.readline()


sys.stderr.write("#filter: %d\n#retain: %d\n"%(filter_count,retain_count))



inf.close()
w1.close()
w2.close()
