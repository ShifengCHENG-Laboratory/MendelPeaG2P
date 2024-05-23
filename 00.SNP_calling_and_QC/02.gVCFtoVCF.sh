ref=/public/home/shiyan/002-pea_ref/03-ZW6/ZW6_ref.fa
gvcf_argument=""
while read -r line;
do
gvcf_argument=$gvcf_argument" -v $line"
done < "zz_697_gvcf.list"
echo "#!/bin/sh
export SENTIEON_LICENSE=10.1.0.2:8999

time /home/songbo/software/sentieon-genomics-202112.01/bin/sentieon driver -r $ref --algo GVCFtyper $gvcf_argument zz.697.raw.vcf.gz" > zz.697.GVCFtyper.sh

