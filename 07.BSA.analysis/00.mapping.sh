############raw data quality control
#!/usr/bin/bash
for sample in `cat sample.list`
do
#mkdir -p $sample
for library in `grep -w $sample sample_lib.list | cut -f 2 -d ","`
do
echo "/vol3/agis/chengshifeng_group/shiyan/0-software/bin/fastp -i /vol3/agis/chengshifeng_group/chenglab/02.pea/00.raw_data/00-raw_data/${library}_1.fq.gz -I /vol3/agis/chengshifeng_group/chenglab/02.pea/00.raw_data/00-raw_data/${library}_2.fq.gz -o ${sample}_${library}_1.fastp.fastq.gz -O ${sample}_${library}_2.fastp.fastq.gz --adapter_sequence=AGTCGGAGGCCAAGCGGTCTTAGGAAGACAA --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA -w 5 -f 5 -F 5 -l 80 -g -j ${sample}_${library}.json -h ${sample}_${library}.html" > ${sample}_${library}.fastq.sh
done

#!/usr/bin/bash
wrfp="/vol3/agis/chengshifeng_group/shiyan/03-pea_WGS/ZW6_ref/01-index/ZW6_ref.fa"

for sample in `cat sample.list`
do
mkdir -p $sample
for library in `grep -w $sample sample_lib.list | cut -f 2 -d ","`
do
echo "/vol3/agis/chengshifeng_group/shiyan/0-software/bin/fastp -i ${sample}_${library}_1.fq.gz -I ${sample}_${library}_2.fq.gz -o ./$sample/${sample}_${library}_1.fastp.fastq.gz -O ./$sample/${sample}_${library}_2.fastp.fastq.gz --adapter_sequence=AGTCGGAGGCCAAGCGGTCTTAGGAAGACAA --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA -w 5 -f 5 -F 5 -l 80 -g -j ${sample}_${library}.json -h ${sample}_${library}.html" > ${sample}_${library}.map.sh
done
done

for sample in `cat sample.list`
do
#mkdir -p $sample
for library in `grep -w $sample sample_lib.list | cut -f 2 -d ","`
do
RGID=$library
echo " source ~/.bashrc
wrfp=\"/vol3/agis/chengshifeng_group/shiyan/03-pea_WGS/ZW6_ref/01-index/ZW6_ref.fa\"
sample=$sample
TMP_DIR=`pwd`/tmp

/public/agis/chengshifeng_group/fengcong/WGRS/software/bwa-0.7.17/bwa mem -t 20 -R \"@RG\\tID:$RGID\\tPL:ILLUMINA\\tLB:$library\\tSM:$sample\" $wrfp ./$sample/${sample}_${library}_1.fastp.fastq.gz ./$sample/${sample}_${library}_2.fastp.fastq.gz | samtools sort -T ${sample}.${library} -m 1G -@ 10 - > ./${sample}/${sample}.${library}.sorted.bam

java -Xmx20g -XX:ParallelGCThreads=20  -Djava.io.tmpdir=`pwd`/tmp -jar /public/agis/chengshifeng_group/fengcong/WGRS/software/picard/picard.jar MarkDuplicates \
I=./${sample}/${sample}.${library}.sorted.bam \
M=./${sample}/${sample}.${library}.markdup_metrics.txt \
O=./${sample}/${sample}.${library}.sorted.markdup.bam

samtools index ./${sample}/${sample}.${library}.sorted.markdup.bam

" >> ${sample}_${library}.map.sh
done
done
