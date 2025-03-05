sentieon_path=sentieon-genomics-202112.01
data_path=./
ref_genome=./ZW6_ref.fa

for sample in `cat sample.list`;
do
mkdir $sample
all=""
nt=8
workdir=`pwd`
exec >$workdir/${i}.log 2>&1
cd ${workdir}

echo "#!/bin/sh
# ******************************************
export SENTIEON_LICENSE=10.1.0.2:8999
# ******************************************

# ******************************************
# 1. Mapping reads with BWA-MEM, sorting
# ******************************************
">${sample}/${sample}.sentieon.sh

for library in `grep $sample sample_lib.list | cut -f 2 -d ","`;
do

fq1=$data_path/$sample/${sample}_${library}_1.fastp.fastq.gz
fq2=$data_path/$sample/${sample}_${library}_2.fastp.fastq.gz
group=$library
platform="ILLUMINA"
bam_option="--bam_compression 1"

echo "time ( $sentieon_path/bin/sentieon bwa mem -M -R '@RG\tID:$group\tPL:$platform\tLB:$library\tSM:$sample' -t $nt \
-K 10000000 $ref_genome $fq1 $fq2 || echo -n 'error' ) | $sentieon_path/bin/sentieon util sort $bam_option \
-r $ref_genome -o ${sample}.${library}.sort.bam -t $nt --sam2bam -i -
">>${sample}/${sample}.sentieon.sh
all="$all ${sample}.${library}.sort.bam"
done

multi=2
runs=$(grep $sample sample_lib.list | cut -f 2 -d ","|wc -l)
if [ $runs -ge $multi ];then

echo "samtools merge -@ $nt ${sample}.sorted.bam $all && rm $all ${all}.bai

samtools index -@ $nt ${sample}.sorted.bam
" >>${sample}/${sample}.sentieon.sh
else

echo "mv ${sample}.${library}.sort.bam ${sample}.sorted.bam
rm ${sample}.${library}.sort.bam.bai
samtools index -@ $nt ${sample}.sorted.bam
" >>${sample}/${sample}.sentieon.sh
fi

echo "
# ******************************************
# 2. Metrics
# ******************************************
time $sentieon_path/bin/sentieon driver -r $ref_genome -t $nt -i ${sample}.sorted.bam --algo MeanQualityByCycle ${sample}_mq_metrics.txt \
--algo QualDistribution ${sample}_qd_metrics.txt --algo GCBias --summary ${sample}_gc_summary.txt ${sample}_gc_metrics.txt --algo AlignmentStat \
--adapter_seq '' ${sample}_aln_metrics.txt --algo InsertSizeMetricAlgo ${sample}_is_metrics.txt

$sentieon_path/bin/sentieon plot GCBias -o ${sample}_gc-report.pdf ${sample}_gc_metrics.txt
$sentieon_path/bin/sentieon plot QualDistribution -o ${sample}_qd-report.pdf ${sample}_qd_metrics.txt
$sentieon_path/bin/sentieon plot MeanQualityByCycle -o ${sample}_mq-report.pdf ${sample}_mq_metrics.txt
$sentieon_path/bin/sentieon plot InsertSizeMetricAlgo -o ${sample}_is-report.pdf ${sample}_is_metrics.txt

# ******************************************
# 3. Remove Duplicate Reads
# ******************************************
time $sentieon_path/bin/sentieon driver -t $nt -i ${sample}.sorted.bam --algo LocusCollector --fun score_info ${sample}_score.txt
time $sentieon_path/bin/sentieon driver -t $nt -i ${sample}.sorted.bam --algo Dedup --rmdup --score_info ${sample}_score.txt --metrics ${sample}_dedup_metrics.txt $bam_option ${sample}.deduped.bam

# ******************************************
# 4. Indel realigner
# ******************************************
#time $sentieon_path/bin/sentieon driver -r $ref_genome -t $nt -i ${sample}.deduped.bam --algo Realigner ${sample}.realigned.bam

# ******************************************
# 5. HC Variant caller
# ******************************************
time $sentieon_path/bin/sentieon driver -r $ref_genome -t $nt -i ${sample}.deduped.bam  --algo Haplotyper --emit_mode gvcf ${sample}.gvcf.gz
samtools view -h -T ./ZW6_ref.fa -C -@ 8 ${sample}.deduped.bam -o ${sample}.deduped.cram
rm ${sample}.deduped.bam ${sample}.deduped.bam.bai
">> ${sample}/${sample}.sentieon.sh
done

