ample in `cat zz_samle.list`
do
echo "/vol3/agis/chengshifeng_group/shiyan/0-software/bin/fastp -i ${sample}_1.fastq.gz -I ${sample}_2.fastq.gz -o ${sample}_1.fastp.fastq.gz -O ${sample}_2.fastp.fastq.gz --adapter_sequence=AGTCGGAGGCCAAGCGGTCTTAGGAAGACAA --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA -w 5 -f 5 -F 5 -l 80 -g -j ${sample}.json -h ${sample}.html" > ${sample}.fastp.sh
done
