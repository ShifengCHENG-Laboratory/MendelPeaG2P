hisat2 -x ZW6_ref  -q -p 20 -1 ../../00.QC/${i}_R1.retain.fq.gz -2 ../../00.QC/${i}_R2.retain.fq.gz -S ../ZW6/${i}.ZW6.sam
samtools view -bS ../ZW6/${i}.ZW6.sam | samtools sort -@ 8 -o ../ZW6/${i}.ZW6.sorted.bam
rm ../ZW6/${i}.ZW6.sam
samtools index -@ 8 ../ZW6/${i}.ZW6.sorted.bam
stringtie ../ZW6/${i}.ZW6.sorted.bam -e -p 8 -G ZW6_chromosome_genes.part.gff3 -o ../ZW6/${i}.ZW6.gtf -l ../ZW6/${i}.ZW6 -A ../ZW6/${i}.ZW6.gene.txt
