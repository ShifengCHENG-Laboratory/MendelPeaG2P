# Step 0 - Input
Each SMRT cell will get a file  "final.bam" and primer sequence "primer.fasta" from company


# Step 1 - Refine
`isoseq3 refine final.bam primer.fasta output.flnc.bam --require-polya -j 5`

# Step 2 - refined transcripts clustering
Merged more than one SMRT cells

`ls output1.lnc.bam  output2.flnc.bam output3.flnc.bam > output.fofn`

Clustering 

`isoseq cluster  output.fofn output.clustered.bam --verbose --use-qvs -j 5`

# Step 3 - align
`pbmm2 align ZW6_ref.mmi output.clustered.bam output.sortede.bam -j 10 --sort --preset ISOSEQ`

# Step 4 - estimated expression 

`collapse -j 4 output.sorted.bam output.gff`

`gffread output.gff -o output.gff3`
