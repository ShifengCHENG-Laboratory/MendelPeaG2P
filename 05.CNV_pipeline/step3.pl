my%hash;
open IN1,"$ARGV[0]"; ### region deep
while(<IN1>){
        chomp;
        my@RD=(split /\t/,$_);
        my$name=$RD[0];
        for (my$i=1;$i<698;$i+=1){
                if($RD[$i]>1.6){
                        push @{$hash{$name}},$i+8;
                }
        }
}

open IN2,"$ARGV[1]"; ### gff
my%gff;
while(<IN2>){
        chomp;
        if(/\tgene\t/){
                my($gene_gff,$start)=(split /\t/,$_)[8,3];
                $gene_gff=(split /;/,$gene_gff)[0];
                $gene_gff=(split /=/,$gene_gff)[1];
                $gff{$gene_gff}=$start;
        }
}

my$case_line="AA\tAA\tAA\tATCG\tA\t.\tPASS\t.\tGT";
for(my$i=0;$i<697;$i+=1){
        $case_line="$case_line\t0|0";
}
my$count;
foreach my$key ( sort keys %hash ){
        if(exists $hash{$key}){
                $count+=1;
                my%reduce;
                my@mutation = grep { ++$reduce{$_} < 2 } @{$hash{$key}};
                my@line=(split /\t/,$case_line);
                my$chr=$key;
                $chr=~s/gene-KIW84_0//;
                $chr=substr($chr, 0, 1);
                $line[0]=$chr;
                if(exists $gff{$key}){
                        $line[1]=$gff{$key};
                }else{
                        $line[1]=$count;
                }
                $line[2]=$key;
                foreach my$num (@mutation){
                        $line[$num]="1|1";
                }
                print join "\t",@line;
                print "\n";
        }
}

