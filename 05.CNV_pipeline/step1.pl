if( $ARGV[0] eq "-h" ){
	print "perl step1.pl duplication-list raw-depth summaryfile\n";
	exit;
}

open IN1,"$ARGV[0]";# duplication list
my%hash_dupgene;
while(<IN1>){
	chomp;
	$_=~s/\.1//g;
	my($id)=(split /\t/,$_)[0];
	$hash_dupgene{$id}=$_;
}
close IN1;

open IN2,"gunzip -dc $ARGV[1] | "; # raw depth list
my%hash_genedepth;
while(<IN2>){
	chomp;
	if(/gene/){
		$gene_num+=1;
		my($gene,$depth)=(split /\t/,$_)[3,4];
		$gene=(split /\./,$gene)[0];
		$hash_genedepth{$gene}=$depth;
	}
}
close IN2;

open IN3,"$ARGV[2]"; # summary file;
my$sample_depth;
while(<IN3>){
	if(/total\t/){
		$sample_depth=(split /\t/,$_)[3];
	}
}
close IN3;

foreach my$gene (sort keys %hash_genedepth){ # cat duplication depth
	my$gene_depth=0;
	if( exists $hash_dupgene{$gene}){
		my@gene=split /\t/,$hash_dupgene{$gene};
		foreach my$key (@gene){
			$gene_depth+=$hash_genedepth{$key};
		}
	}else{
		$gene_depth=$hash_genedepth{$gene};
	}
	my$relative_depth=sprintf("%.3f",$gene_depth/$sample_depth);
	print "$gene\t$gene_depth\t$sample_depth\t$relative_depth\n";
}
