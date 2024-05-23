use warnings;
open IN1,"$ARGV[0]";
# my$header=<IN1>;
# print "$header";
while(<IN1>){
	chomp;
	my@accession=(split /\t/,$_);
	my$num=@accession;
	for (my$i=1;$i<$num;$i+=1){
		if($accession[$i] < 0.1){
			$accession[$i]="0";
		}elsif( 0.1<=$accession[$i] && $accession[$i] <0.6){
			$accession[$i]="0.5";
		}elsif(0.6 <= $accession[$i] && $accession[$i] < 1.6){
			$accession[$i]="1";
		}elsif(1.6 <= $accession[$i] && $accession[$i] < 2.6){
			$accession[$i]="2";
		}elsif(2.6<= $accession[$i] && $accession[$i]< 3.6){
			$accession[$i]="3";
		}elsif(3.6<=$accession[$i]  && $accession[$i]< 4.6){
			$accession[$i]="4";
#		}elsif(4.6<=$accession[$i] && $accession[$i]< 5.6){
#			$accession[$i]="5";
#		}elsif(5.6<=$accession[$i] && $accession[$i]< 6.6){
#			$accession[$i]="6";
#		}elsif(6.6<=$accession[$i]  && $accession[$i]< 7.6){
#			$accession[$i]="7";
#		}elsif(7.6<=$accession[$i] && $accession[$i]< 8.6){
#			$accession[$i]="8";
#		}elsif(8.6<=$accession[$i] && $accession[$i] < 9.6){
#			$accession[$i]="9";
#		}elsif(9.6<=$accession[$i] && $accession[$i]< 10.6){
#			$accession[$i]="10";
		}else{
			$accession[$i]="5";
		}
	}
	print join "\t",@accession;
	print "\n";
}		
