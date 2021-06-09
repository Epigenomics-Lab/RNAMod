#! perl -w
## call high confidence m5C sites between two replicates
## usage: perl call_m5c_highconf.pl -r1 mRNA_IVT_rep1_fail_convert_fisher -r2 mRNA_IVT_rep2_fail_convert_fisher -c 5 -d 10 -r 0.05 -f 0.2 -p 0.05 > m5c_highconf_sites

use Getopt::Long;
GetOptions("r1=s" => \$rep1,"r2=s" => \$rep2, "d=s" => \$dep, "c=s" => \$c,"r=s" => \$rate, "f=s" => \$fpr, "p=s" => \$pvalue) or die ("Invalid arguments.\n");
open (F1, "$rep1") || die "Cant open file: $!";
open (F2, "$rep2") || die "Cant open file: $!";

while(<F1>){
        if($_!~/WARNING/){
        s/\s+$//;
        @tmp1=split/\s+/,$_;
        $line1=join("\t",$tmp1[0],$tmp1[1],$tmp1[2]);
        $an{$line1}=join("\t",$tmp1[3],$tmp1[4]);
        push @ids,$line1;
        if($tmp1[8]>=$c && $tmp1[9]>=$dep && $tmp1[6]>=$dep && $tmp1[10]>=$rate && $tmp1[11]<=$fpr && $tmp1[12]<=$pvalue){
        $hash1{$line1}=$tmp1[10]-$tmp1[7];
	$num1{$line1}=join("\t",$tmp1[5],$tmp1[6],$tmp1[8],$tmp1[9]);
      }
}
}

while(<F2>){
        if($_!~/WARNING/){
        s/\s+$//;
        @tmp2=split/\s+/,$_;
        $line2=join("\t",$tmp2[0],$tmp2[1],$tmp2[2]);
        $an{$line2}=join("\t",$tmp2[3],$tmp2[4]);
        push @ids,$line2;
        if($tmp2[8]>=$c && $tmp2[9]>=$dep && $tmp2[6]>=$dep && $tmp2[10]>=$rate && $tmp2[11]<=$fpr && $tmp2[12]<=$pvalue){
        $hash2{$line2}=$tmp2[10]-$tmp2[7];
	$num2{$line2}=join("\t",$tmp2[5],$tmp2[6],$tmp2[8],$tmp2[9]);
      }
}
}

undef %count;
@uniq_ids = grep { ++$count{ $_ } < 2; } @ids;


for $id(@uniq_ids){
        if(exists $hash1{$id} && exists $hash2{$id}){
		@rep1=split/\s+/,$num1{$id};
		@rep2=split/\s+/,$num2{$id};
		$methy=sprintf "%0.4f",($rep1[2]+$rep2[2])/($rep1[3]+$rep2[3])-($rep1[0]+$rep2[0])/($rep1[1]+$rep2[1]);
        	print "$id\t$an{$id}\t$hash1{$id}\t$hash2{$id}\t$methy\n";
        }
      }
      
close F1;
close F2;

