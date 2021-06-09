#! perl -w
##call high conf. m6A sites

use Getopt::Long;
GetOptions("r1=s" => \$rep1,"r2=s" => \$rep2,"d=s" => \$dep, "f=s" => \$fpr, "p=s" => \$pvalue) or die ("Invalid arguments.\n");
open (F1, "$rep1") || die "Cant open file: $!";
open (F2, "$rep2") || die "Cant open file: $!";

while(<F1>){
        if($_!~/WARNING/){
        s/\s+$//;
        @tmp1=split/\s+/,$_;
        $line1=join("\t",$tmp1[0],$tmp1[1],$tmp1[2]);
        push @ids,$line1;
        if($tmp1[4]>=$dep && $tmp1[7]>=$dep && $tmp1[9]<=$fpr && $tmp1[10]<=$pvalue){
        $methy1{$line1}=$tmp1[8]-$tmp1[5];
        $num1{$line1}=join("\t",$tmp1[3],$tmp1[4],$tmp1[6],$tmp1[7]);
      }
}
}

while(<F2>){
        if($_!~/WARNING/){
        s/\s+$//;
        @tmp2=split/\s+/,$_;
        $line2=join("\t",$tmp2[0],$tmp2[1],$tmp2[2]);
        push @ids,$line2;
        if($tmp2[4]>=$dep && $tmp2[7]>=$dep && $tmp2[9]<=$fpr && $tmp2[10]<=$pvalue){
        $methy2{$line2}=$tmp2[8]-$tmp2[5];
        $num2{$line2}=join("\t",$tmp2[3],$tmp2[4],$tmp2[6],$tmp2[7]);
      }
}
}

undef %count;
@uniq_ids = grep { ++$count{ $_ } < 2; } @ids;


for $id(@uniq_ids){
        if(exists $methy1{$id} && exists $methy2{$id}){
                @rep1=split/\s+/,$num1{$id};
                @rep2=split/\s+/,$num2{$id};
                $methy=sprintf "%0.4f",($rep1[2]+$rep2[2])/($rep1[3]+$rep2[3])-($rep1[0]+$rep2[0])/($rep1[1]+$rep2[1]);
                print "$id\t$methy1{$id}\t$methy2{$id}\t$methy\n";
        }
      }

