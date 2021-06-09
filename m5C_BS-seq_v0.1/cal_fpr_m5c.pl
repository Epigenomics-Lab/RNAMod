#! perl -w
## calculate false positive rate (fpr) using mRNA and IVT RNA samples.
## usage: perl cal_fpr_m5c.pl -i IVT.m5C.pileups.formatted.txt -n mRNA.m5C.pileups.formatted.txt -c 20 -d 1 > mRNA_IVT.txt

use Getopt::Long;
GetOptions("i=s" => \$ivt,"n=s" => \$mrna,"c=s" => \$cutoff,"d=s" => \$dep_cutoff ) or die ("Invalid arguments.\n");
open (F1, "$ivt") || die "Cant open file: $!";
open (F2, "$mrna") || die "Cant open file: $!";

while(<F1>){
        s/\s+$//;
        if($_!~/\#/){
                if($cutoff eq "None"){
                        ($chr,$loc,$strand,$dep,$fail)=(split/\s+/,$_)[0,1,2,9,12];
                }else{
                        ($chr,$loc,$strand)=(split/\s+/,$_)[0,1,2];
                        @alls=(split/\s+/,$_)[14,15,16,17,18,19,20,21,22,23,24,25];
                        for $all(@alls){
                                ($cut,$count)=(split/\;/,$all)[0,1];
                                if($cut == $cutoff){
                                        ($dep,$fail)=(split/\,/,$count)[1,2];
                                        last;
                                }
                        }
                }
                $line1=join("\t",$chr,$loc,$strand);
                if($dep>=$dep_cutoff){
                $v=sprintf "%0.5f",$fail/$dep;
                $line2=join("\t",$fail,$dep,$v);
                $hash{$line1}=$line2;
        }
}
}

while(<F2>){
        s/\s+$//;
        if($_!~/\#/){
                $annot=(split/\s+/,$_)[3];
                $annot_trans=(split/\s+/,$_)[5];
                if($cutoff eq "None"){
                        ($chr2,$loc2,$strand2,$dep2,$fail2)=(split/\s+/,$_)[0,1,2,9,12];
                }else{
                        ($chr2,$loc2,$strand2)=(split/\s+/,$_)[0,1,2];
                        @all2s=(split/\s+/,$_)[14,15,16,17,18,19,20,21,22,23,24,25];
                        for $all2(@all2s){
                                ($cut2,$count2)=(split/\;/,$all2)[0,1];
                                if($cut2 == $cutoff){
                                        ($dep2,$fail2)=(split/\,/,$count2)[1,2];
                                        last;
                                }
                        }
                }
                        $line3=join("\t",$chr2,$loc2,$strand2);
                        if(exists $hash{$line3} && $dep2>=$dep_cutoff && $annot !~/NA/){
                                if($fail2==0){
                                $v2=sprintf "%0.5f",$fail2/$dep2;
                                $fdr="NA";
                                }else{
                                $v2=sprintf "%0.5f",$fail2/$dep2;
                                $v1=(split/\s+/,$hash{$line3})[2];
                                $fdr=sprintf "%0.5f",$v1/$v2;
                                }
                                print "$line3\t$annot\t$annot_trans\t$hash{$line3}\t$fail2\t$dep2\t$v2\t$fdr\n";
                        }
                }
        }

close F1;
close F2;

