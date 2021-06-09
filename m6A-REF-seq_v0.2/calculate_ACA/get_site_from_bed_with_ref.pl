#! perl -w
# get the location of ACA motif based on the bed file and genome reference 
# perl get_site_from_bed_with_ref.pl -f ref_2l.fa -m ref_motif_exon -b file.bed > file_digest

use Getopt::Long;

GetOptions("f=s" => \$ref,"m=s" => \$motif_file,"b=s" => \$bed_file) or die ("Invalid arguments.\n");
open (F1, "$ref") || die "Cant open file: $!";
open (F2, "$motif_file") || die "Cant open file: $!";
open (F3, "$bed_file") || die "Cant open file: $!";
$motif='ACA';

while($id=<F1>,$seq=<F1>){ ## genome fasta with two-line mode
        $id=~s/\s+$//;
        $id=~s/\>//;
        $id_short=(split/\s+/,$id)[0];
        $seq=~s/\s+$//;
        $seq=~tr/atcg/ATCG/;
        $genome{$id_short}=$seq;
}

while(<F2>){ ## the location of ACA motif on the exon region of reference genome
        s/\s+$//;
        @tmp1=split/\s+/,$_;
        $line=join("\t",@tmp1);
        $hash{$line}=1;
}

while(<F3>){ ## bed file
	s/\s+$//;
	$frag="";
	@lens=();
	$all_len=0;
#	($chr,$loc1,$loc2,$strand)=(split/\s+/,$_)[1,2,3,5];
	($chr,$loc1,$loc2,$strand)=(split/\s+/,$_)[0,1,2,4];
	@loc1s=split/\,/,$loc1;
	@loc2s=split/\,/,$loc2;
	for $i(0..$#loc1s){
		$len=$loc2s[$i]-$loc1s[$i];
		$all_len+=$len;
		push @lens,$all_len;
		$frag.=substr($genome{$chr},$loc1s[$i],$len);
		}
		$frag_len=length($frag); ## obtain fragment sequence and length
		if($all_len<=2000){ ## remove too long fragment
			if($strand eq "+"){
        $offset=0;
        @nums=();
        $wz=index($frag,$motif,$offset); ## locate motif
        while($wz != -1){ ## get the location of every ACA motif on a fragment
                push @nums, $wz;
                $offset=$wz+1;
                $wz=index($frag,$motif,$offset);
        }
        foreach $num(@nums){ ## deal with each ACA motif, convert the location of the fragment to the real location on the ref genome
        	for $i(0..$#lens){
        	if($num < $lens[$i]){
        		if($i==0){
        			$in=$num+1;
                }elsif($i>0){
                        $in=$num+1-$lens[$i-1];
                }
                $start=$loc1s[$i]+$in;
                $site_plus=$chr . "\t" .$start;
                if(exists $hash{$site_plus . "\t+"}){ ## check whether the ACA motif located on the exon region of ref genome
                push @plus,$site_plus; ## plus means the plus strand
                if($num+1>5){
                        $inner_plus{$site_plus}++; ## count the ACA motif located at the innner of a fragment, which means Mazf undigest this site
                }elsif($num==0){
                        $outer_plus{$site_plus}++; ## count the ACA motif located at the bounder of a fragment, which means Mazf digest this site
                }
                }
                last;
        	}
        }
      }
    }elsif($strand eq "-"){
        $re_motif=reverse($motif); 
        $re_motif=~tr/ATCG/TAGC/; ### the reverse compliment of ACA, the TGT motif, the modification site is the second T
        $offset2=0; ## similar with above code for ACA
        @num2s=();
        $wz2=index($frag,$re_motif,$offset2);
        while($wz2 != -1){
                push @num2s, $wz2;
                $offset2=$wz2+1;
                $wz2=index($frag,$re_motif,$offset2);
        }
        foreach $num2(@num2s){
        	for $i(0..$#lens){
        	if($num2 < $lens[$i]){
        		if($i==0){
        			$in2=$num2+3;
        			}elsif($i>0){
                	$in2=$num2+3-$lens[$i-1];
                	}
                $start=$loc1s[$i]+$in2;
                $site_minus=$chr . "\t" . $start;
                if(exists $hash{$site_minus . "\t-"}){
                push @minus,$site_minus; ## minus means the minus strand
                if($frag_len-($num2+3)>5){
                        $inner_minus{$site_minus}++;
                }elsif($frag_len-($num2+3)==0){ 
                        $outer_minus{$site_minus}++;
                }
                }
                last;
        }
      }
    }
  }
}
}

undef %count1;
undef %count2;
@uniq_plus = grep { ++$count1{ $_ } < 2; } @plus; ## remove duplicates
@uniq_minus = grep { ++$count2{ $_ } < 2; } @minus; ## remove duplicates
for $site(@uniq_plus){ ## report the results
        if(exists $inner_plus{$site} && exists $outer_plus{$site}){
                print "$site\t$inner_plus{$site}\t$outer_plus{$site}\t+\n";
        }elsif(exists $inner_plus{$site} && !exists $outer_plus{$site}){
                print "$site\t$inner_plus{$site}\t0\t+\n";
        }elsif(!exists $inner_plus{$site} && exists $outer_plus{$site}){
                print "$site\t0\t$outer_plus{$site}\t+\n";
        }
}

for $site(@uniq_minus){ ## report the results
        if(exists $inner_minus{$site} && exists $outer_minus{$site}){
                print "$site\t$inner_minus{$site}\t$outer_minus{$site}\t-\n";
        }elsif(exists $inner_minus{$site} && !exists $outer_minus{$site}){
                print "$site\t$inner_minus{$site}\t0\t-\n";
        }elsif(!exists $inner_minus{$site} && exists $outer_minus{$site}){
                print "$site\t0\t$outer_minus{$site}\t-\n";
        }
}

close F1;
close F2;
close F3;
