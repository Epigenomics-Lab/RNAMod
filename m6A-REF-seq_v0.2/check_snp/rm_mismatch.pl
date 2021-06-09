#! perl -w
# remove mismatch sites
# usage: perl rm_mismatch.pl -p file.mis.pileup -n file.num > file.mis.num

use Getopt::Long;

GetOptions("p=s" => \$pileup,"n=s" => \$num_file) or die ("Invalid arguments.\n");
open (F1, "$pileup") || die "Cant open file: $!";
open (F2, "$num_file") || die "Cant open file: $!";

while(<F1>){
	s/\s+$//;
	@temp=split/\s+/,$_;
	if($temp[2] eq "A"){
		$other=$temp[6]+$temp[7]+$temp[8];
	}elsif($temp[2] eq "T"){
		$other=$temp[5]+$temp[7]+$temp[8];
	}elsif($temp[2] eq "C"){
		$other=$temp[5]+$temp[6]+$temp[8];
	}elsif($temp[2] eq "G"){
		$other=$temp[5]+$temp[6]+$temp[7];
	}
	$line=join("\t",$temp[0],$temp[1]);
	$hash{$line}=$other;
}

## As MazF digested around the motif 'ACA', so the mismatch on the three bases were all counted
while(<F2>){
	s/\s+$//;
	$all_other_count=0;
	$mis=0;
	($chr,$loc,$strand,$in,$out,$all)=(split/\s+/,$_)[0,1,2,3,4,5];
	if($strand eq "+"){
		$site1=join("\t",$chr,$loc);
		$site2=join("\t",$chr,$loc+1);
		$site3=join("\t",$chr,$loc+2);
	}elsif($strand eq "-"){
		$site1=join("\t",$chr,$loc);
		$site2=join("\t",$chr,$loc-1);
		$site3=join("\t",$chr,$loc-2);
	}
	if(exists $hash{$site1}){
		$all_other_count+=$hash{$site1};
	}
	if(exists $hash{$site2}){
		$all_other_count+=$hash{$site2};
	}
	if(exists $hash{$site3}){
		$all_other_count+=$hash{$site3};
	}
	if($all_other_count>0){
		$mis=$all_other_count/$all;
		if($in==0){
			$in2=$in;
			$all2=$all;
		}else{
		$in2=$in-$all_other_count;
		$all2=$all-$all_other_count;
		}
		if($all2>0  && $mis <0.2){
			if($in2<0){
				$in2=0;
			}
		$v=sprintf "%0.4f",$in2/$all2;
		print "$chr\t$loc\t$strand\t$in2\t$out\t$all2\t$v\n";
		}
	}else{
	print "$_\n";
	}
}

close F1;
close F2;
