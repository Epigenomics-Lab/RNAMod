#! perl -w
# calculate the number of A T C G for each sites

while(<>){
	s/\s+$//;
#	next if $_=~/\mpileup/;
	$countA=0;
	$countT=0;
	$countC=0;
	$countG=0;
	$rel=0;
	@chars=();
	@names=();
	undef %hash;
	($chrom,$loc,$ref,$dep,$bases,$line)=(split/\s+/,$_)[0,1,2,3,4,6];
	$bases=~tr/atcg/ATCG/;
	@nums = $bases =~ m/[\-\+][0-9]+[ACGTNacgtn]+/g;
	for $num(@nums){
		if($num=~/^\+/){
		$num=~/([0-9]+)/;
		$bases=~ s/\+[0-9]+[ACGTNacgtn]{$1}//;
		}elsif($num=~/^\-/){
                $num=~/([0-9]+)/;
                $bases=~ s/\-[0-9]+[ACGTNacgtn]{$1}//;
                }
		}
	$bases=~ s/\^[[:ascii:]]//g;
	$bases=~ s/\$//g;
	@chars = $bases =~ m/[ATNGCatncg.,\*\>\<]/g;
	@names=split/\,/,$line;        
	if($ref eq "a" || $ref eq "A"){
		for $i(0..$#names){
			$hash{$names[$i]}=$chars[$i];
		}
		for $key(keys %hash){
			if($hash{$key} eq "," || $hash{$key} eq "."){
				$countA++;
				$rel++;
			}elsif($hash{$key} eq "T"){
				$countT++;
				$rel++;
			}elsif($hash{$key} eq "C"){
				$countC++;
				$rel++;
			}elsif($hash{$key} eq "G"){
				$countG++;
				$rel++;
			}
		}
		if($countT > 0 || $countC > 0 || $countG > 0){
			$other=$countT+$countC+$countG;
			$other_rate=sprintf "%0.4f", $other/$rel;
			print "$chrom\t$loc\t$ref\t$dep\t$rel\t$countA\t$countT\t$countC\t$countG\t$other_rate\n";
		}
	}
	if($ref eq "t" || $ref eq "T"){
		for $i(0..$#names){
			$hash{$names[$i]}=$chars[$i];
		}
		for $key(keys %hash){
			if($hash{$key} eq "A"){
				$countA++;
				$rel++;
			}elsif($hash{$key} eq "," || $hash{$key} eq "."){
				$countT++;
				$rel++;
			}elsif($hash{$key} eq "C"){
				$countC++;
				$rel++;
			}elsif($hash{$key} eq "G"){
				$countG++;
				$rel++;
			}
		}
		if($countA > 0 || $countC > 0 || $countG > 0){
			$other=$countA+$countC+$countG;
			$other_rate=sprintf "%0.4f", $other/$rel;
			print "$chrom\t$loc\t$ref\t$dep\t$rel\t$countA\t$countT\t$countC\t$countG\t$other_rate\n";
		}
	}
	if($ref eq "c" || $ref eq "C"){
		for $i(0..$#names){
			$hash{$names[$i]}=$chars[$i];
		}
		for $key(keys %hash){
			if($hash{$key} eq "A"){
				$countA++;
				$rel++;
			}elsif($hash{$key} eq "T"){
				$countT++;
				$rel++;
			}elsif($hash{$key} eq "," || $hash{$key} eq "."){
				$countC++;
				$rel++;
			}elsif($hash{$key} eq "G"){
				$countG++;
				$rel++;
			}
		}
		if($countA > 0 || $countT > 0 || $countG > 0){
			$other=$countA+$countT+$countG;
			$other_rate=sprintf "%0.4f", $other/$rel;
			print "$chrom\t$loc\t$ref\t$dep\t$rel\t$countA\t$countT\t$countC\t$countG\t$other_rate\n";
		}
	}
	if($ref eq "g" || $ref eq "G"){
		for $i(0..$#names){
			$hash{$names[$i]}=$chars[$i];
		}
		for $key(keys %hash){
			if($hash{$key} eq "A"){
				$countA++;
				$rel++;
			}elsif($hash{$key} eq "T"){
				$countT++;
				$rel++;
			}elsif($hash{$key} eq "C"){
				$countC++;
				$rel++;
			}elsif($hash{$key} eq "," || $hash{$key} eq "."){
				$countG++;
				$rel++;
			}
		}
		if($countA > 0 || $countT > 0 || $countC > 0){
			$other=$countA+$countT+$countC;
			$other_rate=sprintf "%0.4f", $other/$rel;
			print "$chrom\t$loc\t$ref\t$dep\t$rel\t$countA\t$countT\t$countC\t$countG\t$other_rate\n";
		}
	}
}

