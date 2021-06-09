#! perl 
# convert bam file to temp file with start, end and junction information of reads.
# usage: samtools view file.sort.bam | perl sam2bed_with_junction.pl > file.tmp

no warnings 'experimental::smartmatch';
@plus_flags_1=qw/97 145/;
@plus_flags_2=qw/99 147/;
@minus_flags_1=qw/81 161/;
@minus_flags_2=qw/83 163/;

while(<>){
	s/\s+$//;
	($name,$flag,$chr,$loc,$cigar)=(split/\s+/,$_)[0,1,2,3,5];
	$len=0;
	@sites=();
	@matches=();            
	@matches=$cigar=~/(\d+[SMNID])/g;
	if($matches[0]=~/(\d+)S/){
	shift @matches;
	}elsif($matches[-1]=~/(\d+)S/){
	pop @matches;
	}
	## if reads without gaps
	if($cigar !~/N/){
		foreach $match(@matches){
		if($match=~/(\d+)M/){
		$len+=$1;
		}elsif($match=~/(\d+)D/){
		$len+=$1;
		}
		}
		$start=$loc;
		$end=$start+$len-1;
		$line=join("\t",$chr,$start,$end);
		push @{$hash{$name}},$line; ## collect the information of reads without gaps
		push @{$flags{$name}},$flag; ## flag information
	## if reads have gaps
	}elsif($cigar =~/N/){ 
	push @sites, $loc;
	foreach $match(@matches){
		if($match=~/(\d+)M/){
		$len+=$1;
		$mid1=$loc+$len-1;
		}elsif($match=~/(\d+)N/){
		$len+=$1;
		$mid2=$loc+$len;
		$mid=$mid1 . "-" . $mid2; ## start and end of each gap
		push @{$hash2{$name}}, $mid; ## collect the information of gaps
		}elsif($match=~/(\d+)D/){
		$len+=$1;
		}
		}
		$start=$loc;
		$end=$start+$len-1;
		$line=join("\t",$chr,$start,$end);
		push @{$hash{$name}},$line; ## collect the information of reads with gaps
		push @{$flags{$name}},$flag; ## flag information
	}
}

## uniq and sort the gap regions of reads
for $key(keys %hash2){
	undef %uniq;
	@juncs = grep { ++$uniq{$_} < 2 } @{$hash2{$key}};
	if($#juncs==0){
	push @{$junction{$key}},@juncs;
	}elsif($#juncs>0){
	@juncs=sort {$a<=>$b} @juncs;
	push @{$junction{$key}},$juncs[0];
	($pre_a,$pre_b)=(split/\-/,$juncs[0])[0,1];
	for $i(1..$#juncs){
		($a,$b)=(split/\-/,$juncs[$i])[0,1];
		if($a<$pre_b){
			$pre_a=$a;
			$pre_b=$b;
			next;
		}else{
			push @{$junction{$key}},$juncs[$i];
			$pre_a=$a;
			$pre_b=$b;
		}
	}
}
}

## for paired-end mapped reads: 
## 'junction' refers to reads with gaps, R1 and R2 have overlapped region
## 'junction2' refers to reads with gaps, R1 and R2 do not have overlapped region
## 'general' refers to reads without gaps, R1 and R2 have overlapped region
## 'general2' refers to reads without gaps, R1 and R2 do not have overlapped region
for $key(keys %hash){
	@ss=();
	@sort=();
	if($#{$hash{$key}}==1){ ## paired-end mapped reads
		($c1,$s1,$e1)=(split/\s+/,${$hash{$key}}[0])[0,1,2];                
		($c2,$s2,$e2)=(split/\s+/,${$hash{$key}}[1])[0,1,2];
	if($c1 eq $c2){ ## if R1 and R2 mapped to the same chromosome
		push @ss,$s1;
		push @ss,$e1;
		push @ss,$s2;
		push @ss,$e2;
		@flags=@{$flags{$key}};
		@flags=sort {$a<=>$b} @flags;
		if(@flags ~~ @plus_flags_1 || @flags ~~ @plus_flags_2){
			$frag_strand="+";
		}elsif(@flags ~~ @minus_flags_1 || @flags ~~ @minus_flags_2){
			$frag_strand="-";
		}
		if(($s2>=$s1 && $s2<=$e1 && $e1>=$s2 && $e1 <= $e2) || ($s1>=$s2 && $s1<=$e2 && $e2>=$s1 && $e2 <= $e1)){
			## if R1 and R2 have overlapped region
			@sort=sort{$a<=>$b} @ss;
			if(exists ${$junction{$key}}[0]){ ## if paired-end reads have gaps
			$splice_line=join(",",@{$junction{$key}});
			print "junction\t$c1\t$sort[0]\t$sort[-1]\t$key\t$frag_strand\t$splice_line\n"; 
		}else{ ## if paired-end reads do not have gaps
			print "general\t$c1\t$sort[0]\t$sort[-1]\t$key\t$frag_strand\n";
		}
		}elsif($s2>$e1 || $s1>$e2){ ## if R1 and R2 do not have overlapped region
			@sort=sort{$a<=>$b} @ss;
			if(exists ${$junction{$key}}[0]){
			$splice_line=join(",",@{$junction{$key}}); ## if paired-end reads have gaps
			print "junction2\t$c1\t$sort[0]\t$sort[-1]\t$key\t$frag_strand\t$splice_line\n";
		}else{ ## if paired-end reads do not have gaps
			print "general2\t$c1\t$sort[0]\t$sort[-1]\t$key\t$frag_strand\n";
			}
		}
	}
}
}

