#! perl -w
# calculate the number of undigest, digest and sum for each site
# perl cal_num.pl -t file_tmp -d file_digest -m ref_motif_exon > file_num

use Getopt::Long;

GetOptions("t=s" => \$tmp,"d=s" => \$digest_file,"m=s" => \$motif_file) or die ("Invalid arguments.\n");
open (F1, "$tmp") || die "Cant open file: $!";
open (F2, "$digest_file") || die "Cant open file: $!";
open (F3, "$motif_file") || die "Cant open file: $!";

while(<F1>){
	s/\s+$//;
	($chr,$s,$e,$strand)=(split/\s+/,$_,)[1,2,3,5];
	if($strand eq "+"){
		$end{$chr . "\t" . $e}++;
	}elsif($strand eq "-"){
		$start{$chr . "\t" . $s}++;
	}
}

while(<F2>){
	s/\s+$//;
	@b=split/\s+/,$_;
	if($b[4] eq "+"){
		$mid{$b[0] ."\t" . $b[1]}=$b[2];
		$right{$b[0] ."\t" . $b[1]}=$b[3];
	}elsif($b[4] eq "-"){
		$mid{$b[0] ."\t" . $b[1]}=$b[2];
		$left{$b[0] ."\t" . $b[1]}=$b[3];
	}
}

while(<F3>){
	s/\s+$//;
	($chr,$loc,$strand)=(split/\s+/,$_)[0,1,2];
	if($strand eq "+"){
		if(exists $right{$chr . "\t" . $loc}){
			$r=$right{$chr . "\t" . $loc};
			$m=$mid{$chr . "\t" . $loc};
		}else{
			$r=0;
			$m=0;
		}
		$presite=$loc-1;
		if(exists $end{$chr . "\t" . $presite}){
			$l=$end{$chr . "\t" . $presite};
		}else{
			$l=0;
		}
	}elsif($strand eq "-"){
		if(exists $left{$chr . "\t" . $loc}){
			$l=$left{$chr . "\t" . $loc};
			$m=$mid{$chr . "\t" . $loc};
		}else{
			$l=0;
			$m=0;
		}
		$presite=$loc+1;
		if(exists $start{$chr . "\t" . $presite}){
			$r=$start{$chr . "\t" . $presite};
		}else{
			$r=0;
		}
	}
#	$mean=($l+$r)/2;
	if($l>=$r){
		$max=$l;
	}else{
		$max=$r;
	}
	$sum=$m+$max;
	if($sum>0){
	$meth=sprintf "%0.4f", $m/$sum;
	print "$chr\t$loc\t$strand\t$m\t$max\t$sum\t$meth\n";
}
}

