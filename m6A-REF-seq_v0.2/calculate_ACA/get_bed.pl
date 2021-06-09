#! perl
# covert temp.exon file to bed file
# usage: perl get_bed.pl -g gtf_file -j junction_file -t tmp_file > bed_file
# bed_file was 0-based

use Getopt::Long;
GetOptions("g=s" => \$gtf_file,"t=s" => \$tmp_file,"j=s" => \$junction,) or die ("Invalid arguments.\n");
open (F1, "$gtf_file") || die "Cant open file: $!";
open (F2, "$junction") || die "Cant open file: $!";
open (F3, "$tmp_file") || die "Cant open file: $!";

while(<F1>){ ## gtf file of reference
	s/\s+$//;
	next if $_ =~ m/^\#/;
	@temp=split/\t/,$_;
	@as=split/\;/,$temp[8];
	if($as[1] =~ /transcript\_id \"(.*)\"/ || $as[2] =~ /transcript\_id \"(.*)\"/ || $as[3] =~ /transcript\_id \"(.*)\"/){
		$id=$1;
		if($temp[2]=~/exon/){
			for $i($temp[3]..$temp[4]){
				$loc=join("\t",$temp[0],$i,$temp[6]);
				push @{$site{$loc}}, $id; ## hash for each site in exon, so this script takes a lot of RAM, but it runs fastly.
			}
		}
	}
}

while(<F2>){ ## junction sites information
	s/\s+$//;
	@temp=split/\s+/,$_;
	$id=shift @temp;
	push @{$hash{$id}}, @temp;
}

while(<F3>){ # tmp file
	s/\s+$//;
	@all=split/\s+/,$_;
	$status=shift @all;
	### 'general' sites
	if($status=~/general$/){
		$t=join("\t",$all[0],$all[1]-1,$all[2],$all[3],$all[4]);
		print "$t\n";
	#	print "general\t$t\n";
	### nothing to be changed for general sites
	### the following three types has gaps or junction sites
	}else{
		### 'junction' refers to reads with gaps, R1 and R2 have overlapped region
		### As reads had been mapped to ref, the gaps in 'junction' type refer to the skipped introns of fragment
		if($status=~/junction$/){
			@gots=();
			@rels=();
			@juncs=split/\,/,$all[5];
			for $junc(@juncs){
				($l,$r)=(split/\-/,$junc)[0,1];
				push @gots,$l;
				push @gots,$r;
				}
				@gots=sort {$a<=>$b} @gots;
				push @rels,$all[1];
				push @rels,@gots;
				push @rels,$all[2];
		## 'junction2' refers to reads with gaps, R1 and R2 do not have overlapped region
		## this type is complicated, the gaps in reads and skipped introns should be considered
		}elsif($status=~/junction2$/){
			@sites=();
			@gots=();
			@rels=();
			@transp=();
			@trans=();
			$aa=(split/\,/,$all[1])[0];
			$bb=(split/\,/,$all[2])[-1];
			$start=join("\t",$all[0],$aa,$all[4]);
			$end=join("\t",$all[0],$bb,$all[4]);
			if(exists ${$site{$start}}[0] && exists ${$site{$end}}[0]){ ## check if the start and end of each fragment
				@start_ids=@{$site{$start}};
				@end_ids=@{$site{$end}};
				my %hash_a = map{$_=>1} @start_ids;
				my %hash_b = map{$_=>1} @end_ids;
				my @names = grep {$hash_a{$_}} @end_ids;
				if($#names>=0){
				@juncs=split/\,/,$all[5];
				push @sites, $all[1];
				push @sites, $all[2];
				for $name(@names){
				if(exists ${$hash{$name}}[0]){
					$len=0;
					@trans=();
					for $line(@{$hash{$name}}){
						($a,$b)=(split/\,/,$line)[0,1];
						if($sites[1]<=$a && $b<=$sites[2]){
							$len+=$b-$a;
							push @trans,$a;
							push @trans,$b;
						}
					}
					$tran=join("\,",@trans);
					$p=$len ."\=" . $tran;
				}
			}				
			if($#trans >=0){
			push @transp, $p;
			undef %uniq_t;
			@uniq_trans = grep { ++$uniq_t{$_} < 2 } @transp;
			@uniq_trans =sort {$b<=>$a} @uniq_trans;
			$aaa=(split/\=/,$uniq_trans[0])[1];
			@aaas=split/\,/,$aaa;
			push @rels, @aaas;
			}
			for $junc(@juncs){
				($l,$r)=(split/\-/,$junc)[0,1];
				push @gots,$l;
				push @gots,$r;
				}
			push @rels,@gots;
			push @rels,@sites;
			@rels=sort {$a<=>$b} @rels;
			undef %uniq_rels;
			@rels = grep { ++$uniq_rels{$_} < 2 } @rels;
			}
		}
		}elsif($status=~/general2$/){
		## 'general2' refers to reads without gaps, R1 and R2 do not have overlapped region
		## only skipped introns on the ref were calculated
			@sites=();
			@gots=();
			@rels=();
			@transp=();
			@trans=();
			$aa=(split/\,/,$all[1])[0];
			$bb=(split/\,/,$all[2])[-1];
			$start=join("\t",$all[0],$aa,$all[4]);
			$end=join("\t",$all[0],$bb,$all[4]);
			if(exists ${$site{$start}}[0] && exists ${$site{$end}}[0]){ ## check if the start and end of each fragment
				@start_ids=@{$site{$start}};
				@end_ids=@{$site{$end}};
				my %hash_a = map{$_=>1} @start_ids;
				my %hash_b = map{$_=>1} @end_ids;
				my @names = grep {$hash_a{$_}} @end_ids;
			if($#names>=0){
				push @sites, $all[1];
				push @sites, $all[2];
				for $name(@names){
				if(exists ${$hash{$name}}[0]){
					 $len=0;
					 @trans=();
					for $line(@{$hash{$name}}){
					($a,$b)=(split/\,/,$line)[0,1];
					if($sites[1]<=$a && $b<=$sites[2]){
						$len+=$b-$a;
						push @trans,$a;
						push @trans,$b;
						}
					}
					$tran=join("\,",@trans);
					$p=$len ."\=" . $tran;
				}
			}
			if($#trans >=0){
				push @transp, $p;
				undef %uniq_t;
				@uniq_trans = grep { ++$uniq_t{$_} < 2 } @transp;
				@uniq_trans = sort {$b<=>$a} @uniq_trans;
				$aaa=(split/\=/,$uniq_trans[0])[1];
				@aaas=split/\,/,$aaa;
				push @rels, @aaas;
			}
			push @rels, @sites;
			@rels=sort {$a<=>$b} @rels;
		}
		}
		}
		## change the sites to bed formate for 'junction', 'junction2' and 'general2'
		if($#rels>=0){
		@starts=();
		@ends=();
		for $i(0..(($#rels-1)/2)){
			$s=$i*2;
			$e=$i*2+1;
			push @starts, $rels[$s]-1;
			push @ends, $rels[$e];
		}
		$start=join(",", @starts);
		$end=join(",",@ends);
		print "$all[0]\t$start\t$end\t$all[3]\t$all[4]\n";
	#	print "junction\t$all[0]\t$start\t$end\t$all[3]\t$all[4]\n";
		}
	}
}

close F1;
close F2;
close F3;
