#! perl  
# get junction sites of each transcript
# usage: perl get_junction_sites_trans.pl ref.gtf > junction_file

while(<>){
	s/\s+$//;
	if($_!~/^\#/){
		@temp=split/\t/,$_;
		@as=split/\;/,$temp[8];
		if($as[1] =~ /transcript\_id \"(.*)\"/ || $as[2] =~ /transcript\_id \"(.*)\"/ || $as[3] =~ /transcript\_id \"(.*)\"/){ ## revise for mm10
		$id=$1;
		push @ids, $id;
		if($temp[2]=~/exon/){
			push @{$start{$id}}, $temp[3];
			push @{$end{$id}}, $temp[4];
		}
	}
	}
}

undef %uniq;
@names = grep { ++$uniq{$_} < 2 } @ids;

for $id(@names){
	for $i(0..($#{$start{$id}}-1)){
		@ss=sort {$a<=>$b} @{$start{$id}};
		@es=sort {$a<=>$b} @{$end{$id}};
		$junction=$es[$i] . ",". $ss[$i+1];
		push @{$hash{$id}},$junction;
	}
}

for $key(keys %hash){
	print "$key\t@{$hash{$key}}\n";
}

