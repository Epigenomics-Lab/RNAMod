#! perl -w
# keep the site within exon
# usage: perl get_ref_motif.pl -g ref.gtf -f ref_2l.fa > motif_exon_file
use Getopt::Long;

GetOptions("g=s" => \$ref, "f=s" => \$file) or die ("Invalid arguments.\n");
open (F1, "$ref") || die "Cant open file: $!";
open (F2, "$file") || die "Cant open file: $!";
$motif='ACA';

while(<F1>){
	s/\s+$//;
	next if $_=~/\#/;
	@temp=split/\t/,$_;
	if($temp[2]=~/exon/){
		for $i($temp[3]..$temp[4]){
		$hash{$temp[0] . "\t" . $i . "\t" . $temp[6]}=1;        			
		}
	}
}

while($id=<F2>,$seq=<F2>){
	$id=~s/\s+$//;
	$id=~s/\>//;
	$id_short=(split/\s+/,$id)[0];
	$seq=~s/\s+$//;
	$seq=~tr/atcg/ATCG/;
	$offset=0;
	@nums=();
	$wz=index($seq,$motif,$offset);
	while($wz != -1){
		push @nums, $wz;
		$offset=$wz+1;
		$wz=index($seq,$motif,$offset);
	}
	foreach $num(@nums){
		$start=$num+1;
		if(exists $hash{$id_short . "\t" . $start . "\t+"}){
			print "$id_short\t$start\t\+\n";
		}
	}
	$re_motif=reverse($motif);
	$re_motif=~tr/ATCG/TAGC/;
	$offset2=0;
	@num2s=();
	$wz2=index($seq,$re_motif,$offset2);
		while($wz2 != -1){
			push @num2s, $wz2;
			$offset2=$wz2+1;
			$wz2=index($seq,$re_motif,$offset2);
		}
		foreach $num2(@num2s){
			$start=$num2+3;
			if(exists $hash{$id_short . "\t" . $start . "\t-"}){
				print "$id_short\t$start\t\-\n";
			}
		}
}

close F1;
close F2;

