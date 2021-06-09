#! perl -w
# change the fasta sequence to two-line mode
# usage: perl delete.pl seq.fa > seq_2l.fa

use strict;
use warnings;
open (IN,$ARGV[0]);
my @infile = <IN>;
for (my $i=0;$i<@infile;$i++) {
        if ($infile[$i]=~m/>/) {
                if ($i!=0) {print "\n";}
                print "$infile[$i]";
        } else {
                chomp $infile[$i];
                print "$infile[$i]";
        }
}
print "\n";

close IN;

