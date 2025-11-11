#!/usr/bin/perl

#  replace_tip_labels.pl
#
#  Created by Mahwash Jamy on 2018-06-19.
#  

use strict;
use warnings;

die "\nReplaces tip labels in a tree file in newick format based on a user specified table.\n\nUsage: replace_tip_labels.pl <tre file> <replacement tsv> <output>\n\n" unless @ARGV == 3;

my ($tree, $tsv, $output) = @ARGV;


my $find = "";
my $replace = "";
my %new_tip = ();
my $new_newick = "";

# Build a hash where the keys (old tip label) point at the new tip label
open(my $in_tsv, "<$tsv") or die "error opening $tsv for reading";

while (my $line = <$in_tsv>) {
    chomp $line;
    ($find, $replace) = split("\t", $line);
    $new_tip{$find} = $replace;
}

close $in_tsv;


# Replace substrings in newick file
open (my $in_tree, "<$tree") or die "error opening $tree for reading";
open (my $out_tree, ">$output") or die "error opening $output for writing";

while (my $line = <$in_tree>) {
     ($new_newick = $line) =~ s/(@{[join "|", keys %new_tip]})/$new_tip{$1}/g;  
}

print $out_tree "$new_newick\n";

close $in_tree;
close $out_tree;