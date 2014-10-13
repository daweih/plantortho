#!/usr/bin/perl;
use strict;
use warnings;
###############################################################
# Author: yangl
# Date: 2014.9.29
# Description: This program is to get normal .fa file from newick.parser.pl output. 
###############################################################

my $usage = "\tDescription: get normal .fa file from newick.parser.pl output.\n\tUsage: $0 <fastafile>\t<outfile>\t<output_id_index>\n";
if(!@ARGV){
        print STDERR $usage;
        exit 1;
}

open IN,"<","$ARGV[0]";
open OUT,">","$ARGV[1]";
open OUT2,">","$ARGV[2]";
my $i = 1;
while(<IN>){
	chomp;
	next if /^no/;
	my @entry = split /\s/; 
	print OUT ">$i\n$entry[3]\n";
	print OUT2 "$i\t$entry[2]_$entry[1]_\n";
	$i++;
}
close IN;
close OUT;
