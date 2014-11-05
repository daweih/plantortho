#!/usr/bin/perl;
use strict;
use warnings;
###############################################################
# Author: yangl
# Date: 2014.9.29
# Description: This program is to get normal .fa file from newick.parser.pl output. 
###############################################################

my $usage = "\tDescription: get normal .fa file from newick.parser.pl output.\n\tUsage: $0 <fastafile>\t<full_index>\t<outfile>\t<output_id_index>\n";
if(!@ARGV){
        print STDERR $usage;
        exit 1;
}

open IN1,"<","$ARGV[0]";
open IN2,"<","$ARGV[1]";
open OUT,">","$ARGV[2]";
open OUT2,">","$ARGV[3]";

my %hash;
my $i = 1;

while(<IN2>){
	chomp;
	my @fields = split /\t/;
#	print "$fields[3]\n";
	my @entry = split (/_/,$fields[3]);
	my $key =  substr($entry[0],0,1).substr($entry[1],0,3);
#	print $key . "\n";
	$hash{$key}=$fields[4];
}

while(<IN1>){
	chomp;
	next if /^no/;
	my @entry = split /\s/; 
	print OUT ">$i\n$entry[3]\n";
	print OUT2 "$i\t$entry[2]_$hash{$entry[1]}\n";
#	print $hash{$entry[1]} . "\n";
	$i++;
}
close IN1;
close IN2;
close OUT;
close OUT2;
