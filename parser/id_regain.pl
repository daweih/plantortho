#!/usr/bin/perl;
use strict;
use warnings;
###############################################################
 # Author: yangl
 # Date: 2014.9.29
 # Description: This program is to regain long ID for phylip format file. 
###############################################################
my $usage = "\tDescription:This program is to regain long ID for phylip format file.\n\tUsage: $0 <shortname_phylip_file>\t<id_index>\t<longname_phylip_file>\n";
if(3!=@ARGV){
        print STDERR $usage;
        exit 1;
}

open IN1,"<","$ARGV[0]";
open IN2,"<","$ARGV[1]";
open OUT1,">","$ARGV[2]";
my %index;
while(<IN2>){
	chomp;
	my @id_index = split /\s+/;
	$index{$id_index[0]} = $id_index[1];

}
=cut
foreach my $key (sort keys %index){
	my $value = $index{$key};
	print "$key => $value\n";
}
=cut
while(<IN1>){
	chomp;
	if (/^(\d+)(\s+)(.*$)/){
		my $i = $1;
#		s/$i/$index{$i}/;
#		print "$_\n";
#		print "$1\n";
#		print "$2\n";
#		print OUT1;
		my $j = $3;
#		print "$2\n";
		printf OUT1 "%-40s$j\n","$index{$i}";
	}else{
		print OUT1 "$_\n";
	}
}

close IN1;
close IN2;
close OUT1;

