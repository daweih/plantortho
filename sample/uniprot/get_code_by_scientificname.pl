# Contact     : daweih@me.com
# Date        : Fri Oct 10 15:35:14 CST 2014
# Last Update : 
# Reference   : 
# Description : 
#perl phyloxml_leafnode_scientific_names.pl ../../parser/species_tree_rio.Plants_complete.xml | perl get_code_by_scientificname.pl speclist.txt
#perl phyloxml_leafnode_scientific_names.pl ../../parser/species_tree_rio.Plants29.xml | perl get_code_by_scientificname.pl speclist.txt
#===============================================================================================================
use FindBin qw($Bin); use lib "$Bin/../../bin/pm";

use strict;
use warnings;
use Getopt::Long;

use Data::Dumper;
use phyloxml_node_rank;
use phyloxml_node_scientific_name;

while(<STDIN>){
	chomp;
	my $key = $_;
	my @binomial = split /\s+/,$_;
	my $regexp = join ".*", @binomial;
	my $found = 0;
#	print $regexp, "\n";
	open I, "< $ARGV[0]";
	while(<I>){
		chomp;
		if(/.*N=.*/){
			chomp;
			#AADNV V  648330: N=Aedes albopictus densovirus (isolate Boublik/1994)
			my ($code, $TaxonNode, $scientific_name);
			if($_ =~ /^(\w+)\s*(\w*\s*\d+):\s*N=(.*)$/){
				$code            = $1;
				$TaxonNode       = $2;
				$scientific_name = $3;
				if( $scientific_name =~ /^$regexp$/){
					$found = "$scientific_name\t$code";
				}
			}
		}
	}
	if( $found eq 0 ){
		print "\nNOT FOUND!!\t$key\n\n";
	}
	else{
		print $found, "\n";
	}
}
__END__
#perl phyloxml_leafnode_scientific_names.pl ../../parser/species_tree_rio.Plants29.xml | perl get_code_by_scientificname.pl speclist.txt
output
"
Prunus persica	PRUPE
Glycine max	SOYBN
Medicago truncatula	MEDTR
Populus trichocarpa	POPTR
Brassica rapa subsp. rapa	BRARR
Arabidopsis lyrata	ARALY
Arabidopsis thaliana	ARATH
Vitis vinifera	VITVI
Solanum lycopersicum	SOLLC
Solanum tuberosum	SOLTU
Sorghum bicolor	SORBI
Zea mays	MAIZE
Setaria italica	SETIT
Oryza sativa	ORYSA

NOT FOUND!!	Oryza barthii

Oryza brachyantha	ORYBR
Oryza glaberrima	ORYGL

NOT FOUND!!	Oryza glumaepatula

Oryza sativa subsp. indica	ORYSI

NOT FOUND!!	Oryza meridionalis

Oryza nivara	ORYNI
Oryza punctata	ORYPU
Aegilops tauschii	AEGTA
Hordeum vulgare var. distichum	HORVD
Triticum aestivum	WHEAT
Triticum urartu	TRIUA
Brachypodium distachyon	BRADI
Musa acuminata subsp. malaccensis	MUSAM
Amborella trichopoda	AMBTC
"

#awk '{print $1}' ~/git/plantortho/pipeline/OS05G0113900.ortho.pep.tsv | sort | uniq| wc
output
29
