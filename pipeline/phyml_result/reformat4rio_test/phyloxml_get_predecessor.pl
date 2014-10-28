# Contact     : daweih@me.com
# Date        : Thu Oct  9 15:28:46 CST 2014
# Last Update : 
# Reference   : 
# Description : 根据 monocot/dicot/other 的标签来拆分文件

#===============================================================================================================
use lib "/opt/perl5.12.3/lib/site_perl/5.12.3";
use strict;
use warnings;

use Bio::TreeIO;
use Data::Dumper;

my $treeio = Bio::TreeIO->new(-format => 'phyloxml',
#								-file   => '../../../bin/RIO/bcl_2.xml');
								-file   => '../../../parser/species_tree_rio.Plants.xml');
#								-file   => '../../../bin/RIO/ncbi_taxonomy.xml');
while( my $tree = $treeio->next_tree ) {
#	print Dumper($tree);
#	print Dumper($tree->get_root_node);
	foreach my $node ($tree->get_nodes){
		if( $node->is_Leaf ) {
			my $predecessor = &predecessor($node, "Oryza sativa");
			print &get_scientific_name($predecessor), "\n";
			last;
#                   print "node is a leaf ... \n";
#			
			

        }
        else{
#        	print "no\n";
        }
	}
	
}

sub predecessor{
	my ($node,$search_string) = @_;

	my $ac = $node->annotation();
	foreach my $annotation ($ac->get_Annotations("taxonomy")){
		my $name = $annotation->{"_annotation"}->{"scientific_name"}[0]->{"_annotation"}->{"_text"}[0]->{"value"};
		my $rank = $annotation->{"_annotation"}->{"rank"}[0]->{"_annotation"}->{"_text"}[0]->{"value"};
		if($name = $search_string){
			my $predecessor = $node->ancestor;
			return $predecessor;
		}	
	}
}

sub get_scientific_name{
	my ($node) = @_;

	my $ac = $node->annotation();
	foreach my $annotation ($ac->get_Annotations("taxonomy")){
		my $name = $annotation->{"_annotation"}->{"scientific_name"}[0]->{"_annotation"}->{"_text"}[0]->{"value"};
		my $rank = $annotation->{"_annotation"}->{"rank"}[0]->{"_annotation"}->{"_text"}[0]->{"value"};
		return $type.": ".$name;
	}
}

__END__
rank:
	superkingdom
	kingdom
	phylum
	superclass
	class
	subclass
	superorder
	order
	suborder
	infraorder
	superfamily
	family
	subfamily
	tribe
	species

Atau
Brap
Gmax
Obar
Obra
Ogla
Oglu
Oind
Omer
Oniv
Opun
Taes
Tura
