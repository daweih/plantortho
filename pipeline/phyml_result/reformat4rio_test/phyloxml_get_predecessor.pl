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

open I, "< ../../../parser/binomial_full_abbreviation.v22.txt";
my $binomial;
while(<I>){
	chomp;
	#oryza_sativa	Osat	monocot
	my @r = split /\t/,$_;
	$binomial->{$r[0]} = \@r;
}
close I;

my $treeio = Bio::TreeIO->new(-format => 'phyloxml',
#								-file   => '../../../bin/RIO/bcl_2.xml');
								-file   => '../../../parser/species_tree_rio.Plants.xml');
#								-file   => '../../../bin/RIO/ncbi_taxonomy.xml');
while( my $tree = $treeio->next_tree ) {
#	print Dumper($tree);
#	print Dumper($tree->get_root_node);
	foreach my $node ($tree->get_nodes){
		if( $node->is_Leaf ) {
			my $leaf_scientific_name = &get_scientific_name($node);
			my $leaf_scientific_name_formatted = lc($leaf_scientific_name);
			   $leaf_scientific_name_formatted =~s/species:\s//;
			   $leaf_scientific_name_formatted =~s/\s/_/;
			next if( ! exists $binomial->{$leaf_scientific_name_formatted} );
		
			my $predecessor = $node->ancestor;
			my $ac = $predecessor->annotation();
			while( ! $ac->get_Annotations("taxonomy") ){
				$predecessor = $predecessor->ancestor;
				$ac = $predecessor->annotation();

				while( ! $ac->{"_annotation"}->{"taxonomy"}[0]->{"_annotation"}->{"rank"}[0]->{"_annotation"}->{"_text"}[0]->{"value"} ){
					$predecessor = $predecessor->ancestor;
					$ac = $predecessor->annotation();
				}
#				print $ac->{"_annotation"}->{"taxonomy"}[0]->{"_annotation"}->{"rank"}[0]->{"_annotation"}->{"_text"}[0]->{"value"}, "\n";
			}

			if($ac->get_Annotations("taxonomy")){
#				print Dumper($ac);
				print &get_scientific_name($node), "\t";			
				print &get_rank($predecessor), "\t";			
				print &get_scientific_name($predecessor), "\n";			
			}
			else{
				print "No\n";
			}
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
#			$predecessor = $predecessor->ancestor;
#			$predecessor = $predecessor->ancestor;
#			$predecessor = $predecessor->ancestor;
#			print $name, "\t", $rank;
			return $predecessor;
		}	
	}
}

sub get_scientific_name{
	my ($node) = @_;

	my $ac = $node->annotation();
	foreach my $annotation ($ac->get_Annotations("taxonomy")){
		my $name = $annotation->{"_annotation"}->{"scientific_name"}[0]->{"_annotation"}->{"_text"}[0]->{"value"};
		return $name;
	}
}
sub get_rank{
	my ($node) = @_;

	my $ac = $node->annotation();
	foreach my $annotation ($ac->get_Annotations("taxonomy")){
		my $rank = $annotation->{"_annotation"}->{"rank"}[0]->{"_annotation"}->{"_text"}[0]->{"value"};
		return $rank;
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
