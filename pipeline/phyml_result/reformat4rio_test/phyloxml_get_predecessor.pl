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

my $tree_species_abrv;
my $treeio = Bio::TreeIO->new(-format => 'newick',
								-file   => '../accuratemulalign3_longname.phy_phyml_boot_trees.txt');
while( my $tree = $treeio->next_tree ) {
	my @nodes  = $tree->get_leaf_nodes();
	foreach my $node (@nodes){
#		print $node->id, "\n";
		my $name_abbriviation;
		$name_abbriviation = $1 if($node->id =~ /.*_(\w{4})_/);
		$tree_species_abrv->{$name_abbriviation} = 1;
	}
	last;
}


open I, "< ../../../parser/binomial_full_abbreviation.txt";
my $binomial;
while(<I>){
	chomp;
	#oryza_sativa	Osat	monocot
	my @r = split /\t/,$_;
	next if( ! defined $tree_species_abrv->{$r[1]} );
	$binomial->{$r[0]}->{binomial_full_abbreviation} = \@r;
}
close I;

$treeio = Bio::TreeIO->new(-format => 'phyloxml',
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

			$node->set_tag_value("using","1");

			$binomial->{$leaf_scientific_name_formatted}->{leaf_node} = $node;

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
=cut
			print &get_scientific_name($node), "\t";			
			print &get_rank($predecessor), "\t";			
			print &get_scientific_name($predecessor), "\n";
=cut
        }
	}
	
}
__END__
my $ncbi_taxonomy = Bio::TreeIO->new(-format => 'phyloxml',
#								-file   => '../../../bin/RIO/bcl_2.xml');
#								-file   => '../../../parser/species_tree_rio.Plants.xml');
								-file   => '../../../bin/RIO/ncbi_taxonomy.xml');
while( my $tree = $ncbi_taxonomy->next_tree ) {
	foreach my $node ($tree->get_nodes){
		if( $node->is_Leaf ) {
			my $leaf_scientific_name_formatted = lc(&get_scientific_name($node));
			   $leaf_scientific_name_formatted =~s/\s/_/;
			my $need_or_not = 0;
			foreach my $binomial_key (keys %{$binomial}){
				if( ! defined $binomial->{$binomial_key}->{leaf_node} ){
					if( $leaf_scientific_name_formatted eq $binomial_key ){
						$need_or_not = 1;
#						my @tags = $node->get_all_tags();						
#						print $node->get_tag_values("taxonomy");
#						print $node->get_tag_values("using");
#						print "found: $binomial_key\n";
						my $ancestor = $node->ancestor;
						while($ancestor){
							print $ancestor->depth, "\n";
#							$ancestor->set_tag_value("using","1");
							$ancestor = $ancestor->ancestor;
						}
					}
				}
			}
			if($need_or_not == 0){
			}
		}
	}
	print  $tree->number_nodes, "\n";	
=cut
	foreach my $node ($tree->get_nodes){
		if($node->has_tag("using")){
		
		}
		else{
#			$tree->remove_Node($node);
		}
	}
=cut
	print  $tree->number_nodes, "\n";
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
