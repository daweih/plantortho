# Contact     : daweih@me.com
# Date        : Thu Oct  9 15:28:46 CST 2014
# Last Update : 
# Reference   : 
# Description : 根据 monocot/dicot/other 的标签来拆分文件

#===============================================================================================================
#use lib "/opt/perl5.12.3/lib/site_perl/5.12.3";
#!/usr/bin/perl;
use strict;
use warnings;

use Bio::TreeIO;
use Data::Dumper;

my $tree_species_abrv;
my $treeio = Bio::TreeIO->new(-format => 'newick',
								-file   => '../01_04tree_al.phy_phyml_boot_trees.txt');
while( my $tree = $treeio->next_tree ) {
	my @nodes  = $tree->get_leaf_nodes();
	foreach my $node (@nodes){
#		print $node->id, "\n";
		my $name_abbriviation;
		$name_abbriviation = $1 if($node->id =~ /.*_(\w{4})/);
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

foreach my $binomial_key (keys %{$binomial}){
	print $binomial_key, "\n" if( ! defined $binomial->{$binomial_key}->{leaf_node} );
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
