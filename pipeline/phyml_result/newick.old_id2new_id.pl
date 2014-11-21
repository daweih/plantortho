# Contact     : daweih@me.com
# Date        : 
# Last Update : 
# Reference   : 
# Description : 

#===============================================================================================================
use lib "/opt/perl5.12.3/lib/site_perl/5.12.3";
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;


use Bio::TreeIO;
use Data::Dumper;
my $key2barcode;
my $treeio = Bio::TreeIO->new(-format => 'newick',
								-file   => "01_04tree_al3.phy_phyml_tree.txt");
while( my $tree = $treeio->next_tree ) {
	my @nodes  = $tree->get_leaf_nodes();
	foreach my $node (@nodes){
		if( $node->id =~ /^(.*)\_(\w{5})$/ ){
			$key2barcode->{$1} = $2;
		}
	}
}


$treeio = Bio::TreeIO->new(-format => 'newick',
								-file   => "01_04tree_al2.phy_phyml_boot100_tree.txt");
my $rexp;
while( my $tree = $treeio->next_tree ) {
	my @nodes  = $tree->get_leaf_nodes();
	foreach my $node (@nodes){
		if( $node->id =~ /^(.*)\_(\w{5})$/ ){
			my $key = $1;
			my $barcode = $2;
			print $node->id, "\t", $key2barcode->{$1}, "\n";
			$rexp->{$node->id} = $key. "_". $key2barcode->{$1};
		}
	}
}
open I, "< 01_04tree_al2.phy_phyml_boot100_trees.txt";
while(<I>){
	foreach my $node_id (keys %{$rexp} ){
		$_ =~ s/$node_id/$rexp->{$node_id}/;
	}
	print $_;
}
close I;
__END__

 java -Xmx2048m -cp ../bin/RIO/forester_1036.jar org.forester.application.rio ./phyml_result/01_04tree_al2.phy_phyml_boot100_tree.txt species_tree_rio.Plants_complete.xml ./phyml_result/01_04tree_al2.phy_phyml_boot100_tree.rio.txt