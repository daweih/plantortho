# Contact     : daweih@me.com
# Date        : Wed Nov  5 13:58:39 CST 2014
# Last Update : 
# Reference   : 
# Description : 根据 monocot/dicot/other 的标签来拆分文件
#perl phyloxml_leafnode_scientific_names.pl ../../parser/species_tree_rio.Plants29.xml 
#perl phyloxml_leafnode_scientific_names.pl ../../parser/species_tree_rio.Plants_complete.xml

#===============================================================================================================
use FindBin qw($Bin); use lib "$Bin/../../bin/pm";

use strict;
use warnings;
use Getopt::Long;

use Data::Dumper;
use Bio::TreeIO;
use phyloxml_node_scientific_name;
use phyloxml_node_rank;

my $treeio = Bio::TreeIO->new(-format => 'phyloxml',
								-file   => $ARGV[0]);
while( my $tree = $treeio->next_tree ) {
	foreach my $node ($tree->get_nodes){
		if( $node->is_Leaf ) {
			print &phyloxml_node_scientific_name($node), "\n";			
        }
	}
	
}
__END__