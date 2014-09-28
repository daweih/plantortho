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


my %opts;
GetOptions(\%opts,'option:s');
my $usage= <<"USAGE";
	Program: $0
	INPUT:
		-option			some option content
	OUTPUT:
	
USAGE

die $usage unless ($opts{option});
my $startTime=localtime();
################################ Main ##########################################################################
################################ Main ##########################################################################
#===============================================================================================================
my $options="-option $opts{option}";
my $endTime=localtime();
my $Program = $1 if($0 =~ /(.*)\.pl/);
open  LOG,">>$Program\_ProgramRunning.Log";
print LOG "From \<$startTime\> to \<$endTime\>\tperl $0 $options\n";
close LOG;


use Bio::TreeIO;
use Data::Dumper;
my $treeio = Bio::TreeIO->new(-format => 'newick',
								-file   => 'gene_trees_rio.nh');
print $treeio->format(), "\n";
while( my $tree = $treeio->next_tree ) {
	my @nodes  = $tree->get_leaf_nodes();
	foreach my $node (@nodes){
		print $node->id, "\t";
		print $node->branch_length, "\t";
		print "\n";
	}
}
__END__

__END__
print Dumper($treeio);
$VAR1 = bless( {
                 '_handler' => bless( {
                                        '_treelevel' => 0,
                                        'nodetype' => 'Bio::Tree::Node',
                                        '_root_verbose' => 0,
                                        'treetype' => 'Bio::Tree::Tree'
                                      }, 'Bio::TreeIO::TreeEventBuilder' ),
                 '_file' => 'gene_trees_rio.nh',
                 '_params' => {
                                'no_bootstrap_values' => 0,
                                'no_branch_lengths' => 0,
                                'no_internal_node_labels' => 0,
                                'bootstrap_style' => 'traditional',
                                'order_by' => '',
                                'newline_each_node' => 0,
                                'internal_node_id' => 'id'
                              },
                 '_root_cleanup_methods' => [
                                              sub { "DUMMY" }
                                            ],
                 '_flush_on_write' => 1,
                 '_filehandle' => \*Symbol::GEN0,
                 '_root_verbose' => 0,
                 '_print_tree_count' => 0
               }, 'Bio::TreeIO::newick' );
