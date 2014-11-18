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
GetOptions(\%opts,'newick:s');
my $usage= <<"USAGE";
	Program: $0
	INPUT:
		-newick			some option content
	OUTPUT:
	
USAGE

die $usage unless ($opts{newick});
my $startTime=localtime();
################################ Main ##########################################################################
use Bio::TreeIO;
use Data::Dumper;
my $treeio = Bio::TreeIO->new(-format => 'newick',
								-file   => $opts{newick});
print $treeio->format(), "\n";
while( my $tree = $treeio->next_tree ) {
	my @nodes  = $tree->get_leaf_nodes();
	foreach my $node (@nodes){
		print $node->id;
#		print "\t", $node->branch_length;
		print "\n";
	}
#	last;
}

################################ Main ##########################################################################
#===============================================================================================================
my $options="-newick $opts{newick}";
my $endTime=localtime();
my $Program = $1 if($0 =~ /(.*)\.pl/);
open  LOG,">>$Program\_ProgramRunning.Log";
print LOG "From \<$startTime\> to \<$endTime\>\tperl $0 $options\n";
close LOG;


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
