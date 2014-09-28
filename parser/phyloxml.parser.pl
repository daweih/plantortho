use Bio::TreeIO;
use Data::Dumper;
my $treeio = Bio::TreeIO->new(-format => 'phyloxml',
#								-file => 'apaf.xml');
#								-file => 'bcl_2.xml');
								-file => 'ncbi_taxonomy_metazoa.xml');
#								-file => 'species_tree_rio.xml');
#								-file => 'ncbi_taxonomy.xml');
while( my $tree = $treeio->next_tree ) {
	my @nodes  = $tree->get_leaf_nodes();
	foreach my $node (@nodes){
#		print $node->id, "\t";
#		print $node->branch_length, "\t";
#		print "\n";
#		if($node->id eq "Apaf-1_HUMAN"){
#		if($node->id eq "145_XENLA"){
		
		my          $ac = $node->annotation();
		foreach $key ( $ac->get_all_annotation_keys() ) {
              my @values = $ac->get_Annotations($key);
              foreach my $value ( @values ) {
                 # value is an Bio::AnnotationI, and defines a "as_text" method
                 print $key,": ",$value->as_text,"\n";

              }
          }
#		}
	}
}
__END__
print Dumper($treeio);
$VAR1 = bless( {
                 '_handler' => bless( {
                                        '_treelevel' => 0,
                                        'nodetype' => 'Bio::Tree::AnnotatableNode',
                                        '_root_verbose' => 0,
                                        'treetype' => 'Bio::Tree::Tree'
                                      }, 'Bio::TreeIO::TreeEventBuilder' ),
                 '_file' => 'species_tree_rio.xml',
                 '_end_elements' => {
                                      'clade' => sub { "DUMMY" },
                                      'clade_relation' => sub { "DUMMY" },
                                      'phylogeny' => sub { "DUMMY" },
                                      'sequence_relation' => $VAR1->{'_end_elements'}{'clade_relation'}
                                    },
                 '_params' => {},
                 '_root_cleanup_methods' => [
                                              sub { "DUMMY" }
                                            ],
                 '_lastitem' => {},
                 '_flush_on_write' => 1,
                 '_filehandle' => \*Symbol::GEN0,
                 '_mode' => 'r',
                 'nodetype' => 'Bio::Tree::AnnotatableNode',
                 '_reader' => bless( do{\(my $o = '140590063812608')}, 'XML::LibXML::Reader' ),
                 '_root_verbose' => 0,
                 '_start_elements' => {
                                        'clade' => sub { "DUMMY" },
                                        'clade_relation' => sub { "DUMMY" },
                                        'phylogeny' => sub { "DUMMY" },
                                        'sequence_relation' => $VAR1->{'_start_elements'}{'clade_relation'}
                                      },
                 'treetype' => 'Bio::Tree::Tree'
               }, 'Bio::TreeIO::phyloxml' );

				

