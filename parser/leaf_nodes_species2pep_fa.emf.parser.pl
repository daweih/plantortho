# Contact     : daweih@me.com
# Date        : Sep 17, 2014
# Last Update : 
# Reference   : 
# Description : 

#===============================================================================================================
use lib "/opt/perl5.12.3/lib/site_perl/5.12.3";
use strict;
use warnings;
use Getopt::Long;
use Bio::TreeIO;
use Data::Dumper;
use Bio::SeqIO;

my %opts;
GetOptions(\%opts,'newick_tree:s','species_list:s');
my $usage= <<"USAGE";
	Program: $0
	INPUT:
		-newick_tree			
		-species_list
	OUTPUT:
	
USAGE

#die $usage unless ($opts{newick_tree});
#die $usage unless ($opts{species_list});
my $startTime=localtime();
################################ Main ##########################################################################
my $treeio = Bio::TreeIO->new(-format => 'newick',
								-file   => '../sample/ensembl/web_download/OS05G0113900.genetree.newick.tree');
#								-file   => $opts{newick_tree});
my $accession_list;
while( my $tree = $treeio->next_tree ) {
	my @nodes  = $tree->get_leaf_nodes();
	foreach my $node (@nodes){
#		print $node->id, "\t";
		#GLYMA15G04540.1_Gmax_
		my $id = $node->id;
		if($id =~ /^(\S+)_([A-Z]\S+)_$/){
			my $accession = $1;
			$accession_list->{$accession}->{checktree} = 1;
			$accession_list->{$accession}->{species_abr} = $2;
		}
		
	}
}

open I, "< binomial_full_abbreviation.oryza.txt";
#open I, "< ". $opts{species_list};
my $species;
while(<I>){
	chomp;
	my @r = split /\t/, $_;
	$species->{$r[1]} = 1;
}
close I;
foreach my $accession (keys %{$accession_list}){
	next if(!defined $species->{$accession_list->{$accession}->{species_abr}});
	print $accession, "\t", $accession_list->{$accession}->{species_abr}, "\n";
}
__END__
=cut
my $ortholog_in = Bio::SeqIO->new(-file => "OS05G0113900.ortholog.fasta" ,
								-format => 'Fasta');
while ( my $seq = $ortholog_in->next_seq() ) {
	foreach my $accession (keys %{$accession_list}){
		if($seq->id =~ /^($accession).*$/ ){
			print $1, "\t", $seq->seq, "\n";
			$accession_list->{$accession}->{seq} = $seq->seq;
		}
	}
}
=cut
opendir PEP , "./plants_pep" or die "Cannot open plants_pep: $!";
foreach my $file (readdir PEP)
{
	next if($file !~ /^.*.pep.all.fa$/);
	my $check_file_species = 0;
	foreach my $accession (keys %{$accession_list}){
		my @species_abr = split //,$accession_list->{$accession}->{species_abr};
		#Aegilops_tauschii.GCA_000347335.1.22.pep.all.fa
		#Atau
		if($file =~ /^($species_abr[0]\S+_$species_abr[1]$species_abr[2]$species_abr[3])\S+.*$/ ){
			$accession_list->{$accession}->{pep_fasta} = $file;
			$check_file_species = 1;
		}
	}
	if($check_file_species){
#		print $file, "\n";
		my $pep = Bio::SeqIO->new(-file => "plants_pep/$file" ,
										-format => 'Fasta');
		while ( my $seq = $pep->next_seq() ) {
			foreach my $accession (keys %{$accession_list}){
				if($seq->id =~ /^($accession).*$/ ){
#					print $1, "\t", $seq->seq, "\n";
					$accession_list->{$accession}->{seq} = $seq->seq;
				}
			}
		}
	}
}
closedir PEP;

foreach my $accession (keys %{$accession_list}){
	if(defined $accession_list->{$accession}->{seq}){
		print "yes\t", $accession_list->{$accession}->{species_abr}, "\t$accession\t", $accession_list->{$accession}->{seq},"\n";
	}
	else{
		print "no\t", $accession_list->{$accession}->{species_abr}, "\t$accession\n";
	}
}
################################ Main ##########################################################################
#===============================================================================================================
__END__
my $options="-option $opts{option}";
my $endTime=localtime();
my $Program = $1 if($0 =~ /(.*)\.pl/);
open  LOG,">>$Program\_ProgramRunning.Log";
print LOG "From \<$startTime\> to \<$endTime\>\tperl $0 $options\n";
close LOG;
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
"Species
Type
dN/dS
Ensembl identifier & gene name
Compare
Location
Target %id
Query %id"

Amborella trichopoda
Many-to-many
n/a
AMTR_s00057p00215390hypothetical protein
Region ComparisonAlignment (protein)Alignment (cDNA)Gene Tree (image)
AmTr_v1.0_scaffold00057:4261624-4262748:1
82
78
