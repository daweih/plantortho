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
GetOptions(\%opts,'newick_tree:s','species_list:s', 'pep_dir:s');
my $usage= <<"USAGE";
	Program: $0
	INPUT:
		-newick_tree			
		-species_list
		-pep_dir
	OUTPUT:
	
USAGE

#die $usage unless ($opts{newick_tree});
#die $usage unless ($opts{species_list});
#die $usage unless ($opts{pep_dir});
my $startTime=localtime();
################################ Main ##########################################################################
my $treeio = Bio::TreeIO->new(-format => 'newick',
								-file   => '../sample/ensembl/web_download/OS05G0113900.genetree.newick.tree');$opts{newick_tree} = "../sample/ensembl/web_download/OS05G0113900.genetree.newick.tree";
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
open I, "< binomial_full_abbreviation.v22.txt";
$opts{species_list} = "binomial_full_abbreviation.v22.txt";
#open I, "< ". $opts{species_list};
my $species;
while(<I>){
	chomp;
	my @r = split /\t/, $_;
	$species->{$r[1]} = 1;
}
close I;

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
$opts{pep_dir} = "../sample/ensembl/release-22/plants/fasta/pep/ortho_OS05G0113900/";
opendir PEP , $opts{pep_dir} or die "Cannot open plants_pep: $!";

foreach my $file (readdir PEP)
{
	next if($file !~ /^.*.pep.all.fa$/);
	my $check_file_species = 0;
	foreach my $accession (keys %{$accession_list}){

		next if(!defined $species->{$accession_list->{$accession}->{species_abr}});

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
		my $pep = Bio::SeqIO->new(-file => $opts{pep_dir}. "/$file" ,
										-format => 'Fasta');
		while ( my $seq = $pep->next_seq() ) {
			foreach my $accession (keys %{$accession_list}){
				if($seq->id =~ /^($accession).*$/ ){
					$accession_list->{$accession}->{seq} = $seq->seq;
				}
			}
		}
	}
}
closedir PEP;

foreach my $accession (keys %{$accession_list}){

	next if(!defined $species->{$accession_list->{$accession}->{species_abr}});

	if(defined $accession_list->{$accession}->{seq}){
		print $accession_list->{$accession}->{species_abr}, "\t$accession\t", $accession_list->{$accession}->{seq},"\n";
	}
	else{
	}
}

################################ Main ##########################################################################
#===============================================================================================================
my $options="-newick_tree $opts{newick_tree} -species_list $opts{species_list} -pep_dir $opts{pep_dir} ";
my $endTime=localtime();
my $Program = $1 if($0 =~ /(.*)\.pl/);
open  LOG,">>$Program\_ProgramRunning.Log";
print LOG "From \<$startTime\> to \<$endTime\>\tperl $0 $options\n";
close LOG;
__END__