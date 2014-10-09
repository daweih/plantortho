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
my $split_by_kind;
my $kind_column_no = 2;
my $file_name = $ARGV[0];
my $file_extension;

open I, "< binomial_full_abbreviation.txt";
$opts{species_list} = "binomial_full_abbreviation.txt";
#open I, "< ". $opts{species_list};
my $species;
while(<I>){
	chomp;
	my @r = split /\t/, $_;
	$species->{$r[0]} = $r[$kind_column_no];
}
close I;

=cut
print "#Input\n";
print "@ARGV\n";print "\n";
=cut

if(scalar(@ARGV)>1){
	print "input one file.\n";
	exit;
}

if($ARGV[0] =~ /.*\.([c|t]{1}sv)$/){
#	print "csv or tsv\n";
	$file_extension = $1;
	while(<>){
		chomp;
		my @r = split /\t|,/, $_;
		foreach my $species_abbrev (keys %{$species}){
			my @binomial = split /_/,$species_abbrev;
			if(lc($_) =~ /^.*($binomial[0]\s$binomial[1]).*$/ ){
#				print $species->{$species_abbrev},"\n";
				$split_by_kind->{$species->{$species_abbrev}} .= $_. "\n";
			}
		}
	}
}elsif($ARGV[0] =~ /.*\.(f(a|asta|na|aa))$/){
#	print "fasta\n";
	$file_extension = $1;
	my $fa = Bio::SeqIO->new(-file =>  $ARGV[0],
						 -format => 'Fasta');
	while ( my $seq = $fa->next_seq() ) {
		foreach my $species_abbrev (keys %{$species}){
			if(lc($seq->id) =~ /^.*($species_abbrev).*$/ ){
#				print $species->{$species_abbrev},"\n";
				$split_by_kind->{$species->{$species_abbrev}} .= $seq->id ."\n". $seq->seq. "\n";
			}
		}
	}
}

foreach my $kind (keys %{$split_by_kind}){
	$file_name =~ s/$file_extension//g;
	open O, "> $file_name$kind.$file_extension";
	print O $split_by_kind->{$kind};
	close O;
}
__END__