# Contact     : daweih@me.com
# Date        : 
# Last Update : 
# Reference   : 
# Description : 

#===============================================================================================================
use FindBin qw($Bin); use lib "$Bin/../bin/pm";

use lib "/opt/perl5.12.3/lib/site_perl/5.12.3";
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Min;
use arrayMatching;
use round;

use Bio::TreeIO;

my %opts;
GetOptions(\%opts,'option:s');
my $usage= <<"USAGE";
	Program: $0
	INPUT:
		-option			some option content
	OUTPUT:
	
USAGE

#die $usage unless ($opts{option});
my $startTime=localtime();
################################ Main ##########################################################################

my @id;
my $treeio = Bio::TreeIO->new(-format => 'newick',
								-file   => '../pipeline/phyml_result/01_04tree_al2.phy_phyml_boot100_tree.txt');
while( my $tree = $treeio->next_tree ) {
	my @nodes  = $tree->get_leaf_nodes();
	foreach my $node (@nodes){
		my $id = $node->id;
		$id =~ s/_\w{4}_$//g;
#		print "$id\n";
		push @id, $id;
	}
}
foreach(@id){
	print $_, "\n";
}

__END__
my @id;
my $treeio = Bio::TreeIO->new(-format => 'newick',
								-file   => '../pipeline/OS05G0113900.genetree.newick.tree');
while( my $tree = $treeio->next_tree ) {
	my @nodes  = $tree->get_leaf_nodes();
	foreach my $node (@nodes){
		my $id = $node->id;
		$id =~ s/_\w{4}_$//g;
#		print "$id\n";
		push @id, $id;
	}
}

my @title_row = @id;
my $ortholog_info;
open I, "< ../pipeline/orthologues-ComparaOrthologs-Oryza_sativa-Gene-Compara_Ortholog-76-OS05G0113900.csv";
while(<I>){
	chomp;
	next if(/.*Species.*/);
	my @r = split /,/,$_;
	foreach(@r){$_ =~ s/^"//;$_ =~ s/"$//;}

	my $score = 1;
	my $q;
	foreach(@title_row){
		#Solyc01g099410.2.1_SOLLC
		$_ =~ s/_\w{5}$//;
		my $title_row = $1 if( $_ =~ /^([\w|\d]+)[\.]*.*/ );
		my $query = $title_row;
		my @query = split //,$query;
		my $target = $r[3];
		my @target = split //,$target;
		my $editDistance = &arrayMatching(\@query, \@target) / length($title_row);
#		print $editDistance, "\t";
		if( $score >= $editDistance ){
			$score = $editDistance;
			$q = $query;
		}
	}
	if($score < 0.1){
		$ortholog_info->{$r[3]} = $r[0];
	}
#		print "$score\t$r[0]\t$r[3]\t$q\n";
}
close I;

open I, "< ../pipeline/rio_out.txt";
my $line_no = 0;
my @title_column;
while(<I>){
	chomp;
	$_ =~ s/^\s+//;
	my @r = split /\t/,$_;
	if( $line_no == 0 ){
		@title_row = @r;
#		print scalar(@r), "\n";
	}
	else{
		my $title_column = shift @r;
		foreach(0..$#r){
#			print $title_column, "\t", $title_row[$_], "\t", $r[$_], "\n";
		}
	}
	$line_no+=1;
}
close I;

foreach my $ortholog_info_key (keys %{$ortholog_info} ){
	my $score = 1;
	my $q;
	foreach(@title_row){
		#Solyc01g099410.2.1_SOLLC
		$_ =~ s/_\w{5}$//;
		my $title_row = $1 if( $_ =~ /^([\w|\d]+)[\.]*.*/ );
		my $query = $title_row;
		my @query = split //,$query;
		my $target = $ortholog_info_key;
		my @target = split //,$target;
		my $editDistance = &round( &arrayMatching(\@query, \@target) / length($title_row),3 );
#		print $editDistance, "\t";
		if( $score >= $editDistance ){
			$score = $editDistance;
			$q = $query;
		}
	}
	print "$score\t$ortholog_info->{$ortholog_info_key}\t$ortholog_info_key\t$q\n";
}
close I;


################################ Main ##########################################################################
#===============================================================================================================
__END__
my $options="-option $opts{option}";
my $endTime=localtime();
my $Program = $1 if($0 =~ /(.*)\.pl/);
open  LOG,">>$Program\_ProgramRunning.Log";
print LOG "From \<$startTime\> to \<$endTime\>\tperl $0 $options\n";
close LOG;

__END__

0.2	PGSC0003DMG400024748Histone H2A.1 [Source: PGSC_GENE; acc: PGSC0003DMG400024748]	PGSC0003DMT400023064
0.375	TCM_041890Histone H2A 12 isoform 1	EMT03189
0.5	AMTR_s00168p00048800hypothetical protein	ERN10168
0.5	MTR_2g096610Probable histone H2A.1  [Source: UniProtKB/Swiss-Prot; acc: Q2HU68]	Bo2g002630
0.5	MTR_4g063410Histone H2A  [Source: UniProtKB/TrEMBL; acc: G7JTS2]	TRIUR3_34110
0.5	MTR_8g086640Histone H2A  [Source: UniProtKB/TrEMBL; acc: G7LJG8]	EMT10820
0.5	PRUPE_ppa012844mghypothetical protein	EMT12804
0.5	TCM_021990Histone H2A 12	EMJ07281
0.545454545454545	AMTR_s00057p00215390hypothetical protein	Sb02g030950
0.545454545454545	MTR_2g096570Probable histone H2A.2  [Source: UniProtKB/Swiss-Prot; acc: Q2HU65]	Sb02g004970
0.555555555555556	PRUPE_ppa012888mghypothetical protein	Bra028869
0.5625	MTR_4g063280Histone H2A  [Source: UniProtKB/TrEMBL; acc: G7JTR2]	POPTR_0006s08230
0.583333333333333	MTR_4g071150Probable histone H2A.3  [Source: UniProtKB/Swiss-Prot; acc: Q1S053]	TRIUR3_11525
0.636363636363636	Ensembl identifier & gene name	fgenesh1_pg
