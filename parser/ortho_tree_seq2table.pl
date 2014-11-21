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
#		print "$score\t$r[0]\t$r[3]\t$q\n";
	}
}
close I;

open I, "< ../pipeline/bootstrap100.rio.txt";
my $line_no = 0;
my @title_column;
my $title_row_selected;
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
			$q = $_;
		}
	}
	if( $score < 0.05 ){
#		print "$score\t$ortholog_info->{$ortholog_info_key}\t$ortholog_info_key\t$q\n";
		$title_row_selected->{$q} = 1;
=cut
0.000	Arabidopsis lyrata	fgenesh1_pg.C_scaffold_8001992Histone H2A  [Source: UniProtKB/TrEMBL; acc: D7MT39]	fgenesh1_pg
0.000	Arabidopsis lyrata	fgenesh1_pg.C_scaffold_8001991Histone H2A  [Source: UniProtKB/TrEMBL; acc: D7MT38]	fgenesh1_pg
0.000	Arabidopsis thaliana	AT5G59870HTA6 histone H2A 6 [Source: TAIR_LOCUS; acc: AT5G59870]	AT5G59870
0.000	Arabidopsis thaliana	AT5G02560HTA12 histone H2A 12 [Source: TAIR_LOCUS; acc: AT5G02560]	AT5G02560
0.000	Brassica rapa	Bra002523AT5G59870 (E=6e-053) HTA6 | HTA6; DNA binding 	Bra002523
0.000	Brassica rapa	Bra006688AT5G59870 (E=6e-054) HTA6 | HTA6; DNA binding 	Bra006688
0.000	Brassica rapa	Bra020278AT5G59870 (E=9e-032) HTA6 | HTA6; DNA binding 	Bra020278
0.000	Brassica rapa	Bra028869AT5G02560 (E=2e-039) HTA12 | HTA12; DNA binding 	Bra028869
0.000	Glycine max	GLYMA13G40900Histone H2A  [Source: UniProtKB/TrEMBL; acc: I1M4U3]	GLYMA13G40900
0.000	Glycine max	GLYMA15G04530Histone H2A  [Source: UniProtKB/TrEMBL; acc: I1MDH5]	GLYMA15G04530
0.000	Glycine max	GLYMA13G40890Histone H2A  [Source: UniProtKB/TrEMBL; acc: C6TMV8]	GLYMA13G40890
0.000	Glycine max	GLYMA15G04540Histone H2A  [Source: UniProtKB/TrEMBL; acc: C6SWA6]	GLYMA15G04540
0.000	Glycine max	GLYMA13G40940Histone H2A  [Source: UniProtKB/TrEMBL; acc: C6SZZ2]	GLYMA13G40940
0.000	Glycine max	GLYMA15G04520Histone H2A  [Source: UniProtKB/TrEMBL; acc: C6SWA6]	GLYMA15G04520
0.045	Musa acuminata	GSMUA_Achr10G02580_001Probable histone H2A.1 [Source: GMGC_GENE; acc: GSMUA_Achr10G02580_001]	GSMUA_Achr10P02580_001
0.000	Oryza barthii	OBART05G01020No description	OBART05G01020
0.000	Oryza brachyantha	OB05G11010Histone H2A  [Source: UniProtKB/TrEMBL; acc: J3M3C7]	OB05G11010
0.000	Oryza glaberrima	ORGLA05G0010500Histone H2A  [Source: UniProtKB/TrEMBL; acc: I1PRV7]	ORGLA05G0010500
0.000	Oryza glumaepatula	OGLUM05G01030No description	OGLUM05G01030
0.000	Oryza meridionalis	OMERI05G01050No description	OMERI05G01050
0.000	Oryza nivara	ONIVA05G01040No description	ONIVA05G01040
0.000	Oryza punctata	OPUNC05G00910No description	OPUNC05G00910
0.000	Oryza sativa Indica	BGIOSGA019051Probable histone H2A.6  [Source: UniProtKB/Swiss-Prot; acc: A2XZN0]	BGIOSGA019051
0.000	Populus trichocarpa	POPTR_0006s08230Putative uncharacterized protein [Source: UniProtKB/TrEMBL; acc: A9PDX1_POPTR]	POPTR_0006s08230
0.000	Solanum lycopersicum	Solyc01g099410.2Histone H2A.1  [Source: UniProtKB/Swiss-Prot; acc: P25469]	Solyc01g099410
0.000	Vitis vinifera	VIT_00s0259g00040Histone H2A  [Source: UniProtKB/TrEMBL; acc: F6GZ08]	VIT_00s0259g00040
0.000	Vitis vinifera	VIT_06s0004g04270Histone H2A  [Source: UniProtKB/TrEMBL; acc: A5AFR3]	VIT_06s0004g04270
0.000	Vitis vinifera	VIT_00s0753g00020Histone H2A  [Source: UniProtKB/TrEMBL; acc: F6HWV9]	VIT_00s0753g00020
0.000	Vitis vinifera	VIT_00s0259g00020Histone H2A  [Source: UniProtKB/TrEMBL; acc: F6GZ06]	VIT_00s0259g00020
0.000	Vitis vinifera	VIT_00s0194g00360Histone H2A  [Source: UniProtKB/TrEMBL; acc: F6HCX1]	VIT_00s0194g00360
0.000	Vitis vinifera	VIT_08s0040g03300Histone H2A  [Source: UniProtKB/TrEMBL; acc: A5BE97]	VIT_08s0040g03300
=cut
	}
}
close I;


open I, "< ../pipeline/bootstrap100.rio.txt";

$line_no = 0;
while(<I>){
	chomp;
	$_ =~ s/^\s+//;
	my @r = split /\t/,$_;
	if( $line_no == 0 ){
		@title_row = @r;
#		print scalar(@r), "\n";
		foreach(@title_row){
			$_ =~ s/_\w{5}$//;
		}
	}
	else{
		my $title_column = shift @r;
		my $species_id = $1 if( $title_column =~ /.*_(\w{5})$/ );
		$title_column =~ s/_\w{5}$//;
		next if( !defined $title_row_selected->{$title_column} );
		print $species_id, "\t", $title_column;
		foreach(0..$#r){
			if( defined $title_row_selected->{$title_row[$_]} ){
				print "\t", $r[$_];
			}
		}
		print "\n";
	}
	$line_no+=1;
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
