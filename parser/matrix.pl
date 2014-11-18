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
open I, "< ../pipeline/rio_out.txt";
my $line_no = 0;
my @title_row;
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

open I, "< ../pipeline/orthologues-ComparaOrthologs-Oryza_sativa-Gene-Compara_Ortholog-76-OS05G0113900.csv";
while(<I>){
	chomp;
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
	print "$score\t$r[3]\t$q\n";
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
