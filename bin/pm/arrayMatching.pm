# Contact     : daweih@me.com
# Date        :	Wed Jul 25 12:47:37 AST 2012
# Last Update : 
# Reference   : 
# Description : used to get edit distance between two array
# Require     : Min

package arrayMatching;#
use strict;
use warnings;
require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(arrayMatching);#
our $VERSION = '1.0';#
use FindBin qw($Bin); use lib "$Bin/../../bin/pm";
use Min;
#==================================== test data =============================================#
=cut
my $query = 'EIQADEVRL';
my @query = split //,$query;
my $target = 'EIQADEVRL';
my @target = split //,$target;
my $editDistance = &arrayMatching(\@query, \@target);
print $editDistance, "\n";
=cut
#==================================== END test data =========================================#
# USAGE: my $editDistance = &arrayMatching(\@query, \@target);

sub arrayMatching{
	my ($q_ref, $t_ref) = @_;
	
	my $TLEN = scalar(@{$t_ref});
	my $QLEN = scalar(@{$q_ref});
	
	my $distance_matrix_ref = [];
	
	#the first row and first column of distance matrix
	for (my $q=0;$q <= $QLEN ;++$q)
	{
		$distance_matrix_ref->[$q][0] = $q;
	}
	for (my $t=0;$t <= $TLEN ;++$t)
	{
		$distance_matrix_ref->[0][$t] = 0;	
	}	
	
	#compute the edit distances.
	for (my $q=1; $q <= $QLEN ; ++$q)
	{
		for (my $t=1; $t <= $TLEN ; ++$t)
		{
			$distance_matrix_ref->[$q][$t] = &Min
			(
				${$t_ref}[$t-1] eq ${$q_ref}[$q-1]
				? $distance_matrix_ref->[$q-1][$t-1] : $distance_matrix_ref->[$q-1][$t-1] + 1,

				$distance_matrix_ref->[$q][$t-1] + 1,

				$distance_matrix_ref->[$q-1][$t] + 1
			)
		}
	}
=cut
	# Print resulting edit distance array
	for (my $q=0; $q <= $QLEN ; ++$q) {
		for (my $t=0; $t <= $TLEN ; ++$t) {
			print $distance_matrix_ref->[$q][$t], " ";
	}
	print "\n"; }
=cut
	my $editDistance = &Min(@{$distance_matrix_ref->[$QLEN]});
	return	$editDistance;
}
1;
