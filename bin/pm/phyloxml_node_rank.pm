# Contact     : daweimhuang@gmail.com
# Date        : 
# Last Update : Wed Nov  5 13:51:48 CST 2014
# Reference   : 
# Description : 
# Require     :

package phyloxml_node_rank;
use strict;
use warnings;
require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(phyloxml_node_rank);
our $VERSION = '0.01';

#Usage: $result=round($number, $decimalDigits);
#=======round()=similar to C language function==================================================================
sub phyloxml_node_rank{
	my ($node) = @_;

	my $ac = $node->annotation();
	foreach my $annotation ($ac->get_Annotations("taxonomy")){
		my $rank = $annotation->{"_annotation"}->{"rank"}[0]->{"_annotation"}->{"_text"}[0]->{"value"};
		return $rank;
	}
}
1;
__END__
