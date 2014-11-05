# Contact     : daweimhuang@gmail.com
# Date        : Wed Nov  5 13:51:48 CST 2014
# Last Update : 
# Reference   : 
# Description : 
# Require     :

package phyloxml_node_scientific_name;
use strict;
use warnings;
require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(phyloxml_node_scientific_name);
our $VERSION = '0.01';

#Usage: $result=round($number, $decimalDigits);
#=======round()=similar to C language function==================================================================
sub phyloxml_node_scientific_name{
	my ($node) = @_;

	my $ac = $node->annotation();
	foreach my $annotation ($ac->get_Annotations("taxonomy")){
		my $name = $annotation->{"_annotation"}->{"scientific_name"}[0]->{"_annotation"}->{"_text"}[0]->{"value"};
		return $name;
	}
}
1;
__END__
