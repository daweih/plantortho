# Contact     : daweimhuang@gmail.com
# Date        : 
# Last Update : Sat Nov 19 10:28:34 AST 2011 by daweimhuang@gmail.com: changeFormat2pm
# Reference   : 
# Description : 

package Uniq;
use strict;
use warnings;
require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(Uniq);
our $VERSION = '2.0';

#对数组取uniq
#Useage: $result=&Uniq($string);
#version: 1
#update : ?
#===============================================================================================================
sub Uniq{
	my @strings = @_;
	my %hash;
	foreach(@strings){
		$hash{$_} += 1;
	}
	my @uniq = keys(%hash);
	return @uniq;
}	
1;
__END__
