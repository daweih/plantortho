# Contact     : daweimhuang@gmail.com
# Date        : 
# Last Update : Sat Nov 19 10:28:34 AST 2011 by daweimhuang@gmail.com: changeFormat2pm
# Reference   : 
# Description : 
# Require     :

package round;
use strict;
use warnings;
require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(round);
our $VERSION = '0.01';

#Usage: $result=round($number, $decimalDigits);
#=======round()=similar to C language function==================================================================
sub round{
	my ($number,$decimalDigits)=@_;

	return $number if($number !~/\d+/);

	my $times=1;
	for (my $i=1; $i<=$decimalDigits;$i++) { $times*=10; }
	$number=($number*$times + 0.5)/$times;

	$number.="." if ($number !~/\./ && $decimalDigits !=0);
	for (my $i=1; $i<=$decimalDigits;$i++) { $number.="0"; }
	$number=($number=~/^([\+\-]*\d+\.\d{$decimalDigits})/) ? $1 : $number;

	return $number;
}
1;
__END__
