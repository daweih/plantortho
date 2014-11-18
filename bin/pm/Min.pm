# Contact     : daweimhuang@gmail.com
# Date        : 
# Last Update : Sat Nov 19 10:28:34 AST 2011 by daweimhuang@gmail.com: changeFormat2pm
# Reference   : 
# Description : 

package Min;#
use strict;
use warnings;
require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(Min);#
our $VERSION = '2.0';

#get the min in @x
#Useage: $MinItm = &min(@x);
#version: 1
#update : ?
#===============================================================================================================
sub Min { 
    my($min_so_far) = shift @_; 
    foreach (@_){ 
        if($_<$min_so_far){ 
            $min_so_far=$_; 
        } 
    } 
    $min_so_far; 
} 
1;
__END__
