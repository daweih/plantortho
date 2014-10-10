# Contact     : 
# Date        : Fri Oct 10 15:01:39 CST 2014
# Last Update : 
# Reference   : 
# Description : 
#===============================================================================================================
use FindBin qw($Bin); use lib "$Bin/../bin/pm";
use strict;
use warnings;

use aa_nt_info;
my $test =  &aa_nt_info("*", "aa") ;
print $test->{"full"};

__END__
