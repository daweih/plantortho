# Contact     : daweih@me.com
# Date        : 
# Last Update : 
# Reference   : 
# Description : 

#===============================================================================================================
use lib "/opt/perl5.12.3/lib/site_perl/5.12.3";
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;


my %opts;
GetOptions(\%opts,'newick:s');
my $usage= <<"USAGE";
	Program: $0
	INPUT:
		-newick			some option content
	OUTPUT:
	
USAGE

die $usage unless ($opts{newick});
my $startTime=localtime();
################################ Main ##########################################################################
open I, "< id_index";
my $id_index;
while(<I>){
	chomp;
	my @r = split /\t/,$_;
	$id_index->{$r[0]} = $r[1];
}
close I;

open I, "< ". $opts{newick};
while(<I>){
	foreach my $index (keys %{$id_index}){
		$_ =~ s/([\(|\)|,]{1})$index:/$1$id_index->{$index}:/;
#		$_ =~ s/([\(|\)|,]{1})$index:/$1:/;
#		print $1, "\t" if( $_ =~ /[\(|\)|,]{1}($index):/ );
#		print $index, "\t", $id_index->{$index}, "\n";
		
	}
	print $_;
}
close I;

################################ Main ##########################################################################
#===============================================================================================================
my $options="-newick $opts{newick}";
my $endTime=localtime();
my $Program = $1 if($0 =~ /(.*)\.pl/);
open  LOG,">>$Program\_ProgramRunning.Log";
print LOG "From \<$startTime\> to \<$endTime\>\tperl $0 $options\n";
close LOG;
__END__
