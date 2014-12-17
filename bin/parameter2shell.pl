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
GetOptions(\%opts,'pignore:s','shell:s','name:s');
my $usage= <<"USAGE";
	Program: $0
	INPUT:
		-option			some option content
	OUTPUT:
	
USAGE

die $usage unless ($opts{shell});
die $usage unless ($opts{name});
my $startTime=localtime();
################################ Main ##########################################################################
#perl parameter2shell.pl -pignore d -shell shell.sh -name phyml

my @parameter_ignore = split /,/,$opts{pignore};
my $parameter_ignore;
foreach(@parameter_ignore){
	$parameter_ignore->{ "-". $_ } = 1;
}
open I, "< ". $opts{shell};
while(<I>){
	chomp;
	next if($_ eq "");
	my @shellname;
	my @r = split /\s+/,$_;
	my $filename;
	foreach(0..$#r-1){
		next if( $r[$_] !~ /^-/ );
		next if( $parameter_ignore->{$r[$_]} );
		$r[$_+1] = $1 if($r[$_+1] =~ /.*\/(.*)$/);
		my $parameter = $r[$_].$r[$_+1];
		$parameter =~ s/^-//;
		unshift @shellname, $parameter;
		$filename = $r[$_+1] if( $r[$_] eq "-i" );
	}
	my $input_name = join ".", @shellname;
	push @shellname, "sh";
	my $shellname = join ".", @shellname;
	print "cp $filename $input_name\n";
	$_ =~ s/$filename/$input_name/;
	$shellname = $opts{name}. "-". $shellname;
	open O, "> ". $shellname;
	print O "#!/bin/sh\n";
	print O "#PBS -N $shellname\n";
	print O "#PBS -q bioque\n";
	print O "#PBS -l mem=8gb,walltime=120:00:00\n";
	print O "#HSCHED -s hschedd\n";
	print O "cd /leofs/zhangz_group/huangdw/plantortho/pipeline/OS09G0133200-01\n";
	print O "$_\n";
	close O;
}
close I;

################################ Main ##########################################################################
#===============================================================================================================
my $options="-pignore $opts{pignore} -shell $opts{shell} -name $opts{name} ";
my $endTime=localtime();
my $Program = $1 if($0 =~ /(.*)\.pl/);
open  LOG,">>$Program\_ProgramRunning.Log";
print LOG "From \<$startTime\> to \<$endTime\>\tperl $0 $options\n";
close LOG;
__END__
