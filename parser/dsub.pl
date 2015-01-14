# Contact     : daweih@me.com
# Date        : 
# Last Update : 
# Reference   : 
# Description : 

#===============================================================================================================
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

my %opts;
GetOptions(\%opts,'usr:s','max_job:s');
my $usage= <<"USAGE";
	Program: $0
	INPUT:
		-option			some option content
	OUTPUT:
	
USAGE

die $usage unless ($opts{usr});
#die $usage unless ($opts{max_job});
$opts{max_job} = 100;
my $startTime=localtime();
################################ Main ##########################################################################
# for i in `qstat -u huangdw| grep huangdw| awk -F'.' '{print $1}'`; do qdel $i; done
# for i in `find -name *gfp*`; do  rm $i; done

my $pwd = `pwd`;
chomp($pwd);
$pwd =~s/\/$//;
print $pwd, "\n";

my $files = `find -name shell.sh`;
my @files = split /\n/,$files;
foreach(@files){
	$_ =~ s/\.\///;
}
foreach(0..$#files){
	last if($_ > $opts{max_job} - 1);
	my $shell = shift @files;
	my $cmd = "dsub $pwd/". $shell;
	print $cmd, "\n";
	system $cmd;
}
sleep 100;
print "Remains ", scalar(@files), " jobs.\n";

my $qstat = `qstat -u $opts{usr} | grep $opts{usr}`;
my @qstat = split /\n/,$qstat;
while($qstat){
	last if(scalar(@files)==0);

	print scalar(@qstat), "\n";
	if( scalar(@qstat)<$opts{max_job} && scalar(@files)>0 ){
		foreach( 1..($opts{max_job}-scalar(@qstat)) ){
			last if( scalar(@files)==0 );

			my $shell = shift @files;
			my $cmd = "dsub $pwd/". $shell;	
			print $_, "\t", $cmd, "\n";
			system $cmd;
		}
	}
	print scalar(@files),"\n";
	sleep 50;
	$qstat = `qstat -u $opts{usr} | grep $opts{usr}`;
	@qstat = split /\n/,$qstat;
}

################################ Main ##########################################################################
#===============================================================================================================
my $options="-usr $opts{usr}";
my $endTime=localtime();
my $Program = $1 if($0 =~ /(.*)\.pl/);
open  LOG,">>$Program\_ProgramRunning.Log";
print LOG "From \<$startTime\> to \<$endTime\>\tperl $0 $options\n";
close LOG;

