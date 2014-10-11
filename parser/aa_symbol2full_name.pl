# Contact     : daweih@me.com
# Date        : Fri Oct 10 15:35:14 CST 2014
# Last Update : 
# Reference   : 
# Description : 

#===============================================================================================================
use FindBin qw($Bin); use lib "$Bin/../bin/pm";

use strict;
use warnings;
use Getopt::Long;
use Bio::SeqIO;
use Data::Dumper;
use Bio::AlignIO;

use aa_nt_info;


my $file_name = $ARGV[0];
=cut
my $file_extension;

        my $fa = Bio::SeqIO->new(-file =>  $file_name,
                                                 -format => 'Fasta');
        while ( my $seq = $fa->next_seq() ) {
			
			my $index = 0;
			foreach my $aa (split //,$seq->seq){
				$index++;
#				print $aa, "\n";
				my $aa_nt_info = &aa_nt_info($aa, "aa");
#				print $index, "\t", $aa_nt_info->{abb}, "\n";
			}
			last;
        }
use Bio::Structure::IO;
my $in  = Bio::Structure::IO->new(-file => $ARGV[1],
                                          -format => 'pdb');
           while ( my $struc = $in->next_structure() ) {
              print "Structure ", $struc->id, " number of models: ",
                    Dumper( $struc->residue),"\n";
           }
=cut           
my $clustal = Bio::AlignIO->new( -file   => $file_name,
                                 -format => "clustalw" );
while(my $aln = $clustal->next_aln) {
	 foreach my $seq ($aln->each_seq) {
	 		print $seq->seq, "\n";
#             $res = $seq->subseq($pos, $pos);
#            $count{$res}++;
         }
}                                 
__END__