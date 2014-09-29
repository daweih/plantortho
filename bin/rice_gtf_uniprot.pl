# Contact     : daweih@me.com
# Date        : may 21, 2014
# Last Update : Mon Sep 29 14:09:06 CST 2014
# Reference   : daweih@me.com
# Description : generate data for IC4R seq database, rice.

#===============================================================================================================
use lib "/opt/perl5.12.3/lib/site_perl/5.12.3";
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

use Bio::SeqIO;
use Data::Dumper;
use Array::Unique;
#gff3
use Bio::GFF3::LowLevel qw/ gff3_parse_feature /;
#uniprot
use XML::XPath;
use XML::XPath::XMLParser;
use XML::XPath::Node::Attribute;

my %opts;
GetOptions(\%opts,'dir_ensembl_fasta_dna:s','rapdb_locus_gtf:s','ensembl_gtf:s','uniprot_xml:s');
my $usage= <<"USAGE";
	Program: $0
	INPUT:
			-dir_ensembl_fasta_dna
			-rapdb_locus_gtf
			-ensembl_gtf
			-uniprot_xml
	OUTPUT:
	
USAGE

die $usage unless ($opts{dir_ensembl_fasta_dna});
die $usage unless ($opts{rapdb_locus_gtf});
die $usage unless ($opts{ensembl_gtf});
die $usage unless ($opts{uniprot_xml});
my $startTime=localtime();
################################ Main ##########################################################################

#1
my $dirname = "< ". $opts{dir_ensembl_fasta_dna};
#my $dirname = "../sample/ensembl/release-22/oryza_sativa/fasta/dna";
my $genome = &genome_seq($dirname);#print length($genome->{1}->seq);
sub genome_seq {
	my $dirname = pop @_;
	my $genome_seq;
	opendir ( DIR, $dirname ) || die "Error in opening dir $dirname\n";
	while( (my $filename = readdir(DIR))){
		if( $filename =~ /.*\.fa|.*\.fasta$/ ){#&& ( $filename eq "10.fasta" || $filename eq "11.fasta" ) 
	#		print("$filename\n");
			use Bio::SeqIO;
			my $in  = Bio::SeqIO->new(-file => $dirname. "/". $filename , -format => 'Fasta');
			while ( my $seq = $in->next_seq() ) {
				$genome_seq->{$seq->primary_id} = $seq;
			}
		}
	}
	closedir(DIR);
	return $genome_seq;
}
sub SeqFetch{
	my ($loc, $genome) = @_;
	my $chr;my $strand;my $st;my $ed;
	if($loc =~ /^(\w*\d*)([+|-])(\d+)\.\.(\d+)$/){
		$chr    = $1;
		$strand = $2;
		$st     = +$3;
		$ed     = +$4;
	}
	@_=();
	push @_, $chr    ;
	push @_, $strand ;
	push @_, $st     ;
	push @_, $ed     ;
	#print "my $chr;my $strand;my $st;my $ed;";
	my $seq = $genome->{$chr}->seq;
	my $str;
	if($_[1]=~/\+/){
		if($_[2]<$_[3]){
			$str = substr($seq,$_[2]-1,$_[3]-$_[2]+1);
#			$str = $fa_seq_obj->subseq($_[2],$_[3]);
		}
		else{
			$str = substr($seq,$_[3]-1,$_[2]-$_[3]+1);
			$str = reverse($str);
		}
	}
	else{
#		$seq = reverse $seq;
        	$seq =~ tr/agctAGCT/tcgaTCGA/;
		if($_[2]<$_[3]){
			$str = substr($seq,$_[2]-1,$_[3]-$_[2]+1);
			$str = reverse($str);
		}
		else{
			$str = substr($seq,$_[3]-1,$_[2]-$_[3]+1);
			$str = reverse($str);
		}
	}
return $str;
}
sub gc_content{
	my ($string)=@_;

	$string=uc $string;

	my $G=$string; $G=~s/[^G]//g; $G=length $G;
	my $C=$string; $C=~s/[^C]//g; $C=length $C;
	my $L=$string; $L=~s/N//g; $L=length $L;

	my $GC=($L==0) ? 0 : 100*($G+$C)/$L;
	   $GC = round($GC, 3);
	return ($GC);
}
sub purine_content {
	my ($string)=@_;

	$string=uc $string;

	my $G=$string; $G=~s/[^G]//g; $G=length $G;
	my $A=$string; $A=~s/[^A]//g; $A=length $A;
	my $L=$string; $L=~s/N//g; $L=length $L;

	my $GA=($L==0) ? 0 : 100*($G+$A)/$L;
	   $GA = round($GA, 3);
	return ($GA);
}
sub round{
	my ($number,$decimalDigits)=@_;

	return $number if($number !~/\d+/);

	my $times=1;
	for (my $i=1; $i<=$decimalDigits;$i++) { $times*=10; }
	$number=($number*$times + 0.5)/$times;

	$number.="." if ($number !~/\./ && $decimalDigits !=0);
	for (my $i=1; $i<=$decimalDigits;$i++) { $number.="0"; }
	$number=($number=~/^([\+\-]*\d+\.\d{$decimalDigits})/) ? $1 : $number;
	my $num += $number;
	return $num;
}
sub uniq_array {
	my $array_ref = pop(@_);
	if($array_ref){
		my @uniq_array;
		tie @uniq_array, 'Array::Unique';
		push @uniq_array, @{$array_ref};
		return @uniq_array;
	}
	else{
		return "-";
	}
}

#2
my $rapdb;
open I, "< ". $opts{rapdb_locus_gtf};
#open I, "< ../sample/rapdb/IRGSP-1.0_representative_2014-03-05/locus.gff";
while(<I>) {
	my @r = split /\t/, $_;
	#chr01	irgsp1_locus	gene	2983	10815	.	+	.	ID=Os01g0100100;Name=Os01g0100100;Note=RabGAP/TBC domain containing protein. (Os01t0100100-01);Transcript variants=Os01t0100100-01
	my $chr       = +$1 if($r[0] =~ /^chr(\d+)$/);
	my $st        = $r[3];
	my $ed        = $r[4];
	my $strand    = $r[6];
	my @attribute = split /;/,$r[8];
	my $gene_id;
	foreach my $attribute ( @attribute ) {
		if( $attribute =~ /^ID=(\S+)/ ) {
			$gene_id = uc($1);
		}
		if( $attribute =~ /^Transcript variants=(\S+)/ ) {
			#Transcript variants=Os01t0101600-01,Os01t0101600-02,Os01t0101600-03
			my @transcript_id = split /,/,$1;
			foreach my $transcript_id ( @transcript_id ) {
				$rapdb->{$gene_id}->{uc($transcript_id)}->{check} = 1;
			}
		}
	}
}
close I;

#3
open I, "< ". $opts{ensembl_gtf};
#open I, "< ../sample/ensembl/release-22/oryza_sativa/gtf/test.gtf";#chr9
#open I, "< /Users/daweih/Archieve/database/ensembl/release-22/oryza_sativa/gtf/Oryza_sativa.IRGSP-1.0.22.gtf";
=cut
#split_start_codon
#open I, "< 92_human_gene_with_2_cds_start_site_exon_162_transcripts/ENSG00000110958.gtf";#ENST00000262033
#open I, "< 92_human_gene_with_2_cds_start_site_exon_162_transcripts/ENSG00000015479.gtf";#ENST00000504203
#open I, "< 92_human_gene_with_2_cds_start_site_exon_162_transcripts/ENSG00000091106.gtf";#ENST00000360906
#open I, "< 92_human_gene_with_2_cds_start_site_exon_162_transcripts/ENSG00000004866.gtf";#ENST00000393443
#split_stop_codon
#open I, "< 91_human_gene_with_split_stop_codon_exon_133_transcripts/ENSG00000087274.gtf";#ENST00000398125
#open I, "< 91_human_gene_with_split_stop_codon_exon_133_transcripts/ENSG00000048405.gtf";#ENST00000393313
#open I, "< 91_human_gene_with_split_stop_codon_exon_133_transcripts/ENSG00000178922.gtf";#ENST00000372432
#open I, "< 91_human_gene_with_split_stop_codon_exon_133_transcripts/ENSG00000144031.gtf";#ENST00000272421
=cut
my $gtf;
while(<I>){
	chomp;
	my @r = split /\t/,$_;
	my $chr       = $r[0];
	my $st        = $r[3];
	my $ed        = $r[4];
	my $strand    = $r[6];
	my $frame     = $r[7];
	my $phase;
	my @attribute = split /;\s*/,$r[8];

	if( $r[1] eq "protein_coding" ) {
		my $gene_id = $1       if($attribute[0] =~ /gene_id\s\"(\S+)\"/);
		my $transcript_id = $1 if($attribute[1] =~ /transcript_id\s\"(\S+)\"/);
		my $exon_no += $1 if($attribute[2] =~ /exon_number\s*\"(\d+)\"/);

		next if( !defined $rapdb->{$gene_id}->{$transcript_id}->{check} );
		$rapdb->{$gene_id}->{$transcript_id}->{ensembl} = 1;
		
		
		if( $r[2] eq "exon" ) {

			$gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{chr}    = $chr;	
			$gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{strand} = $strand;	
			$gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{st}     = $st;	
			$gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{ed}     = $ed;	

			my $loc = $chr. $strand. $st ."..". $ed;
			$gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{loc} = $loc;
		}
		elsif( $r[2] eq "CDS" ) {
			next if(!defined $gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{chr} );

			$gtf->{$gene_id}->{CDS}->{$transcript_id}->{$exon_no}->{chr}    = $chr;	
			$gtf->{$gene_id}->{CDS}->{$transcript_id}->{$exon_no}->{strand} = $strand;	
			$gtf->{$gene_id}->{CDS}->{$transcript_id}->{$exon_no}->{st}     = $st;	
			$gtf->{$gene_id}->{CDS}->{$transcript_id}->{$exon_no}->{ed}     = $ed;	
			$gtf->{$gene_id}->{CDS}->{$transcript_id}->{$exon_no}->{frame}  = $frame;	

			my $loc = $chr. $strand. $st ."..". $ed;
			$gtf->{$gene_id}->{CDS}->{$transcript_id}->{$exon_no}->{loc} = $loc;
		}
		elsif( $r[2] eq "start_codon" ){
			next if(!defined $gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{chr} );
			if( $frame == 0 ){
				$gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{start_codon} = 1;
				$gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{split_start_codon} = 1 if( $ed-$st < 2 );
				#human
				#27434 0
				#62645 1
				# 162 2
				#0: has no start_codon information; 1: has one start_codon information for one transcript; 2: start_codon was split at exon junction
			
				$gtf->{$gene_id}->{start_codon}->{$transcript_id}->{chr}    = $chr;	
				$gtf->{$gene_id}->{start_codon}->{$transcript_id}->{strand} = $strand;	
				$gtf->{$gene_id}->{start_codon}->{$transcript_id}->{st}     = $st;	
				$gtf->{$gene_id}->{start_codon}->{$transcript_id}->{ed}     = $ed;	
			}
		}
		elsif( $r[2] eq "stop_codon" ){
			next if(!defined $gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{chr} );
			if( $ed - $st == 2 || $frame != 0 ){
				$gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{stop_codon} = 1;

				$gtf->{$gene_id}->{stop_codon}->{$transcript_id}->{chr}    = $chr;	
				$gtf->{$gene_id}->{stop_codon}->{$transcript_id}->{strand} = $strand;	
				$gtf->{$gene_id}->{stop_codon}->{$transcript_id}->{st}     = $st;	
				$gtf->{$gene_id}->{stop_codon}->{$transcript_id}->{ed}     = $ed;
			}
		}
	}
}
close I;
=cut
open OR, "> rapdb.only.txt";
open OE, "> ensembl.only.txt";
open ORE, "> rapdb.ensembl.txt";
foreach my $gene_id (keys %{$rapdb}) {
	foreach my $transcript_id (keys %{$rapdb->{$gene_id}}) {
#		print Dumper($rapdb->{$gene_id}->{$transcript_id});
		if( $rapdb->{$gene_id}->{$transcript_id}->{ensembl} && $rapdb->{$gene_id}->{$transcript_id}->{check} ){
#			print "both\n";
			print ORE "$gene_id\t$transcript_id\n";
		}
		elsif( !$rapdb->{$gene_id}->{$transcript_id}->{ensembl} && $rapdb->{$gene_id}->{$transcript_id}->{check} ) {
#			print "rapdb\n";		
			print OR "$gene_id\t$transcript_id\n";
		}
		elsif( $rapdb->{$gene_id}->{$transcript_id}->{ensembl} && !$rapdb->{$gene_id}->{$transcript_id}->{check} ) {
#			print "ensembl\n";
			print OE "$gene_id\t$transcript_id\n";
		}
	
}}
close OR;
close OE; 
close ORE;
=cut

foreach my $gene_id (keys %{$gtf}) {
	foreach my $transcript_id (keys %{$gtf->{$gene_id}->{exon}}) {
		my $start_codon_exon = 0;
		my $stop_codon_exon  = 0;
		my $last_exon_no     = 0;
		foreach my $exon_no (sort {$a<=>$b} keys %{$gtf->{$gene_id}->{exon}->{$transcript_id}}) {
			$start_codon_exon = $exon_no if($gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{start_codon});
			$stop_codon_exon  = $exon_no if($gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{stop_codon});
			$last_exon_no     = $exon_no;
		}
		my $strand = $gtf->{$gene_id}->{exon}->{$transcript_id}->{$last_exon_no}->{strand};
		my $chr = $gtf->{$gene_id}->{exon}->{$transcript_id}->{$last_exon_no}->{chr};
		if( $start_codon_exon && $stop_codon_exon ) {
			$gtf->{$gene_id}->{transcript}->{$transcript_id}->{full_cds} = 1;
			foreach my $exon_no (sort {$a<=>$b} keys %{$gtf->{$gene_id}->{exon}->{$transcript_id}}) {
				if( defined $gtf->{$gene_id}->{CDS}->{$transcript_id}->{$exon_no}->{chr} ) {
					$gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{cds} = 1;
				}
				else{
					if( $exon_no < $start_codon_exon ) {
						$gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{five_prime_UTR}  = 1;
						$gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{cds}             = 0;
						$gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{three_prime_UTR} = 0;
					}
					if( $exon_no > $stop_codon_exon ) {
						$gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{five_prime_UTR}  = 0;
						$gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{cds}             = 0;
						$gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{three_prime_UTR} = 1;
					}
				}

				if($gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{start_codon}){
					if( $strand eq "+" ) {
						if( $gtf->{$gene_id}->{start_codon}->{$transcript_id}->{st} == $gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{st} ) {
							$gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{five_prime_UTR} = 0;
						}
						elsif( $gtf->{$gene_id}->{start_codon}->{$transcript_id}->{st} > $gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{st} ) {
							$gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{five_prime_UTR} = 1;
						
							$gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{sub_structure}->{five_prime_UTR}->{st} = $gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{st};
							$gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{sub_structure}->{five_prime_UTR}->{ed} = $gtf->{$gene_id}->{start_codon}->{$transcript_id}->{st} - 1 ;
							my $loc = "$chr$strand". $gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{sub_structure}->{five_prime_UTR}->{st}."..". $gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{sub_structure}->{five_prime_UTR}->{ed};
							$gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{sub_structure}->{five_prime_UTR}->{loc} = $loc;
						}
					}
					elsif( $strand eq "-" ) {
						if( $gtf->{$gene_id}->{start_codon}->{$transcript_id}->{ed} == $gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{ed} ) {
							$gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{five_prime_UTR} = 0;
						}
						elsif( $gtf->{$gene_id}->{start_codon}->{$transcript_id}->{ed} < $gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{ed} ) {
							$gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{five_prime_UTR} = 1;

							$gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{sub_structure}->{five_prime_UTR}->{st} = $gtf->{$gene_id}->{start_codon}->{$transcript_id}->{ed} + 1 ;
							$gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{sub_structure}->{five_prime_UTR}->{ed} = $gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{ed};
							my $loc = "$chr$strand". $gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{sub_structure}->{five_prime_UTR}->{st}."..". $gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{sub_structure}->{five_prime_UTR}->{ed};
							$gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{sub_structure}->{five_prime_UTR}->{loc} = $loc;
						}
					}
				}
			
				if($gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{stop_codon}) {

					if( $strand eq "+" ) {
						if( $gtf->{$gene_id}->{stop_codon}->{$transcript_id}->{ed} == $gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{ed} ) {
							$gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{three_prime_UTR} = 0;
						}
						elsif( $gtf->{$gene_id}->{stop_codon}->{$transcript_id}->{ed} < $gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{ed} ) {
							$gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{three_prime_UTR} = 1;
						
							$gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{sub_structure}->{three_prime_UTR}->{st} = $gtf->{$gene_id}->{stop_codon}->{$transcript_id}->{ed} - 2;
							$gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{sub_structure}->{three_prime_UTR}->{ed} = $gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{ed};
							my $loc = "$chr$strand". $gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{sub_structure}->{three_prime_UTR}->{st}."..". $gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{sub_structure}->{three_prime_UTR}->{ed};
							$gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{sub_structure}->{three_prime_UTR}->{loc} = $loc;
						}
					}
					elsif( $strand eq "-" ) {
						if( $gtf->{$gene_id}->{stop_codon}->{$transcript_id}->{st} == $gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{st} ) {
							$gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{three_prime_UTR} = 0;
						}
						elsif( $gtf->{$gene_id}->{stop_codon}->{$transcript_id}->{st} > $gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{st} ) {
							$gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{three_prime_UTR} = 1;

							$gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{sub_structure}->{three_prime_UTR}->{st} = $gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{st};
							$gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{sub_structure}->{three_prime_UTR}->{ed} = $gtf->{$gene_id}->{stop_codon}->{$transcript_id}->{st} + 2;
							my $loc = "$chr$strand". $gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{sub_structure}->{three_prime_UTR}->{st}."..". $gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{sub_structure}->{three_prime_UTR}->{ed};
							$gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{sub_structure}->{three_prime_UTR}->{loc} = $loc;
						}
					}
				}
				$gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{five_prime_UTR}  = 0 if( !exists $gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{five_prime_UTR}  );
				$gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{cds}             = 0 if( !exists $gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{cds}             );
				$gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{three_prime_UTR} = 0 if( !exists $gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{three_prime_UTR} );

				my $phase;
				if( defined $gtf->{$gene_id}->{CDS}->{$transcript_id}->{$exon_no}->{frame} ) {
					my $frame = $gtf->{$gene_id}->{CDS}->{$transcript_id}->{$exon_no}->{frame};
					if( !$gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{start_codon} ){
						if( $frame == 0){
							$phase = 0;
						}
						elsif( $frame == 2 ){
							$phase = 1;
						}
						elsif( $frame == 1 ){
							$phase = 2;
						}
						$gtf->{$gene_id}->{CDS}->{$transcript_id}->{$exon_no}->{phase} = $phase;
						$gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{phase} = $phase;
					}
				}
			}
			foreach my $exon_no (sort {$a<=>$b} keys %{$gtf->{$gene_id}->{exon}->{$transcript_id}}) {
				$gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{phase} = -1 if( !defined $gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{phase} );
			}
		}
		else {
			$gtf->{$gene_id}->{transcript}->{$transcript_id}->{full_cds} = 0;
		}
	}
}
foreach my $gene_id (keys %{$gtf}) {
	#print $gene_id, "\n";
	#transcript
	foreach my $transcript_id (keys %{$gtf->{$gene_id}->{exon}}) {
		my $last_exon_no;
		foreach my $exon_no (sort {$a<=>$b} keys %{$gtf->{$gene_id}->{exon}->{$transcript_id}}) {
			$last_exon_no = $exon_no;
		}
		
		$gtf->{$gene_id}->{transcript}->{$transcript_id}->{chr}    = $gtf->{$gene_id}->{exon}->{$transcript_id}->{$last_exon_no}->{chr};
		$gtf->{$gene_id}->{transcript}->{$transcript_id}->{strand} = $gtf->{$gene_id}->{exon}->{$transcript_id}->{$last_exon_no}->{strand};
		
		if( $gtf->{$gene_id}->{exon}->{$transcript_id}->{$last_exon_no}->{strand} eq "+" ){
			$gtf->{$gene_id}->{transcript}->{$transcript_id}->{st}     = $gtf->{$gene_id}->{exon}->{$transcript_id}->{1}->{st};	
			$gtf->{$gene_id}->{transcript}->{$transcript_id}->{ed}     = $gtf->{$gene_id}->{exon}->{$transcript_id}->{$last_exon_no}->{ed};	
		}
		else{
			$gtf->{$gene_id}->{transcript}->{$transcript_id}->{st}     = $gtf->{$gene_id}->{exon}->{$transcript_id}->{$last_exon_no}->{st};
			$gtf->{$gene_id}->{transcript}->{$transcript_id}->{ed}     = $gtf->{$gene_id}->{exon}->{$transcript_id}->{1}->{ed};	
		}
		
		my $loc = $gtf->{$gene_id}->{transcript}->{$transcript_id}->{chr}. $gtf->{$gene_id}->{transcript}->{$transcript_id}->{strand}. $gtf->{$gene_id}->{transcript}->{$transcript_id}->{st} ."..". $gtf->{$gene_id}->{transcript}->{$transcript_id}->{ed};
		$gtf->{$gene_id}->{transcript}->{$transcript_id}->{loc} = $loc;
	}
	
	#gene
	foreach my $transcript_id (keys %{$gtf->{$gene_id}->{exon}}) {
		$gtf->{$gene_id}->{gene}->{chr}    = $gtf->{$gene_id}->{transcript}->{$transcript_id}->{chr};
		$gtf->{$gene_id}->{gene}->{strand} = $gtf->{$gene_id}->{transcript}->{$transcript_id}->{strand};

		if( !defined $gtf->{$gene_id}->{gene}->{st}){
			$gtf->{$gene_id}->{gene}->{st} = $gtf->{$gene_id}->{transcript}->{$transcript_id}->{st};
		}
		elsif( $gtf->{$gene_id}->{gene}->{st} > $gtf->{$gene_id}->{transcript}->{$transcript_id}->{st} ){
			$gtf->{$gene_id}->{gene}->{st} = $gtf->{$gene_id}->{transcript}->{$transcript_id}->{st};
		}
		
		if( !defined $gtf->{$gene_id}->{gene}->{ed}){
			$gtf->{$gene_id}->{gene}->{ed} = $gtf->{$gene_id}->{transcript}->{$transcript_id}->{ed};
		}
		elsif( $gtf->{$gene_id}->{gene}->{ed} < $gtf->{$gene_id}->{transcript}->{$transcript_id}->{ed} ){
			$gtf->{$gene_id}->{gene}->{ed} = $gtf->{$gene_id}->{transcript}->{$transcript_id}->{ed};
		}
	}
	my $loc      =  $gtf->{$gene_id}->{gene}->{chr}. $gtf->{$gene_id}->{gene}->{strand}. $gtf->{$gene_id}->{gene}->{st} ."..". $gtf->{$gene_id}->{gene}->{ed};
	$gtf->{$gene_id}->{gene}->{loc} = $loc;
	
	#intron/exon
	foreach my $transcript_id (keys %{$gtf->{$gene_id}->{exon}}) {
		my $last_exon_no;
		foreach my $exon_no (sort {$a<=>$b} keys %{$gtf->{$gene_id}->{exon}->{$transcript_id}}) {
			$last_exon_no = $exon_no;
		}
		foreach my $exon_no (sort {$a<=>$b} keys %{$gtf->{$gene_id}->{exon}->{$transcript_id}}) {
		}
		foreach my $exon_no (sort {$a<=>$b} keys %{$gtf->{$gene_id}->{exon}->{$transcript_id}}) {
			my $next_exon_no = $exon_no + 1;
			next if( $exon_no == $last_exon_no );
			$gtf->{$gene_id}->{intron}->{$transcript_id}->{$exon_no}->{chr}    = $gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{chr}   ;
			$gtf->{$gene_id}->{intron}->{$transcript_id}->{$exon_no}->{strand} = $gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{strand};
			if( $gtf->{$gene_id}->{intron}->{$transcript_id}->{$exon_no}->{strand} eq "+" ){
				$gtf->{$gene_id}->{intron}->{$transcript_id}->{$exon_no}->{st}     = $gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{ed} + 1     ;
				$gtf->{$gene_id}->{intron}->{$transcript_id}->{$exon_no}->{ed}     = $gtf->{$gene_id}->{exon}->{$transcript_id}->{$next_exon_no}->{st} - 1;
			}else{
				$gtf->{$gene_id}->{intron}->{$transcript_id}->{$exon_no}->{st}     = $gtf->{$gene_id}->{exon}->{$transcript_id}->{$next_exon_no}->{ed} + 1     ;
				$gtf->{$gene_id}->{intron}->{$transcript_id}->{$exon_no}->{ed}     = $gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{st} - 1;
			}
			my $loc = $gtf->{$gene_id}->{intron}->{$transcript_id}->{$exon_no}->{chr}. $gtf->{$gene_id}->{intron}->{$transcript_id}->{$exon_no}->{strand}. $gtf->{$gene_id}->{intron}->{$transcript_id}->{$exon_no}->{st} ."..". $gtf->{$gene_id}->{intron}->{$transcript_id}->{$exon_no}->{ed};
			$gtf->{$gene_id}->{intron}->{$transcript_id}->{$exon_no}->{loc} = $loc;
			
			$gtf->{$gene_id}->{intron}->{$transcript_id}->{$exon_no}->{phase} = $gtf->{$gene_id}->{exon}->{$transcript_id}->{$next_exon_no}->{phase};
			
		}
	}	
}


my $gene_feature = $gtf;

open PF, "> protein.fa";

#4
my $ref;
open I, "< ". $opts{uniprot_xml};
#open I, "< ../sample/uniprot/test.uniprot.xml";
#open I, "< Archieve/database/uniprot/39947_uniprot_sprot.select_oryza_sativa_japonica.xml";
my $xml;
while(<I>){
	chomp;
	if(/^<entry.*/){
		$xml = "<?xml version='1.0' encoding='UTF-8'?>\n<uniprot xmlns=\"http://uniprot.org/uniprot\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://uniprot.org/uniprot http://www.uniprot.org/support/docs/uniprot.xsd\">\n";
	}
	
	$xml .= $_ ."\n";
	
	if($_ eq "</entry>"){
		$xml .= "<copyright>\nCopyrighted by the UniProt Consortium, see http://www.uniprot.org/terms\nDistributed under the Creative Commons Attribution-NoDerivs License\n</copyright>\n</uniprot>";
		# Initial
		my $whole_xml = XML::XPath->new(xml=>$xml);
		# Using XPath to find out the node
		my $entry = $whole_xml->find('/uniprot/entry');
 		# Iterator through the node list
		foreach my $node ($entry->get_nodelist) {
			my $check_modify_gene_name = 0;
			my $gene_name = $whole_xml->find('gene', $node);my @gene_name_l = split /\n/,$gene_name;my @gene_name;
			foreach(@gene_name_l){
				next if($_ eq "");
				push @gene_name, $_;
				if($_ =~ /^Os\d+g\d+$/){
#				if($modify_gene_name->{$_}){
					$check_modify_gene_name = 1;
#					$gene_name = $_;$gene_name =~ s/Os/OS/g;$gene_name =~ s/g/G/g;
					$gene_name = uc($_);
				}
			}
			
			
			next if(!defined $gene_feature->{$gene_name});
			#print $gene_name, "\n";
			my $accession = $whole_xml->find('accession', $node);
			for(my $i = 0;$i<length($accession);$i+=6){
				my $accession_id = substr($accession, $i, 6);
				push @{$gene_feature->{$gene_name}->{uniprot_accession}}, $accession_id;
			}
			
			
			#gene_name_synomous################
			my $gene_name_synomous = join ";", @gene_name;#print $gene_name_synomous, "\n";
			push @{$gene_feature->{$gene_name}->{synomous}}, @gene_name;

			if($check_modify_gene_name){
				my $chr_id += $1 if($gene_name =~ /^Os(\d+)g\d+$/);

				my $uniprot_name = $whole_xml->find('name', $node);
			my $primary_identify = $uniprot_name;#print "\n", $uniprot_name, "\n";

			push @{$gene_feature->{$gene_name}->{uniprot_name}}, $uniprot_name;

				#sequence################
				my $sequence = $whole_xml->find('sequence', $node);$sequence =~ s/\n//g;
				print PF ">$uniprot_name\n";
				print PF "$sequence\n";

				my $uniprot_proteinname = $whole_xml->find('/uniprot/entry/protein', $node);my @uniprot_proteinname_l = split /\n/,$uniprot_proteinname;my @uniprot_proteinname;
				foreach(@uniprot_proteinname_l){
					next if($_ eq "");
					push @uniprot_proteinname, $_;
					push @{$gene_feature->{$gene_name}->{uniprot_proteinname_synomous}}, $_;

				}
				#################
				my $uniprot_proteinname_synomous = join ";", @uniprot_proteinname;#print $uniprot_proteinname_synomous, "\n";
				
				
				my $uniprot_organism = $whole_xml->find('/uniprot/entry/organism/lineage', $node);my @uniprot_organism_l = split /\n/,$uniprot_organism;my @uniprot_organism;
				foreach(@uniprot_organism_l){
					next if($_ eq "");
					push @uniprot_organism, $_;
				}
				#lineage
				$uniprot_organism = join "->", @uniprot_organism;#print $uniprot_organism, "\n";
				$gene_feature->{$gene_name}->{organism} = $uniprot_organism;

				#dbReference
				my @dbReference = $node->findnodes('dbReference'); 
				my $dbReference_type;
				my $dbReference_id;
				foreach my $element (@dbReference) {
					foreach my $attribute ($element->getAttributes) {
#				        print $element->getName, ", ", $attribute->getName, ":", $attribute->getData, "\n";
				        $dbReference_type = $attribute->getData if($attribute->getName eq "type");
						$dbReference_id = $attribute->getData   if($attribute->getName eq "id");
					}
					#annotation type				annotation value, stored in an array(@)
					#print $dbReference_type, ": ", $dbReference_id, "\n"; 
					if($dbReference_type eq "UniGene"){
						push @{$gene_feature->{$gene_name}->{annotation}->{UniGene}}, $dbReference_id;
					}
					if($dbReference_type eq "GeneID"){
						push @{$gene_feature->{$gene_name}->{annotation}->{GeneID}}, $dbReference_id;
					}
					if($dbReference_type eq "GO"){
						push @{$gene_feature->{$gene_name}->{annotation}->{GO}}, $dbReference_id;
					}


					if($dbReference_type eq "RefSeq"){
						push @{$gene_feature->{$gene_name}->{annotation}->{$uniprot_name}->{RefSeq}}, $dbReference_id;
					}
					if($dbReference_type eq "InterPro"){
						push @{$gene_feature->{$gene_name}->{annotation}->{$uniprot_name}->{InterPro}}, $dbReference_id;
					}
					if($dbReference_type eq "Pfam"){
						push @{$gene_feature->{$gene_name}->{annotation}->{$uniprot_name}->{Pfam}}, $dbReference_id;
					}
					if($dbReference_type eq "EC"){
						push @{$gene_feature->{$gene_name}->{annotation}->{$uniprot_name}->{EC}}, $dbReference_id;
					}
					if($dbReference_type eq "KEGG"){
						push @{$gene_feature->{$gene_name}->{annotation}->{$uniprot_name}->{KEGG}}, $dbReference_id;
					}
				}
				
				#reference
				my @reference = $node->findnodes('reference/citation');
				foreach my $element (@reference) {
					##
					my $ref_type = "-";
					my $date     = "-";
					my $journal  = "-";
					my $v        = "-";
					my $v_st     = "-";
					my $v_ed     = "-";
					foreach my $attribute ($element->getAttributes) {
						#print $element->getName, ", ", $attribute->getName, ":", $attribute->getData, "\n";
				        $ref_type = $attribute->getData if($attribute->getName eq "type");
				        $date     = $attribute->getData if($attribute->getName eq "date");
				        $journal  = $attribute->getData if($attribute->getName eq "name");
				        $v        = $attribute->getData if($attribute->getName eq "volume");
				        $v_st     = $attribute->getData if($attribute->getName eq "first");
				        $v_ed     = $attribute->getData if($attribute->getName eq "last");
					}
					##
					my $title = $element->find('title');#print $title, "\n";
					my @person = $element->findnodes('authorList/person');
					my @name;
					foreach my $person (@person) {
						foreach my $attribute ($person->getAttributes) {
#							print $attribute->getName, ":", $attribute->getData, "\n";
							push @name, $attribute->getData;
						}
					}
					my $names = join ", ", @name;#print $names, "\n";
					
					####
					my $consortium; 
					my @consortium = $element->findnodes('authorList/consortium');
					foreach my $consortium (@consortium){
						foreach my $attribute ($consortium->getAttributes) {
							#print $attribute->getData, "\n";
						}
					}
					
					my $id_type;my $pubmed_id = "";my $doi_id = "";
					if($ref_type eq "journal article"){
						my @ref_dbReference = $element->findnodes('dbReference');
						foreach my $ref_dbReference (@ref_dbReference){
							foreach my $attribute ($ref_dbReference->getAttributes) {
								$id_type = $attribute->getData if($attribute->getName eq "type");
								$id_type = $attribute->getData if($attribute->getName eq "type");
								$pubmed_id = $attribute->getData if($id_type eq "PubMed");
								$doi_id = $attribute->getData    if($id_type eq "DOI");
								#print $ref_dbReference->getName, "\t", $attribute->getName, ":", $attribute->getData, "\n";
							}
						}
					}
					if($ref_type eq "journal article"){
						next if($pubmed_id eq "");
						my $journal_article = $ref_type. "\t". $title. "\t". $names. "\t". $date. "\t". $journal. "\t". $v. "\t". $v_st. "\t". $v_ed. "\t". $pubmed_id. "\t". $doi_id;
						$gene_feature->{$gene_name}->{pub}->{$pubmed_id} = $journal_article;


						$ref->{$ref_type}->{$pubmed_id}->{ref_type } = $ref_type ;
						$ref->{$ref_type}->{$pubmed_id}->{title    } = $title    ;
						$ref->{$ref_type}->{$pubmed_id}->{names    } = $names    ;
						$ref->{$ref_type}->{$pubmed_id}->{date     } = $date     ;
						$ref->{$ref_type}->{$pubmed_id}->{journal  } = $journal  ;
						$ref->{$ref_type}->{$pubmed_id}->{v        } = $v        ;
						$ref->{$ref_type}->{$pubmed_id}->{v_st     } = $v_st     ;
						$ref->{$ref_type}->{$pubmed_id}->{v_ed     } = $v_ed     ;
						$ref->{$ref_type}->{$pubmed_id}->{pubmed_id} = $pubmed_id;
						$ref->{$ref_type}->{$pubmed_id}->{doi_id   } = $doi_id   ;
					}
					else{
						my $title_u = lc($title);
						$ref->{$ref_type}->{$title_u}->{ref_type} = $ref_type;
						$ref->{$ref_type}->{$title_u}->{title   } = $title   ;
						$ref->{$ref_type}->{$title_u}->{names   } = $names   ;
						$ref->{$ref_type}->{$title_u}->{date    } = $date    ;
#						print R  $ref_type ,"\t";
#						print R  $title    ,"\t";
#						print R  $names    ,"\t";
#						print R  $date     ,"\t";
#						print R  "-"  ,"\t";
#						print R  "-"        ,"\t";
#						print R  "-"     ,"\t";
#						print R  "-"     ,"\t";
#						print R  "-","\t";
#						print R  "-"   ,"\n";
					}
					
				}
#				my $entry_xml = XML::XPath::XMLParser::as_string($node);
			}
		}

		$xml = "";
	}
}
close I;
close PF;
print "#Linkout\n";
#print "##Use PubMed link: http://www.ncbi.nlm.nih.gov/pubmed/?term=pubmed_id\n";
print "##Use UniProt link: http://www.uniprot.org/uniprot/uniprot_accession\n";
print "##Use UniProt link: http://www.uniprot.org/uniprot/uniprot_entry_name2seq\n";
print "##Use NCBI link: http://www.ncbi.nlm.nih.gov/gene/?term=GeneID\n";
print "##Use KEGG link: http://www.genome.jp/dbget-bin/www_bget?KEGG_id\n";#osa:4326459
print "##Use brenda-enzymes link: http://www.brenda-enzymes.info/php/result_flat.php4?ecno=EC_no\n";#4.1.2.27
print "##Use PDB link: http://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/enzymes/GetPage.pl?ec_number=EC_no\n";#4.1.2.27
print "##Use GO link: http://amigo.geneontology.org/amigo/term/GO_no\n";#GO:0016831
print "##Use Gramene link: http://ensembl.gramene.org/Oryza_sativa/Gene/Summary?db=core;g=primaryIdentifier\n";#OS01G0101200
print "##Use Pfam link: http://pfam.xfam.org/family/Pfam_id\n";#PF13419

print "#primaryIdentifier|Gene symbol\t";
print "gene_loc|Genome location\t";
print "gene_seq|Gene full length sequence\t";
print "transcript_loc|Transcript location\t";
print "transcript_stru|Transcript structure\t";
print "transcript_seq|Transcript sequence (mRNA sequence)\t";
print "exon_strus|Exon/Intron structure\t";
print "exon_phase|Exon/Intron phase\t";
print "exon_gc|Exon/Intron GC content\t";
print "exon_purine|Exon/Intron purine content\t";
print "GeneID|NCBI Gene ID\t";#4326460
print "UniGene_id|NCBI UniGene ID\t";#Os.19507
print "synomous|Synomous\t";
print "RefSeq|NCBI Reference Sequence\t";
print "GO|Gene Ontology number\t";
print "uniprot_accession|UniProt Accession\t";
print "uniprot_proteinname_synomous|UniProt Protein Synomous\t";
print "uniprot_entry_name2seq|UniProt Entry name\t";
print "InterPro_id|InterPro Domain ID\t";
print "Pfam_id|Pfam Domain ID\t";
print "EC|EC Number\t";
print "KEGG_protein|KEGG Protein\t";
#print "pubmed_id\t";
print "organism\n";


my $gene_synomous;
foreach my $gene_id (keys %{$gtf}) {
	#print $gene_id, "\n";
	next if(!defined $gtf->{$gene_id}->{gene}->{loc} );
	my $gene_seq = &SeqFetch( $gtf->{$gene_id}->{gene}->{loc}, $genome );
	my $strand = $gtf->{$gene_id}->{gene}->{strand};
	my @transcript_stru_summary;
	my @transcript_phase_summary;
	my @transcript_seq;
	foreach my $transcript_id (keys %{$gtf->{$gene_id}->{exon}}) {
		#print $gene_id, "\t$transcript_id\n";
		my $last_exon_no;foreach my $exon_no (sort {$a<=>$b} keys %{$gtf->{$gene_id}->{exon}->{$transcript_id}}) {$last_exon_no = $exon_no;}
		my @transcript_stru = ();
		my @transcript_phase = ();
		my $transcript_seq;
		if( $gtf->{$gene_id}->{transcript}->{$transcript_id}->{full_cds} ) {
		#含有完整cds结构的
			foreach my $exon_no (sort {$a<=>$b} keys %{$gtf->{$gene_id}->{exon}->{$transcript_id}}) {
				my $next_exon_no = $exon_no + 1;
				my $unit;
				my $unit_seq;
				#five_prime_UTR
				if( $gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{five_prime_UTR} && !$gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{cds} ) {
					$unit_seq = &SeqFetch($gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{loc}, $genome );
					$unit = "five_prime_UTR_$exon_no:". $gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{loc}. ":$unit_seq";#
					push @transcript_stru, $unit;
				}
				#cds
				if( exists $gtf->{$gene_id}->{CDS}->{$transcript_id}->{$exon_no}->{loc}){
					if( $gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{five_prime_UTR} ) {
						$unit_seq = &SeqFetch($gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{sub_structure}->{five_prime_UTR}->{loc}, $genome );
						$unit = "five_prime_UTR_$exon_no:". $gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{sub_structure}->{five_prime_UTR}->{loc}. ":$unit_seq";#
						push @transcript_stru, $unit;
					}
					my $unit_seq = &SeqFetch($gtf->{$gene_id}->{CDS}->{$transcript_id}->{$exon_no}->{loc}, $genome );
					$transcript_seq .= $unit_seq;
					$unit = "cds_$exon_no:". $gtf->{$gene_id}->{CDS}->{$transcript_id}->{$exon_no}->{loc}. ":$unit_seq";#
					push @transcript_stru, $unit;
					if( $gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{three_prime_UTR} ) {
						$unit_seq = &SeqFetch($gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{sub_structure}->{three_prime_UTR}->{loc}, $genome );
						$unit = "three_prime_UTR_$exon_no:". $gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{sub_structure}->{three_prime_UTR}->{loc}. ":$unit_seq";#
						push @transcript_stru, $unit;
					}
				}
				#three_prime_UTR
				if( $gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{three_prime_UTR} && !$gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{cds}  ) {
					$unit_seq = &SeqFetch($gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{loc}, $genome );
					$unit = "three_prime_UTR_$exon_no:". $gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{loc}. ":$unit_seq";#
					push @transcript_stru, $unit;
				}
				
				push @transcript_phase, "exon_$exon_no:". $gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{phase}. ",-1"if( $exon_no == $last_exon_no );
				next if( $exon_no == $last_exon_no );
				push @transcript_phase, "exon_$exon_no:". $gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{phase}. ",". $gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no+1}->{phase}. "|intron_$exon_no:". $gtf->{$gene_id}->{exon}->{$transcript_id}->{$exon_no+1}->{phase};
				
				$unit_seq = &SeqFetch($gtf->{$gene_id}->{intron}->{$transcript_id}->{$exon_no}->{loc}, $genome );
				my $intron_stru = "intron_$exon_no:". $gtf->{$gene_id}->{intron}->{$transcript_id}->{$exon_no}->{loc}. ":$unit_seq";#
				push @transcript_stru, $intron_stru;
			}
		}
		else {
		#不具有完整cds结构的
			foreach my $exon_no (sort {$a<=>$b} keys %{$gtf->{$gene_id}->{exon}->{$transcript_id}}) {
				if( exists $gtf->{$gene_id}->{CDS}->{$transcript_id}->{$exon_no}->{loc}){
					my $seq = &SeqFetch($gtf->{$gene_id}->{CDS}->{$transcript_id}->{$exon_no}->{loc}, $genome );
					$transcript_seq .= $seq;
				}
			}
			push @transcript_stru, "-";
			push @transcript_phase, "-";
		}

		my $transcript_phase_summary = ( scalar(@transcript_phase) > 1) ? join "|", @transcript_phase : $transcript_phase[0];
		$transcript_phase_summary = "$transcript_id:". $transcript_phase_summary;
		push @transcript_phase_summary, $transcript_phase_summary;

		my $transcript_stru_summary = ( scalar(@transcript_stru) > 1) ? join "|", @transcript_stru : $transcript_stru[0];
		$transcript_stru_summary = "$transcript_id:". $transcript_stru_summary;
		push @transcript_stru_summary, $transcript_stru_summary;
		
		$transcript_seq = "-" if(!defined $transcript_seq);
		$transcript_seq = "$transcript_id:". $transcript_seq;
		push @transcript_seq, $transcript_seq;
	}
	my $print_transcript_stru = (scalar(@transcript_stru_summary) > 1) ? join ";", @transcript_stru_summary : $transcript_stru_summary[0];
	my $print_transcript_phase = (scalar(@transcript_phase_summary) > 1) ? join ";", @transcript_phase_summary : $transcript_phase_summary[0];
	my $print_transcript_seq = (scalar(@transcript_seq) > 1) ? join ";", @transcript_seq : $transcript_seq[0];
	
	my @transcript;
	my @exon_strus;
	my @exon_gc;
	my @exon_purine;
	foreach my $transcript_id (keys %{$gene_feature->{$gene_id}->{transcript}}){
		my $transcript_loc = $transcript_id. ":". $gene_feature->{$gene_id}->{transcript}->{$transcript_id}->{loc};
		push @transcript, $transcript_loc;
		my @exons;
		my @gcs;
		my @purines;
		foreach my $exon_no (sort {$a<=>$b} keys %{$gtf->{$gene_id}->{exon}->{$transcript_id}}) {
			my $unit_seq;
			$unit_seq = &SeqFetch($gene_feature->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{loc}, $genome );
			my $exon = "exon_". $exon_no. ":". $gene_feature->{$gene_id}->{exon}->{$transcript_id}->{$exon_no}->{loc}. ":$unit_seq";#
			push @exons, $exon;
			push @gcs,     "exon_". $exon_no. ":". &gc_content    ($unit_seq);
			push @purines, "exon_". $exon_no. ":". &purine_content($unit_seq);
			next if(!defined $gene_feature->{$gene_id}->{intron}->{$transcript_id}->{$exon_no}->{loc} );
			$unit_seq = &SeqFetch($gene_feature->{$gene_id}->{intron}->{$transcript_id}->{$exon_no}->{loc}, $genome );
			my $intron = "intron_". $exon_no. ":". $gene_feature->{$gene_id}->{intron}->{$transcript_id}->{$exon_no}->{loc}. ":$unit_seq";#
			push @exons, $intron;
			push @gcs    , "intron_". $exon_no. ":". &gc_content    ($unit_seq);
			push @purines, "intron_". $exon_no. ":". &purine_content($unit_seq);
		}
		my $exon_stru = join "|", @exons;
		$exon_stru = $transcript_id. ":". $exon_stru;
		push @exon_strus, $exon_stru;
		
		my $exon_gc = join "|", @gcs;
		$exon_gc = $transcript_id. ":". $exon_gc;
		push @exon_gc, $exon_gc;
		
		my $exon_purine = join "|", @purines;
		$exon_purine = $transcript_id. ":". $exon_purine;
		push @exon_purine, $exon_purine;
	}
	my $print_transcript      = (scalar(@transcript)      > 1) ? join ";", @transcript      : $transcript[0];
	my $print_exon_strus      = (scalar(@exon_strus)      > 1) ? join ";", @exon_strus      : $exon_strus[0];
	my $print_exon_gc         = (scalar(@exon_gc)         > 1) ? join ";", @exon_gc         : $exon_gc[0];
	my $print_exon_purine     = (scalar(@exon_purine)     > 1) ? join ";", @exon_purine     : $exon_purine[0];

	my @pubmed_id = keys %{$gene_feature->{$gene_id}->{pub}} if(defined $gene_feature->{$gene_id}->{pub});
	my $pubmed_id = (defined $gene_feature->{$gene_id}->{pub}) ? join ";", @pubmed_id : "-";

	my $uniprot_name = (defined $gene_feature->{$gene_id}->{uniprot_name}) ? join ";", @{$gene_feature->{$gene_id}->{uniprot_name}} : "-";
	my $uniprot_accession = (defined $gene_feature->{$gene_id}->{uniprot_name}) ? join ";", @{$gene_feature->{$gene_id}->{uniprot_accession}} : "-";

	my @print_refseq_id = ();
	my @print_InterPro_id = ();
	my @print_Pfam_id = ();
	my @print_EC_id = ();
	my @print_KEGG_id = ();
	if(defined $gene_feature->{$gene_id}->{uniprot_name}){
		foreach my $uniprot_name (@{$gene_feature->{$gene_id}->{uniprot_name}}){
	
			my @refseq_id = &uniq_array($gene_feature->{$gene_id}->{annotation}->{$uniprot_name}->{RefSeq});
			my @InterPro  = &uniq_array($gene_feature->{$gene_id}->{annotation}->{$uniprot_name}->{InterPro});
			my @Pfam      = &uniq_array($gene_feature->{$gene_id}->{annotation}->{$uniprot_name}->{Pfam});
			my @EC        = &uniq_array($gene_feature->{$gene_id}->{annotation}->{$uniprot_name}->{EC});
			my @KEGG      = &uniq_array($gene_feature->{$gene_id}->{annotation}->{$uniprot_name}->{KEGG});

			foreach(@refseq_id){$gene_synomous->{$gene_id}->{$_} = 1;}
			foreach(@InterPro ){$gene_synomous->{$gene_id}->{$_} = 1;}
			foreach(@Pfam     ){$gene_synomous->{$gene_id}->{$_} = 1;}
			foreach(@EC       ){$gene_synomous->{$gene_id}->{$_} = 1;}
			foreach(@KEGG     ){$gene_synomous->{$gene_id}->{$_} = 1;}

			my $refseq_id = (scalar(@refseq_id) > 1) ?  join "|", @refseq_id : $refseq_id[0];
			$refseq_id = $uniprot_name. ":". $refseq_id;
			push @print_refseq_id, $refseq_id;

			my $InterPro = (scalar(@InterPro) > 1) ?  join "|", @InterPro : $InterPro[0];
			my $Pfam = (scalar(@Pfam) > 1) ?  join "|", @Pfam : $Pfam[0];
			$InterPro = $uniprot_name. ":". $InterPro;
			$Pfam     = $uniprot_name. ":". $Pfam;
			push @print_InterPro_id, $InterPro;
			push @print_Pfam_id, $Pfam;

			my $EC = (scalar(@EC) > 1) ?  join "|", @EC : $EC[0];
			$EC = $uniprot_name. ":". $EC;
			push @print_EC_id, $EC;

			my $KEGG = (scalar(@KEGG) > 1) ?  join "|", @KEGG : $KEGG[0];
			$KEGG = $uniprot_name. ":". $KEGG;
			push @print_KEGG_id, $KEGG;
		}
	}
	
	my @UniGene = &uniq_array($gene_feature->{$gene_id}->{annotation}->{UniGene});
	my @GeneID  = &uniq_array($gene_feature->{$gene_id}->{annotation}->{GeneID});
	my @GO      = &uniq_array($gene_feature->{$gene_id}->{annotation}->{GO});
	foreach(@UniGene){$gene_synomous->{$gene_id}->{$_} = 1;}
	foreach(@GeneID ){$gene_synomous->{$gene_id}->{$_} = 1;}
	foreach(@GO     ){$gene_synomous->{$gene_id}->{$_} = 1;}

	my $print_GeneID = join ";", @GeneID;
	my $print_UniGene = join ";", @UniGene;
	my $print_GO = join ";", @GO;


	my $print_refseq_id   = join ";", @print_refseq_id;
	my $print_InterPro_id = join ";", @print_InterPro_id;
	my $print_Pfam_id     = join ";", @print_Pfam_id;
	my $print_EC_id       = join ";", @print_EC_id;
	my $print_KEGG_id     = join ";", @print_KEGG_id;

	#去掉
	my @synomous;
	tie @synomous, 'Array::Unique';
	push @synomous, @{$gene_feature->{$gene_id}->{synomous}} if(defined $gene_feature->{$gene_id}->{synomous});
	my @synomous_rm;
	foreach(@synomous){
		next if(defined $gene_synomous->{$gene_id}->{$_});
		push @synomous_rm, $_;
	}
#	my $synomous = (defined $gene_feature->{$gene_id}->{synomous}) ? join ";", @synomous : "-";
	my $synomous = ( scalar(@synomous_rm) > 1 ) ? join ";", @synomous_rm : "-";

	my @uniprot_proteinname_synomous;
	tie @uniprot_proteinname_synomous, 'Array::Unique';
	push @uniprot_proteinname_synomous, @{$gene_feature->{$gene_id}->{uniprot_proteinname_synomous}} if(defined $gene_feature->{$gene_id}->{uniprot_proteinname_synomous});
	my @uniprot_proteinname_synomous_rm;
	foreach(@uniprot_proteinname_synomous){
		next if( defined $gene_synomous->{$gene_id}->{$_});
		push @uniprot_proteinname_synomous_rm, $_;
	}
#	my $uniprot_proteinname_synomous = (defined $gene_feature->{$gene_id}->{uniprot_proteinname_synomous}) ? join ";", @uniprot_proteinname_synomous : "-";
	my $uniprot_proteinname_synomous = ( scalar(@uniprot_proteinname_synomous_rm) > 1) ? join ";", @uniprot_proteinname_synomous_rm : "-";

	
	my $print_gene_loc      = (defined $gene_feature->{$gene_id}->{gene}->{loc}) ? $gene_feature->{$gene_id}->{gene}->{loc} : "-";
	my $print_gene_organism = (defined $gene_feature->{$gene_id}->{organism})    ? $gene_feature->{$gene_id}->{organism} : "-";
	print $gene_id,               "\t";
	print $print_gene_loc, "\t";
	print $gene_seq, "\t";
	print $print_transcript,      "\t";
	print $print_transcript_stru, "\t";
	print $print_transcript_seq, "\t";
	print $print_exon_strus,      "\t";
	print $print_transcript_phase, "\t";
	print $print_exon_gc,      "\t";
	print $print_exon_purine,      "\t";
	print $print_GeneID       ,"\t";
	print $print_UniGene      ,"\t";
	print $synomous, "\t";
	print $print_refseq_id       ,"\t";
	print $print_GO       ,"\t";
	print $uniprot_accession, "\t";
	print $uniprot_proteinname_synomous,      "\t";
	print $uniprot_name, "\t";
	print $print_InterPro_id, "\t";
	print $print_Pfam_id, "\t";
	print $print_EC_id      , "\t";
	print $print_KEGG_id    , "\t";
#	print $pubmed_id, "\t";
	print $print_gene_organism, "\n";
}

################################ Main ##########################################################################
#===============================================================================================================
my $options="-dir_ensembl_fasta_dna $opts{dir_ensembl_fasta_dna} -rapdb_locus_gtf $opts{rapdb_locus_gtf} -ensembl_gtf $opts{ensembl_gtf} -uniprot_xml $opts{uniprot_xml} ";
my $endTime=localtime();
my $Program = $1 if($0 =~ /(.*)\.pl/);
open  LOG,">>$Program\_ProgramRunning.Log";
print LOG "From \<$startTime\> to \<$endTime\>\tperl $0 $options\n";
close LOG;
__END__





