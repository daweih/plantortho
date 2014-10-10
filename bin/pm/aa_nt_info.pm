# Contact     : daweimhuang@gmail.com
# Date        : Fri Oct 10 14:30:14 CST 2014
# Last Update : 
# Reference   : 
# Description : 
# Require     :

package aa_nt_info;
use strict;
use warnings;
require Exporter;
use Data::Dumper;

our @ISA = qw(Exporter);
our @EXPORT = qw(aa_nt_info);
our $VERSION = '0.01';


my $amino_acids;
$amino_acids->{"*"} = {"abb" => "*",   "full" => "STOP",          "polarity" => "NA"          ,"atomic_mass" => "NA",};
$amino_acids->{"A"} = {"abb" => "Ala", "full" => "Alanine",       "polarity" => "nonpolar"    ,"atomic_mass" => "89.09"};
$amino_acids->{"C"} = {"abb" => "Cys", "full" => "Cysteine",      "polarity" => "nonpolar"    ,"atomic_mass" => "121.16"};
$amino_acids->{"D"} = {"abb" => "Asp", "full" => "Aspartic",      "polarity" => "acidic polar","atomic_mass" => "133.10"};
$amino_acids->{"E"} = {"abb" => "Glu", "full" => "Glutamic",      "polarity" => "acidic polar","atomic_mass" => "147.13"};
$amino_acids->{"F"} = {"abb" => "Phe", "full" => "Phenylalanine", "polarity" => "nonpolar"    ,"atomic_mass" => "165.19"};
$amino_acids->{"G"} = {"abb" => "Gly", "full" => "Glycine",       "polarity" => "nonpolar"    ,"atomic_mass" => "75.07"};
$amino_acids->{"H"} = {"abb" => "His", "full" => "Histidine",     "polarity" => "Basic polar" ,"atomic_mass" => "155.16"};
$amino_acids->{"I"} = {"abb" => "Ile", "full" => "Isoleucine",    "polarity" => "nonpolar"    ,"atomic_mass" => "131.17"};
$amino_acids->{"K"} = {"abb" => "Lys", "full" => "Lysine",        "polarity" => "Basic polar" ,"atomic_mass" => "146.19"};
$amino_acids->{"L"} = {"abb" => "Leu", "full" => "Leucine",       "polarity" => "nonpolar"    ,"atomic_mass" => "131.17"};
$amino_acids->{"M"} = {"abb" => "Met", "full" => "Methionine",    "polarity" => "nonpolar"    ,"atomic_mass" => "149.21"};
$amino_acids->{"N"} = {"abb" => "Asn", "full" => "Asparagine",    "polarity" => "polar"       ,"atomic_mass" => "132.12"};
$amino_acids->{"P"} = {"abb" => "Pro", "full" => "Proline",       "polarity" => "nonpolar"    ,"atomic_mass" => "115.13"};
$amino_acids->{"Q"} = {"abb" => "Gln", "full" => "Glutamine",     "polarity" => "polar"       ,"atomic_mass" => "146.15"};
$amino_acids->{"R"} = {"abb" => "Arg", "full" => "Arginine",      "polarity" => "Basic polar" ,"atomic_mass" => "174.20"};
$amino_acids->{"S"} = {"abb" => "Ser", "full" => "Serine",        "polarity" => "polar"       ,"atomic_mass" => "105.09"}; 
$amino_acids->{"T"} = {"abb" => "Thr", "full" => "Threonine",     "polarity" => "polar"       ,"atomic_mass" => "119.12"};#Threonine	Thr	T	polar	neutral
$amino_acids->{"V"} = {"abb" => "Val", "full" => "Valine",        "polarity" => "nonpolar"    ,"atomic_mass" => "117.15"};
$amino_acids->{"W"} = {"abb" => "Trp", "full" => "Tryptophan",    "polarity" => "nonpolar"    ,"atomic_mass" => "204.23"};
$amino_acids->{"Y"} = {"abb" => "Tyr", "full" => "Tyrosine",      "polarity" => "polar"       ,"atomic_mass" => "181.19"};#Tyrosine	Tyr	Y	polar
foreach my $symbol (keys %{$amino_acids}){
	$amino_acids->{$symbol}->{"symbol"} = $symbol;
	$amino_acids->{lc($amino_acids->{$symbol}->{abb})} = $amino_acids->{$symbol};
	$amino_acids->{lc($amino_acids->{$symbol}->{full})} = $amino_acids->{$symbol};
}
my $nucleic_acids = {
"a" => "Adenosine",
"c" => "Cytidine",
"g" => "Guanosine",
"t" => "Thymidine"
};
foreach my $symbol (keys %{$nucleic_acids}){
	$nucleic_acids->{lc($nucleic_acids->{$symbol})} = $symbol;
}
my $nucleobase = {
"r" => "Purine",
"n" => "Purine and Pyrimidine",
"y" => "Pyrimidine"
};
foreach my $symbol (keys %{$nucleobase}){
	$nucleobase->{lc($nucleobase->{$symbol})} = $symbol;
}


sub aa_nt_info{
	my ($str, $type) = @_;
	if($type eq "nb"){
		return $nucleobase->{lc($str)};
	}
	elsif($type eq "na"){
		return $nucleic_acids->{lc($str)};
	}
	elsif($type eq "aa"){
		return $amino_acids->{lc($str)};
	}
}
1;
__END__
my $test =  &aa_nt_info("y", "nb") ;
   $test =  &aa_nt_info("A", "na") ;
   $test =  &aa_nt_info("*", "aa") ;
print Dumper($test);
