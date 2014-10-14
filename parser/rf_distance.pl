# Contact     : daweih@me.com
# Date        : Tue Oct 14 21:53:02 CST 2014
# Last Update : 
# Reference   : 
# Description : perl rf_distance.pl ~/plantortho/pipeline/phyml_result/accuratemulalign3_longname.phy_phyml_boot_trees.txt ~/plantortho/pipeline/phyml_result/accuratemulalign3_longname.phy_phyml_tree.txt

#===============================================================================================================
use lib "/opt/perl5.12.3/lib/site_perl/5.12.3";
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Bio::Phylo::IO;


################################ Main ##########################################################################
my $tree_list;
open I, "< $ARGV[0]";
while(<I>){
	chomp;
	$tree_list->{$_} = 1;
}
print "# to bootstrap $ARGV[0]:\n";
foreach my $tree_a (keys %{$tree_list}){
	my $forest_a = Bio::Phylo::IO->parse(-format => 'newick', -string => $tree_a);
	my $tree_x = $forest_a->first;
	foreach my $tree_b (keys %{$tree_list}){
		my $forest_b = Bio::Phylo::IO->parse(-format => 'newick', -string => $tree_b);
		my $tree_y = $forest_b->first;
		my $symdiff = $tree_x->calc_symdiff($tree_y);
		print $symdiff, "\t";
	}
	print "\n";
}

print "\n# to tree $ARGV[1]:\n";
open I, "< $ARGV[1]";
while(<I>){
	chomp;
	my $forest_a = Bio::Phylo::IO->parse(-format => 'newick', -string => $_);
	my $tree_x = $forest_a->first;
	foreach my $tree_b (keys %{$tree_list}){
		my $forest_b = Bio::Phylo::IO->parse(-format => 'newick', -string => $tree_b);
		my $tree_y = $forest_b->first;
		my $symdiff = $tree_x->calc_symdiff($tree_y);
		print $symdiff, "\t";
	}
	print "\n";
}

__END__
