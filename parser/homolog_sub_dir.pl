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
use Bio::SeqIO;

my %opts;
GetOptions(\%opts,'id:s', 'pep_dir:s');
my $usage= <<"USAGE";
	Program: $0
	INPUT:
		-option			some option content
	OUTPUT:
	
USAGE

#die $usage unless ($opts{id});
my $startTime=localtime();
################################ Main ##########################################################################
my $pwd = `pwd`;
$pwd =~ s/\/$//;
chomp($pwd);

my $rec;
my $index = 0;
my $check = 0;
my $check_species_uniqc = 1;
my $included = 0;
my @check_species_uniqc = qw / oryza_barthii oryza_brachyantha oryza_glaberrima oryza_glumaepatula oryza_indica oryza_meridionalis oryza_nivara oryza_punctata oryza_sativa arabidopsis_thaliana /;
my @check_species = qw /  arabidopsis_thaliana brassica_rapa glycine_max oryza_barthii oryza_brachyantha oryza_glaberrima oryza_glumaepatula oryza_meridionalis oryza_nivara oryza_punctata oryza_sativa populus_trichocarpa solanum_lycopersicum sorghum_bicolor zea_mays  /;

my @species17 = qw / arabidopsis_thaliana brassica_rapa glycine_max oryza_barthii oryza_brachyantha oryza_glaberrima oryza_glumaepatula oryza_indica oryza_meridionalis oryza_nivara oryza_punctata oryza_rufipogon oryza_sativa populus_trichocarpa solanum_lycopersicum sorghum_bicolor zea_mays /;
my $species;
foreach(@species17){
	$species->{$_} = 1;
}

open I, "< Compara.newick_trees.24.emf";
my $emf;
my $accession_list;
while(<I>){
	$emf .= $_;
	if( $_ eq "//\n" ){
		my @rec = split /\n/,$emf;
		foreach my $rec (@rec){
			next if($rec !~ /^SEQ.*/);
			#SEQ oryza_barthii OBART12G20230.1 12 20479112 20479483 -1 OBART12G20230 NULL
			my @column = split /\s+/,$rec;
			next if(!defined $species->{$column[1]});
			$accession_list->{$column[1]}->{$column[2]}->{gene} = $column[7];
		}
		$emf = "";
	}
}

$opts{pep_dir} = "/leofs/zhangz_group/huangdw/plantortho/release24_pep_s15/";
opendir PEP , $opts{pep_dir} or die "Cannot open plants_pep: $!";
foreach my $file (readdir PEP)
{
	next if($file !~ /^.*.pep.all.fa$/);
	my $check_file_species = 0;
	my $species_name;
	foreach(@species17){
		my $file_speciesname = uc($1). $2 if($_ =~ /^(\S)(.*)$/);
		#Aegilops_tauschii.GCA_000347335.1.22.pep.all.fa
		#Atau
		if($file =~ /^$file_speciesname\..*$/ ){
			$check_file_species = 1;
			$species_name = $_;
			print $file, "\t$species_name\n";
		}
	}
	if($check_file_species){
#		print $file, "\n";
		my $pep = Bio::SeqIO->new(-file => $opts{pep_dir}. "$file" ,
										-format => 'Fasta');
#										print "$pwd/../release24_pep_s15/$file\n";
		while ( my $seq = $pep->next_seq() ) {
			next if(!defined $accession_list->{$species_name}->{$seq->id});
			$accession_list->{$species_name}->{$seq->id}->{seq} = $seq->seq;
		}
	}
}
closedir PEP;
#foreach my $species (keys %{$accession_list}){print $species, "\t", scalar(keys %{$accession_list->{$species}}), "\n";}


my $oryza_sativa_id;
my $sub_dir = 0;
my $sub_dir_index = 0;
open I, "< Compara.newick_trees.24.emf";
while(<I>){
	$rec .= $_;
	if( $_ eq "//\n" ){
		my @rec = split /\n/,$rec;
		my $rec_species_uniqc;
		my $newick;
		foreach my $rec (@rec){
if($rec =~ /^\(.*/){
$newick = $rec;
}
			next if($rec !~ /^SEQ.*/);
			my @column = split /\s+/,$rec;
			if( $column[1] eq "oryza_sativa" ){
				$check = 1;
				$oryza_sativa_id = $column[2];
			}
			$rec_species_uniqc->{$column[1]} += 1;
		}

		foreach my $species ( @check_species_uniqc ){
			next if( !defined $rec_species_uniqc->{$species} );
			$check_species_uniqc = 0 if( $rec_species_uniqc->{$species} > 1 );
		}

		foreach my $species ( @check_species ){
			next if( !defined $rec_species_uniqc->{$species} );
			$included += 1;
		}
		
		if( $check == 1 && $included > 2 ){
			if($sub_dir_index==0){
				system "mkdir sub_$sub_dir";
			}
			system "mkdir sub_$sub_dir/$oryza_sativa_id";

			open  O, "> sub_$sub_dir/". $oryza_sativa_id. "/out.emf.txt";
			print O $rec;
			close O;

			open  O, "> sub_$sub_dir/". $oryza_sativa_id. "/out.newick.txt";
			print O  $newick;
			close O;

			open  O, "> sub_$sub_dir/". $oryza_sativa_id. "/out.fasta";
			open  OI, "> sub_$sub_dir/". $oryza_sativa_id. "/out.index.txt";
			my $accession_cnt = 0;
			foreach my $rec (@rec){
				next if($rec !~ /^SEQ.*/);
				my @column = split /\s+/,$rec;
				#SEQ oryza_barthii OBART12G20230.1 12 20479112 20479483 -1 OBART12G20230 NULL

				if(defined $species->{$column[1]}){
					print O ">$accession_cnt\t", "$column[2]\t$column[1]\n", $accession_list->{$column[1]}->{$column[2]}->{seq},"\n";
					print OI "$accession_cnt\t$column[2]\n";
				}
				else{
				}
				$accession_cnt+=1;
			}
			close O;
			close OI;


			open O, "> sub_$sub_dir/". $oryza_sativa_id. "/shell.sh";
			print O "#PBS -N $oryza_sativa_id\n#PBS -q bioque\n#PBS -l mem=1gb,walltime=24:00:00\n#HSCHED -s hschedd\n";

			print O "/software/biosoft/bin/mafft --genafpair --maxiterate 16 --phylipout --reorder $pwd/sub_$sub_dir/$oryza_sativa_id/out.fasta >$pwd/sub_$sub_dir/$oryza_sativa_id/gfp.phy\n";

			print O "/leofs/zhangz_group/huangdw/bin/trimAl/source/trimal -in $pwd/sub_$sub_dir/$oryza_sativa_id/gfp.phy -gt 0.8 -st 0.001 -cons 60 -phylip -out $pwd/sub_$sub_dir/$oryza_sativa_id/tri_gfp.phy\n";
			print O "/leofs/zhangz_group/huangdw/bin/trimAl/source/trimal -in $pwd/sub_$sub_dir/$oryza_sativa_id/gfp.phy -gt 0.2 -phylip                    -out $pwd/sub_$sub_dir/$oryza_sativa_id/tri_gt02_gfp.phy\n";

			print O "cp $pwd/sub_$sub_dir/$oryza_sativa_id/tri_gfp.phy      $pwd/sub_$sub_dir/$oryza_sativa_id/tri_gfp.b10.jtt.nni.phy\n";
			print O "cp $pwd/sub_$sub_dir/$oryza_sativa_id/tri_gfp.phy      $pwd/sub_$sub_dir/$oryza_sativa_id/tri_gfp.b10.jtt.spr.phy\n";
			print O "cp $pwd/sub_$sub_dir/$oryza_sativa_id/tri_gfp.phy      $pwd/sub_$sub_dir/$oryza_sativa_id/tri_gfp.b10.jtt.bst.phy\n";
			print O "cp $pwd/sub_$sub_dir/$oryza_sativa_id/tri_gfp.phy      $pwd/sub_$sub_dir/$oryza_sativa_id/tri_gfp.b10.lg..nni.phy\n";
			print O "cp $pwd/sub_$sub_dir/$oryza_sativa_id/tri_gfp.phy      $pwd/sub_$sub_dir/$oryza_sativa_id/tri_gfp.b10.lg..spr.phy\n";
			print O "cp $pwd/sub_$sub_dir/$oryza_sativa_id/tri_gfp.phy      $pwd/sub_$sub_dir/$oryza_sativa_id/tri_gfp.b10.lg..bst.phy\n";

			print O "cp $pwd/sub_$sub_dir/$oryza_sativa_id/tri_gt02_gfp.phy $pwd/sub_$sub_dir/$oryza_sativa_id/tri_gt02_gfp.b10.jtt.nni.phy\n";
			print O "cp $pwd/sub_$sub_dir/$oryza_sativa_id/tri_gt02_gfp.phy $pwd/sub_$sub_dir/$oryza_sativa_id/tri_gt02_gfp.b10.jtt.spr.phy\n";
			print O "cp $pwd/sub_$sub_dir/$oryza_sativa_id/tri_gt02_gfp.phy $pwd/sub_$sub_dir/$oryza_sativa_id/tri_gt02_gfp.b10.jtt.bst.phy\n";
			print O "cp $pwd/sub_$sub_dir/$oryza_sativa_id/tri_gt02_gfp.phy $pwd/sub_$sub_dir/$oryza_sativa_id/tri_gt02_gfp.b10.lg..nni.phy\n";
			print O "cp $pwd/sub_$sub_dir/$oryza_sativa_id/tri_gt02_gfp.phy $pwd/sub_$sub_dir/$oryza_sativa_id/tri_gt02_gfp.b10.lg..spr.phy\n";
			print O "cp $pwd/sub_$sub_dir/$oryza_sativa_id/tri_gt02_gfp.phy $pwd/sub_$sub_dir/$oryza_sativa_id/tri_gt02_gfp.b10.lg..bst.phy\n";

			print O "/leofs/zhangz_group/huangdw/bin/PhyML-3.1/PhyML-3.1_linux64 -d aa -b 10 -m JTT  -s NNI  -i $pwd/sub_$sub_dir/$oryza_sativa_id/tri_gfp.b10.jtt.nni.phy\n";
			print O "/leofs/zhangz_group/huangdw/bin/PhyML-3.1/PhyML-3.1_linux64 -d aa -b 10 -m JTT  -s SPR  -i $pwd/sub_$sub_dir/$oryza_sativa_id/tri_gfp.b10.jtt.spr.phy\n";
			print O "/leofs/zhangz_group/huangdw/bin/PhyML-3.1/PhyML-3.1_linux64 -d aa -b 10 -m JTT  -s BEST -i $pwd/sub_$sub_dir/$oryza_sativa_id/tri_gfp.b10.jtt.bst.phy\n";
			print O "/leofs/zhangz_group/huangdw/bin/PhyML-3.1/PhyML-3.1_linux64 -d aa -b 10 -m LG   -s NNI  -i $pwd/sub_$sub_dir/$oryza_sativa_id/tri_gfp.b10.lg..nni.phy\n";
			print O "/leofs/zhangz_group/huangdw/bin/PhyML-3.1/PhyML-3.1_linux64 -d aa -b 10 -m LG   -s SPR  -i $pwd/sub_$sub_dir/$oryza_sativa_id/tri_gfp.b10.lg..spr.phy\n";
			print O "/leofs/zhangz_group/huangdw/bin/PhyML-3.1/PhyML-3.1_linux64 -d aa -b 10 -m LG   -s BEST -i $pwd/sub_$sub_dir/$oryza_sativa_id/tri_gfp.b10.lg..bst.phy\n";

			print O "/leofs/zhangz_group/huangdw/bin/PhyML-3.1/PhyML-3.1_linux64 -d aa -b 10 -m JTT  -s NNI  -i $pwd/sub_$sub_dir/$oryza_sativa_id/tri_gt02_gfp.b10.jtt.nni.phy\n";
			print O "/leofs/zhangz_group/huangdw/bin/PhyML-3.1/PhyML-3.1_linux64 -d aa -b 10 -m JTT  -s SPR  -i $pwd/sub_$sub_dir/$oryza_sativa_id/tri_gt02_gfp.b10.jtt.spr.phy\n";
			print O "/leofs/zhangz_group/huangdw/bin/PhyML-3.1/PhyML-3.1_linux64 -d aa -b 10 -m JTT  -s BEST -i $pwd/sub_$sub_dir/$oryza_sativa_id/tri_gt02_gfp.b10.jtt.bst.phy\n";
			print O "/leofs/zhangz_group/huangdw/bin/PhyML-3.1/PhyML-3.1_linux64 -d aa -b 10 -m LG   -s NNI  -i $pwd/sub_$sub_dir/$oryza_sativa_id/tri_gt02_gfp.b10.lg..nni.phy\n";
			print O "/leofs/zhangz_group/huangdw/bin/PhyML-3.1/PhyML-3.1_linux64 -d aa -b 10 -m LG   -s SPR  -i $pwd/sub_$sub_dir/$oryza_sativa_id/tri_gt02_gfp.b10.lg..spr.phy\n";
			print O "/leofs/zhangz_group/huangdw/bin/PhyML-3.1/PhyML-3.1_linux64 -d aa -b 10 -m LG   -s BEST -i $pwd/sub_$sub_dir/$oryza_sativa_id/tri_gt02_gfp.b10.lg..bst.phy\n\n";

			print O "cd $pwd\n";
			print O "perl /leofs/zhangz_group/huangdw/bin/id_index.pl -dir sub_$sub_dir/$oryza_sativa_id\n";


			print O "/usr/bin/java  -Xmx2048m -cp /leofs/zhangz_group/huangdw/bin/forester_1036.jar org.forester.application.rio $pwd/sub_$sub_dir/$oryza_sativa_id/tri_gfp.b10.jtt.nni.phy_phyml_boot_trees.c.txt /leofs/zhangz_group/huangdw/bin/species_tree_rio.Plants_complete.xml $pwd/sub_$sub_dir/$oryza_sativa_id/tri_gfp.b10.jtt.nni.phy_phyml_boot_trees.c.rio.txt\n";
			print O "/usr/bin/java  -Xmx2048m -cp /leofs/zhangz_group/huangdw/bin/forester_1036.jar org.forester.application.rio $pwd/sub_$sub_dir/$oryza_sativa_id/tri_gfp.b10.jtt.spr.phy_phyml_boot_trees.c.txt /leofs/zhangz_group/huangdw/bin/species_tree_rio.Plants_complete.xml $pwd/sub_$sub_dir/$oryza_sativa_id/tri_gfp.b10.jtt.spr.phy_phyml_boot_trees.c.rio.txt\n";
			print O "/usr/bin/java  -Xmx2048m -cp /leofs/zhangz_group/huangdw/bin/forester_1036.jar org.forester.application.rio $pwd/sub_$sub_dir/$oryza_sativa_id/tri_gfp.b10.jtt.bst.phy_phyml_boot_trees.c.txt /leofs/zhangz_group/huangdw/bin/species_tree_rio.Plants_complete.xml $pwd/sub_$sub_dir/$oryza_sativa_id/tri_gfp.b10.jtt.bst.phy_phyml_boot_trees.c.rio.txt\n";
			print O "/usr/bin/java  -Xmx2048m -cp /leofs/zhangz_group/huangdw/bin/forester_1036.jar org.forester.application.rio $pwd/sub_$sub_dir/$oryza_sativa_id/tri_gfp.b10.lg..nni.phy_phyml_boot_trees.c.txt /leofs/zhangz_group/huangdw/bin/species_tree_rio.Plants_complete.xml $pwd/sub_$sub_dir/$oryza_sativa_id/tri_gfp.b10.lg..nni.phy_phyml_boot_trees.c.rio.txt\n";
			print O "/usr/bin/java  -Xmx2048m -cp /leofs/zhangz_group/huangdw/bin/forester_1036.jar org.forester.application.rio $pwd/sub_$sub_dir/$oryza_sativa_id/tri_gfp.b10.lg..spr.phy_phyml_boot_trees.c.txt /leofs/zhangz_group/huangdw/bin/species_tree_rio.Plants_complete.xml $pwd/sub_$sub_dir/$oryza_sativa_id/tri_gfp.b10.lg..spr.phy_phyml_boot_trees.c.rio.txt\n";
			print O "/usr/bin/java  -Xmx2048m -cp /leofs/zhangz_group/huangdw/bin/forester_1036.jar org.forester.application.rio $pwd/sub_$sub_dir/$oryza_sativa_id/tri_gfp.b10.lg..bst.phy_phyml_boot_trees.c.txt /leofs/zhangz_group/huangdw/bin/species_tree_rio.Plants_complete.xml $pwd/sub_$sub_dir/$oryza_sativa_id/tri_gfp.b10.lg..bst.phy_phyml_boot_trees.c.rio.txt\n";

			print O "/usr/bin/java  -Xmx2048m -cp /leofs/zhangz_group/huangdw/bin/forester_1036.jar org.forester.application.rio $pwd/sub_$sub_dir/$oryza_sativa_id/tri_gt02_gfp.b10.jtt.nni.phy_phyml_boot_trees.c.txt /leofs/zhangz_group/huangdw/bin/species_tree_rio.Plants_complete.xml $pwd/sub_$sub_dir/$oryza_sativa_id/tri_gt02_gfp.b10.jtt.nni.phy_phyml_boot_trees.c.rio.txt\n";
			print O "/usr/bin/java  -Xmx2048m -cp /leofs/zhangz_group/huangdw/bin/forester_1036.jar org.forester.application.rio $pwd/sub_$sub_dir/$oryza_sativa_id/tri_gt02_gfp.b10.jtt.spr.phy_phyml_boot_trees.c.txt /leofs/zhangz_group/huangdw/bin/species_tree_rio.Plants_complete.xml $pwd/sub_$sub_dir/$oryza_sativa_id/tri_gt02_gfp.b10.jtt.spr.phy_phyml_boot_trees.c.rio.txt\n";
			print O "/usr/bin/java  -Xmx2048m -cp /leofs/zhangz_group/huangdw/bin/forester_1036.jar org.forester.application.rio $pwd/sub_$sub_dir/$oryza_sativa_id/tri_gt02_gfp.b10.jtt.bst.phy_phyml_boot_trees.c.txt /leofs/zhangz_group/huangdw/bin/species_tree_rio.Plants_complete.xml $pwd/sub_$sub_dir/$oryza_sativa_id/tri_gt02_gfp.b10.jtt.bst.phy_phyml_boot_trees.c.rio.txt\n";
			print O "/usr/bin/java  -Xmx2048m -cp /leofs/zhangz_group/huangdw/bin/forester_1036.jar org.forester.application.rio $pwd/sub_$sub_dir/$oryza_sativa_id/tri_gt02_gfp.b10.lg..nni.phy_phyml_boot_trees.c.txt /leofs/zhangz_group/huangdw/bin/species_tree_rio.Plants_complete.xml $pwd/sub_$sub_dir/$oryza_sativa_id/tri_gt02_gfp.b10.lg..nni.phy_phyml_boot_trees.c.rio.txt\n";
			print O "/usr/bin/java  -Xmx2048m -cp /leofs/zhangz_group/huangdw/bin/forester_1036.jar org.forester.application.rio $pwd/sub_$sub_dir/$oryza_sativa_id/tri_gt02_gfp.b10.lg..spr.phy_phyml_boot_trees.c.txt /leofs/zhangz_group/huangdw/bin/species_tree_rio.Plants_complete.xml $pwd/sub_$sub_dir/$oryza_sativa_id/tri_gt02_gfp.b10.lg..spr.phy_phyml_boot_trees.c.rio.txt\n";
			print O "/usr/bin/java  -Xmx2048m -cp /leofs/zhangz_group/huangdw/bin/forester_1036.jar org.forester.application.rio $pwd/sub_$sub_dir/$oryza_sativa_id/tri_gt02_gfp.b10.lg..bst.phy_phyml_boot_trees.c.txt /leofs/zhangz_group/huangdw/bin/species_tree_rio.Plants_complete.xml $pwd/sub_$sub_dir/$oryza_sativa_id/tri_gt02_gfp.b10.lg..bst.phy_phyml_boot_trees.c.rio.txt\n";


			close O;




			$sub_dir_index += 1;
			if($sub_dir_index==100){
				$sub_dir_index = 0;
				$sub_dir += 1;
			}
		}
		$rec = "";		
		$index+=1;
		$check = 0;
		$included = 0;
		$check_species_uniqc = 1;
	}
}

__END__

#		print $rec;
#		foreach my $key (keys %{$rec_species_uniqc}){
#			$rec_species_uniqc->{$key} = 0;
#			print $key, $rec_species_uniqc->{$key}, "\n";
#		}
		if( $check == 1 && $check_species_uniqc == 0 ){
			print $rec, "\n";
		}



oryza_barthii
oryza_brachyantha
oryza_glaberrima
oryza_glumaepatula
oryza_indica
oryza_meridionalis
oryza_nivara
oryza_punctata
oryza_rufipogon
oryza_sativa


#ortholog
wget "http://ensembl.gramene.org/Oryza_sativa/Gene/Compara_Ortholog/pan_compara?db=core;g=OS02T0762800-00;_format=Excel" -O OS02T0762800-00.csv

for i in `ls| awk -F'/' '/OS/{print $1}'`; do wget "http://ensembl.gramene.org/Oryza_sativa/Gene/Compara_Ortholog/pan_compara?db=core;g=$i;_format=Excel" -O $i/$i.csv; sleep 2; done

for i in `ls| awk -F'/' '/OS/{print $1}'| awk -F'-' '{print $1}'`; do wget "http://ensembl.gramene.org/Oryza_sativa/Gene/Compara_Ortholog/pan_compara?db=core;g=$i;_format=Excel" -O $i/$i.csv; sleep 2; done
http://ensembl.gramene.org/Oryza_sativa/Gene/Compara_Ortholog/pan_compara?db=core;g=OS01T0133200;_format=Excel

#tree:
http://ensembl.gramene.org/Oryza_sativa/Export/Output/Gene?db=core;flank3_display=0;flank5_display=0;g=OS03G0168700;output=phyloxml;t=OS03T0168700-00;_format=Text
oryza_barthii
oryza_brachyantha
oryza_glaberrima
oryza_glumaepatula
oryza_meridionalis
oryza_nivara
oryza_punctata
oryza_sativa
zea_mays
arabidopsis_thaliana
brassica_rapa
glycine_max
solanum_lycopersicum
sorghum_bicolor
populus_trichocarpa
http://ensembl.gramene.org/Oryza_sativa/Gene/Compara_Ortholog?db=core;g=OS02G0762800;r=2:32125904-32134832;t=OS02T0762800-00;filename=ComparaOrthologs-Oryza_sativa-Gene-Compara_Ortholog-77-OS02G0762800;_format=Excel


http://ensembl.gramene.org/Oryza_sativa/Gene/Compara_Ortholog/pan_compara?db=core;g=OS11G0181200;_format=Excel

OS12G0618000;_format=Excel


OS11T0181200-01/OS11T0181200-01.csv
