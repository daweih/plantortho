‹#Pipeline mark    

##forester_1036.jar  

for i in `ls OS02T0762800_00/itrimal_*.b10_phyml_boot_trees.txt| awk -F'.txt' '{print $1}'`; do java -Xmx2048m -cp ../../bin/RIO/forester_1036.jar org.forester.application.rio $i.txt ../species_tree_rio.Plants_complete.xml $i.rio.txt; done >rio.log.txt
for i in `ls OS05G0113900/itrimal_*.b10_phyml_boot_trees.txt| awk -F'.txt' '{print $1}'`; do java -Xmx2048m -cp ../../bin/RIO/forester_1036.jar org.forester.application.rio $i.txt ../species_tree_rio.Plants_complete.xml $i.rio.txt; done >>rio.log.txt
for i in `ls OS09G0133200-01/itrimal_*.b10_phyml_boot_trees.txt| awk -F'.txt' '{print $1}'`; do java -Xmx2048m -cp ../../bin/RIO/forester_1036.jar org.forester.application.rio $i.txt ../species_tree_rio.Plants_complete.xml $i.rio.txt; done >>rio.log.txt
 
## select ortholog record
for i in `ls ../pipeline/parameter_mining/OS02T0762800_00/*.rio.txt`; do  perl ortho_tree_seq2table.pl  -e_tree ../pipeline/OS02T0762800-00/OS02T0762800-00.genetree.newick.tree -e_ortholog ../pipeline/OS02T0762800-00/orthologues-ComparaOrthologs-Oryza_sativa-Gene-Compara_Ortholog-76-OS02G0762800.csv -rio_out     $i; done
for i in `ls ../pipeline/parameter_mining/OS05G0113900/*.rio.txt`; do  perl ortho_tree_seq2table.pl  -e_tree ../pipeline/OS05G0113900.genetree.newick.tree -e_ortholog ../pipeline/orthologues-ComparaOrthologs-Oryza_sativa-Gene-Compara_Ortholog-76-OS05G0113900.csv -rio_out     $i; done
for i in `ls ../pipeline/parameter_mining/OS09G0133200-01/*.rio.txt`; do  perl ortho_tree_seq2table.pl -e_tree ../pipeline/OS09T0133200-01/OS09G0133200-01.genetree.newick.tree -e_ortholog ../pipeline/OS09T0133200-01/orthologues-ComparaOrthologs-Oryza_sativa-Gene-Compara_Ortholog-76-OS09G0133200.csv -rio_out     $i; done


 -rio_out ../pipeline/OS05G0113900/rio_mafft_tri30_out 