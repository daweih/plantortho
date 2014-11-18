‹#Pipeline mark    

##files:     
out.txt: output of newick.parser.species.pl*[the fasta file of our gene family]*   
full_index: the index file to keep entry identity of speciestree,genetree,id_index,ensembl pep file,speclist.   
species_tree_rio.Plants_complete.xml: we need add new leaf node and check code before we want to add new specie.  

##workflow:  

##fiter_newick.parser.out.pl  
**1. description:** get normal .fa file with short name from newick.parser.pl output.     
**2. parameter:**    
perl ../parser/filter_newick.parser.out.pl out.txt full_index out id_index    

##mafft   
**1. description:**  multiple allignment.                                                              
**2. parameter:**  
- /usr/bin/mafft  --genafpair  --maxiterate 16 --phylipout --reorder out > accuratemulalign.phy *[note: E-INS-i]*  

##trimal   
**1. description:** multiple allignment mask.   
**2. parameter:**   
trimal -in accuratemulalign.phy -gt 0.8 -st 0.001 -cons 60 -phylip -out 01_04tree_a.phy  
mv 01_04tree_a.phy phyml_result/    
cd phyml_result/  

##id_regain.pl  
**1.description:**This program is to regain long ID for phylip format file.                      
**2. parameter:**     
perl ../../parser/id_regain.pl 01_04tree_a.phy ../id_index3 01_04tre_al3.phy   

##Phyml   
**1.description:** phylogenetic construction.       
**2. parameter:**  
01_04tree_al.phy  d: AA  
phyml result: 01_04tree_al3.phy_phyml_tree.txt  
cd ../  

##forester_1036.jar  
**1.description:** ortholog inference.      
**2. parameter:**   
java -Xmx2048m -cp ../bin/RIO/forester_1036.jar org.forester.application.rio ./phyml_result/01_04tree_al3.phy_phyml_tree.txt species_tree_rio.Plants_complete.xml rio_out
