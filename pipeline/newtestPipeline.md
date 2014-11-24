#文件说明:  
最新物种树species_tree_rio.Plants_complete.xml(用的code是：speclist_yangl.txt)  
~~最新的trimal过的multialign文件：01_04tree_al3.phy(对于OS05G0113900而言)~~    

#原始测试样本(OS05G0113900)： 
```
please see pipeline.md
```
#new test
##OS02T0762800-00
```java
perl ../leaf_nodes_species2pep_fa.emf.parser.pl -newick_tree /home/yangl/perl_src/OS02T0762800-00/OS02T0762800-00.genetree.newick.tree -species_list /media/yangl/TOOLS/plantortho/parser/binomial_full_abbreviation.v22.txt -pep_dir /media/yangl/TOOLS/RIO/test_ensembl_tree/plants_pep/ > OS02T0762800-00.out
perl ../filter_newick.parser.out.v2.pl OS02T0762800-00.out ../full_index out id_index
/usr/bin/mafft --genafpair --maxiterate 16 --phylipout --reorder out > accuratemulalign.phy
trimal -in accuratemulalign.phy -gt 0.8 -st 0.001 -cons 60 -phylip -out trimal_mulalign.phy
perl ../id_regain.pl trimal_mulalign.phy id_index trimal_mulalign_longname.phy
/home/yangl/bin/PhyML-3.1_linux64 -i /leofs/zhangz_group/yangl/phyml/OS02T0762800-00/trimal_mulalign_longname.phy -d aa -b 100
/*java -Xmx2048m -cp /media/yangl/TOOLS/RIO/forester_1036.jar org.forester.application.rio trimal_mulalign_longname.phy_phyml_tree.txt /media/yangl/TOOLS/RIO/species_tree_rio.Plants_complete.xml OS02T0762800-00.rio.out*/
java -Xmx2048m -cp /media/yangl/TOOLS/RIO/forester_1036.jar org.forester.application.rio trimal_mulalign_longname.phy_phyml_boot_trees.txt /media/yangl/TOOLS/RIO/species_tree_rio.Plants_complete.xml OS02T0762800-00.rio.out
```

##OS09T0133200-01
```java
perl ../leaf_nodes_species2pep_fa.emf.parser.pl -newick_tree /home/yangl/perl_src/OS09T0133200-01/OS09G0133200-01.genetree.newick.tree -species_list /media/yangl/TOOLS/plantortho/parser/binomial_full_abbreviation.v22.txt -pep_dir /media/yangl/TOOLS/RIO/test_ensembl_tree/plants_pep/ > OS09G0133200-01.out
perl ../filter_newick.parser.out.v2.pl OS09G0133200-01.out ../full_index out id_index
/usr/bin/mafft --genafpair --maxiterate 16 --phylipout --reorder out > accuratemulalign.phy
trimal -in accuratemulalign.phy -gt 0.8 -st 0.001 -cons 60 -phylip -out trimal_mulalign.phy
perl ../id_regain.pl trimal_mulalign.phy id_index trimal_mulalign_longname.phy
/home/yangl/bin/PhyML-3.1_linux64 -i /leofs/zhangz_group/yangl/phyml/OS09T0133200-01/trimal_mulalign_longname.phy -d aa -b 100
/*java -Xmx2048m -cp /media/yangl/TOOLS/RIO/forester_1036.jar org.forester.application.rio trimal_mulalign_longname.phy_phyml_tree.txt /media/yangl/TOOLS/RIO/species_tree_rio.Plants_complete.xml OS09G0133200-01.rio.out*/
java -Xmx2048m -cp /media/yangl/TOOLS/RIO/forester_1036.jar org.forester.application.rio trimal_mulalign_longname.phy_phyml_boot_trees.txt /media/yangl/TOOLS/RIO/species_tree_rio.Plants_complete.xml OS09G0133200-01.rio.out
```
