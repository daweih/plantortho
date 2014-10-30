#Pipeline mark
***(file path is based on my ubuntu)***
##fiter_newick.parser.out.pl
**1. description:** get normal .fa file with short name from newick.parser.pl output.     
**2. parameter:**                                                        

- perl filter_newick.parser.out.pl out.txt  out id_index 
                 
##mafft  
**1. description:**  multiple allignment.                                                            
**2. parameter:**           
                               
- /usr/bin/mafft  --auto --phylipout --reorder ../out2 > accuratemulalign01.phy *[note: --auto = FFT-NS-i]*     
- /usr/bin/mafft --lop --auto --phylipout --reorder ../out2> accuratemulalign02.phy *[note: --auto =FFT-NS-2]*
- /usr/bin/mafft  --genafpair  --maxiterate 16 --phylipout --reorder ../out2 > accuratemulalign03.phy *[note: E-INS-i]*
- /usr/bin/mafft  --globalpair --maxiterate 16 --phylipout --reorder ../out2 > accuratemulalign04.phy *[note: G-INS-i]*

##trimal
**1. description:** multiple allignment mask.  
**2. parameter:**   

***
$cat 01_04tree.phy  
../01tree/accuratemulalign01.phy   
../02tree/accuratemulalign02.phy  
../03tree/accuratemulalign03.phy  
../04tree/accuratemulalign04.phy  
***
a: trimal -compareset 01_04tree.phy -gt 0.8 -st 0.001 -cons 60 -phylip -out 01_04tree_a.phy    *[I perfer]*                     
b: trimal -compareset 01_04tree.phy -phylip -out 01_04tree_b.phy                       
c: trimal -compareset 01_04tree.phy -gt 0.8 -st 0.001  -phylip -out 01_04tree_c.phy            
d: trimal -compareset 01_04tree.phy  -phylip -out 01_04tree_d.phy -automated1       
*[note: -compareset    File Selected:	../03tree/accuratemulalign03.phy]*

##id_regain.pl
**1.description:**This program is to regain long ID for phylip format file.                    
**2. parameter:**    
 perl ../id_regain.pl 01_04tree_a.phy ../id_index 01_04tree_al.phy    
 perl ../id_regain.pl 01_04tree_b.phy ../id_index 01_04tree_bl.phy    
 perl ../id_regain.pl 01_04tree_c.phy ../id_index 01_04tree_cl.phy   
 perl ../id_regain.pl 01_04tree_d.phy ../id_index 01_04tree_dl.phy
 
##Phyml
**1.description:** phylogenetic construction.      
**2. parameter:**  
d: AA

  
 
 
           
 

 



