#!/bin/sh
#PBS -N phyml_run_OS09T0133200-01 
#PBS -q bioque 
#PBS -l mem=8gb,walltime=120:00:00
#HSCHED -s hschedd
cd /leofs/zhangz_group/huangdw/plantortho/pipeline/OS09G0133200-01 
/leofs/zhangz_group/huangdw/bin/PhyML-3.1/PhyML-3.1_linux64 -d aa -b 10 -m LG  -s NNI  -i/leofs/zhangz_group/huangdw/plantortho/pipeline/OS09G0133200-01/trimal_localpair.phy
/leofs/zhangz_group/huangdw/bin/PhyML-3.1/PhyML-3.1_linux64 -d aa -b 10 -m JTT -s SPR  -i /leofs/zhangz_group/huangdw/plantortho/pipeline/OS09G0133200-01/trimal_localpair.phy


/leofs/zhangz_group/huangdw/bin/PhyML-3.1/PhyML-3.1_linux64 -d aa -b 10 -i /leofs/zhangz_group/huangdw/plantortho/pipeline/OS09G0133200-01/trimal_genafpair.phy
/leofs/zhangz_group/huangdw/bin/PhyML-3.1/PhyML-3.1_linux64 -d aa -b 10 -i /leofs/zhangz_group/huangdw/plantortho/pipeline/OS09G0133200-01/trimal_auto.phy

