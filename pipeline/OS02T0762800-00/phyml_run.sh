#!/bin/sh
#PBS -N phyml_run_OS02T0762800-00 
#PBS -q bioque 
#PBS -l mem=8gb,walltime=120:00:00
#HSCHED -s hschedd
cd /leofs/zhangz_group/yangl/phyml/OS02T0762800-00
/home/yangl/bin/PhyML-3.1_linux64 -i /leofs/zhangz_group/yangl/phyml/OS02T0762800-00/trimal_mulalign_longname.phy -d aa -b 100
