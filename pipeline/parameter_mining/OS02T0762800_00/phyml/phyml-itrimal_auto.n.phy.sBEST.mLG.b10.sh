#!/bin/sh
#PBS -N phyml-itrimal_auto.n.phy.sBEST.mLG.b10.sh
#PBS -q bioque
#PBS -l mem=8gb,walltime=120:00:00
#HSCHED -s hschedd
cd /leofs/zhangz_group/huangdw/plantortho/pipeline/OS09G0133200-01
/leofs/zhangz_group/huangdw/bin/PhyML-3.1/PhyML-3.1_linux64 -d aa -b 10 -m LG  -s BEST  -i /leofs/zhangz_group/huangdw/plantortho/pipeline/OS02T0762800_00/itrimal_auto.n.phy.sBEST.mLG.b10
