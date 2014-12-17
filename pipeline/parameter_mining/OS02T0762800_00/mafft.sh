#PBS -N OS02T0762800
#PBS -q bioque
#PBS -l mem=1gb,walltime=24:00:00 
#HSCHED -s hschedd
/software/biosoft/bin/mafft --genafpair --maxiterate 16 --phylipout --reorder /leofs/zhangz_group/huangdw/plantortho/pipeline/OS02T0762800_00/out > /leofs/zhangz_group/huangdw/plantortho/pipeline/OS02T0762800_00/genafpair.phy
/software/biosoft/bin/mafft --localpair --maxiterate 16 --phylipout --reorder /leofs/zhangz_group/huangdw/plantortho/pipeline/OS02T0762800_00/out > /leofs/zhangz_group/huangdw/plantortho/pipeline/OS02T0762800_00/localpair.phy
/software/biosoft/bin/mafft --auto --phylipout --reorder /leofs/zhangz_group/huangdw/plantortho/pipeline/OS02T0762800_00/out > /leofs/zhangz_group/huangdw/plantortho/pipeline/OS02T0762800_00/auto.phy

