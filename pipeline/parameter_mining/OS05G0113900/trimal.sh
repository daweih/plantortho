#PBS -q bioque
#PBS -l mem=1gb,walltime=24:00:00 
#HSCHED -s hschedd
#OK ，但是你的trimal 参数用的是我以前的版本，按我的理解用这个参数会好一点： trimal -in multialignFile.phy -gt 0.2 -out multialignMasked.txt   上次你测试后结果会好一点点
/leofs/zhangz_group/huangdw/bin/trimAl/source/trimal -in /leofs/zhangz_group/huangdw/plantortho/pipeline/OS05G0113900/genafpair.phy -gt 0.8 -st 0.001 -cons 60 -phylip -out /leofs/zhangz_group/huangdw/plantortho/pipeline/OS05G0113900/trimal_genafpair.phy
/leofs/zhangz_group/huangdw/bin/trimAl/source/trimal -in /leofs/zhangz_group/huangdw/plantortho/pipeline/OS05G0113900/localpair.phy -gt 0.8 -st 0.001 -cons 60 -phylip -out /leofs/zhangz_group/huangdw/plantortho/pipeline/OS05G0113900/trimal_localpair.phy
/leofs/zhangz_group/huangdw/bin/trimAl/source/trimal -in /leofs/zhangz_group/huangdw/plantortho/pipeline/OS05G0113900/auto.phy      -gt 0.8 -st 0.001 -cons 60 -phylip -out /leofs/zhangz_group/huangdw/plantortho/pipeline/OS05G0113900/trimal_auto.phy     

