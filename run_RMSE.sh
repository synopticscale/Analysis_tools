#!/bin/bash
# ###################
# PBS -N ens_FSS 
# PBS -o /glade/campaign/univ/umcp0011/colab_scripts/run_fss.log
# PBS -e /glade/campaign/univ/umcp0011/colab_scripts/run_fss.err
# PBS -l select=1:ncpus=1
# PBS -l walltime=06:00:00
# PBS -q casper 
# PBS -A UMCP0011
# PBS -J 1-15
# ####################
rm -rf process_fss_*.log
module load python
module load ncarenv
ncar_pylib
################
#settings 
export scratch=/glade/scratch/jmccurry/WOF/realtime
#export scratch=/glade/campaign/univ/umcp0011/OSSE_CYCLE_ARCHIVE/storage_20230711
#export scratch=/glade/scratch/jmccurry/WRF-DART/DATA/DATA/TEST_20230514

export campaign_nature_store=/glade/campaign/univ/umcp0011/OSSE_NATURE_ARCHIVE/natureruns
variables='PREC_ACC_NC,U10,V10,T2,U_500,U_850,T_500,T_850'
#variables='U10,V10,T2,U_500,U_850,T_500,T_850'

#variables='PREC_ACC_NC,T2'
#events='20200413 20200812 20200903 20210717 20220608 20220702 20220716'
#events='20200812'
#events='20210717'
events='20220716'
cases='default case1 case2 case3 case4'
#cases='NO_BIAS_case1 enkf_NO_BIAS'
#case='default_hum'
#cases='default case1'

################
cd /glade/campaign/univ/umcp0011/colab_scripts
for event in $events;do
    for case in $cases;do
            echo $event.$case
             for i in {0..20};do
            

                    python process_RMSE.py $event $i $scratch $variables $case

done
done
done

#for event in $events;do
#for case in $cases;do
#echo $event.$case
#python process_RMSE_cycle.py $event $case $variables
#done 
#done
