#!/bin/bash
####################
#PBS -N ens_FSS 
#PBS -o /glade/work/jmccurry/ANALYSIS/batch_scripts/run_fss6.log
#PBS -e /glade/work/jmccurry/ANALYSIS/batch_scripts/run_fss6.err
#PBS -l select=1:ncpus=1
#PBS -l walltime=06:00:00
#PBS -q regular
#PBS -k oed
#PBS -A UMCP0011
#PBS -J 1-20 
#####################

#NOTES:
#changed on 05/03/2023 to grid nature obs online instead of using gridded obs sequence files 
#'process_fss_refl_ensemble' and 'process_fss_real_precip_ensemble' still work on old system
#'process_fss_refl_ensemble' uses 'cropped_osse' obs_type for nature reflectivity
#'process_fss_real_precip_ensemble' uses ????? re-gridded MRMS precip obs 
#new script 'process_fss_ensemble' uses obs_type to refer to WRF field of interest 

module load python
module load ncarenv
module load cdo 
ncar_pylib
################
#settings 
mems=20
scratch=/glade/scratch/jmccurry/WRF-DART/DATA/DATA/TEST_20230514
#cases='default case2 case1 case3 case4'
cases='case2'
#cases='case1'

#cases='case1_newfcst default_newfcst case2_newfcst case3_newfcst case4_newfcst'
#cases='case1_hum default_hum'

#event=20220608
#events='20200413 20200812 20200903'
#event='20210717'
events='20220608'
#events='20220716'
#events='20220608'


#events='20200812 20200903 20220608 20220716 20210717'
#events='20200812 20200903'
#events='20220608 20220716'
#events='
#event='20210717'

for event in $events;do

for case in $cases;do
echo $event.$case
################
#reflectivity fss
###############
#cd /glade/campaign/univ/umcp0011/colab_scripts
cd /glade/work/jmccurry/ANALYSIS/batch_scripts
#obstype='cropped_osse'
obstype='REFL_10CM'

python process_fss_refl_ensemble_kent.py $mems $event `expr $PBS_ARRAY_INDEX - 1` $scratch $obstype 25 $case
python process_fss_refl_ensemble_kent.py $mems $event `expr $PBS_ARRAY_INDEX - 1` $scratch $obstype 35 $case
#python process_fss_refl_ensemble_kent.py $mems $event 9 $scratch $obstype 25 $case
#python process_fss_refl_ensemble_kent.py $mems $event 9 $scratch $obstype 35 $case

  

#################
#precip fss
#################
obstype='PREC_ACC_NC'
#python process_fss_2d_ensemble_kent.py $mems $event 9 $scratch $obstype 2.5  $case
#python process_fss_2d_ensemble_kent.py $mems $event 9 $scratch $obstype 1.25 $case
#python process_fss_2d_ensemble_kent.py $mems $event 9 $scratch $obstype 0.625 $case
#python process_fss_2d_ensemble_kent.py $mems $event `expr $PBS_ARRAY_INDEX - 1` $scratch $obstype 2.5  $case
#python process_fss_2d_ensemble_kent.py $mems $event `expr $PBS_ARRAY_INDEX - 1` $scratch $obstype 1.25 $case
#python process_fss_2d_ensemble_kent.py $mems $event `expr $PBS_ARRAY_INDEX - 1` $scratch $obstype 0.625 $case

done
done

