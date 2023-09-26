#!/bin/bash
#note* need to modify to handle different obs type numbering between DT2 and Cheyenne implementations of filter
#DOPPLER_RADIAL_VELOCITY is 164 on DT2, 36 on Cheyenne
#multimem #output mems currently fixed at 20. Need to modify to dynamically find #output mems
begin=`awk '$1=="first:"{begin_row=NR} END{print begin_row}' ${1}`
vr_num=`awk '$0 ~ "DOPPLER_RADIAL_VELOCITY" {begin_row=$1} END{print begin_row}' ${1}`
fields=`awk '$1=="num_obs:"{begin_row=NR+1} $1=="first:"{end_row=NR} END{print end_row-begin_row}' ${1}`
contains_truth=`grep "truth" ${1} | wc -l`

echo $fields
echo $3
echo $contains_truth
if [ $contains_truth -eq 1 ];then
if [[ "$3" == "prior" ]];then

awk -v begin="$begin" -v vr_num="$vr_num" -v fields="$fields" 'BEGIN{OFS=","}  NR==0{last=NA;increment=0} increment!=0{increment+=1}  last=="OBS"{value=$1;increment=1} increment==7{val_qc=$1} increment==2{prior=$1} increment==4{post=$1} increment==5{prior_sp=$1} increment==6{post_sp=$1} increment==(fields+6){type=$1} increment==(fields+4){x_loc=$1; y_loc=$2;z_loc=$3} {if(increment==(fields+8) && type!=vr_num){obs_err=$1;print type,value,x_loc,y_loc,z_loc,val_qc,ob_num,prior,post,prior_sp,post_sp,obs_err}} {if(increment==(fields+15) && type==vr_num){obs_err=$1;print type,value,x_loc,y_loc,z_loc,val_qc,ob_num,prior,post,prior_sp,post_sp,obs_err}} NR>begin{last=$1} last=="OBS"{ob_num=$2}' ${1} > ${2}

elif [[ "$3" == "truth" ]];then

awk -v begin="$begin" -v vr_num="$vr_num" -v fields="$fields" 'BEGIN{OFS=","}  NR==0{last=NA;increment=0} increment!=0{increment+=1}  last=="OBS"{value=$1;increment=1} increment==7{val_qc=$1} increment==3{prior=$1} increment==4{post=$1} increment==5{prior_sp=$1} increment==6{post_sp=$1} increment==(fields+6){type=$1} increment==(fields+4){x_loc=$1; y_loc=$2;z_loc=$3} {if(increment==(fields+8) && type!=vr_num){obs_err=$1;print type,value,x_loc,y_loc,z_loc,val_qc,ob_num,prior,post,prior_sp,post_sp,obs_err}} {if(increment==(fields+15) && type==vr_num){obs_err=$1;print type,value,x_loc,y_loc,z_loc,val_qc,ob_num,prior,post,prior_sp,post_sp,obs_err}} NR>begin{last=$1} last=="OBS"{ob_num=$2}' ${1} > ${2}
fi
else

awk -v begin="$begin" -v vr_num="$vr_num" -v fields="$fields" 'BEGIN{OFS=","}  NR==0{last=NA;increment=0} increment!=0{increment+=1}  last=="OBS"{value=$1;increment=1} increment==7{val_qc=$1} increment==2{prior=$1} increment==4{post=$1} increment==5{prior_sp=$1} increment==6{post_sp=$1} increment==(fields+6){type=$1} increment==(fields+4){x_loc=$1; y_loc=$2;z_loc=$3} {if(increment==(fields+8) && type!=vr_num){obs_err=$1;print type,value,x_loc,y_loc,z_loc,val_qc,ob_num,prior,post,prior_sp,post_sp,obs_err}} {if(increment==(fields+15) && type==vr_num){obs_err=$1;print type,value,x_loc,y_loc,z_loc,val_qc,ob_num,prior,post,prior_sp,post_sp,obs_err}} NR>begin{last=$1} last=="OBS"{ob_num=$2}' ${1} > ${2}

fi
# last=="kind"{flag=3}
