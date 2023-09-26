#!/bin/bash
######
#output a summary space delimited file containing type,value, x,y,z locs ,QC val,and ob number
#currently only works for radial velocity obs_seq files
#####
begin=`awk '$1=="first:"{begin_row=NR} END{print begin_row}' ${1}`
vr_num=`awk '$0 ~ "DOPPLER_RADIAL_VELOCITY" {begin_row=$1} NR==40{print begin_row}' ${1}`
echo $begin
echo $vr_num
awk -v begin="$begin" -v vr_num="$vr_num" 'BEGIN{OFS=",";temp=-888888;val_qc=-888888}  NR==0{flag=0;last=NA;increment=0} increment!=0{increment+=1} {flag=0} last=="OBS"{increment=1;value2=$1} {if(increment==2){value=$1}} last=="kind"{type=$1}  {if(last=="loc3d" && increment==7){x_loc=$1; y_loc=$2;z_loc=$3}} {if (increment==11 && type!=vr_num){oerror=$1}} {if (increment==18 && type==vr_num){oerror=$1}} {if (increment==10 && type!=vr_num){day=$1; second=$2}} {if (increment==17 && type==vr_num){day=$1; second=$2}} {if(increment==11 && type!=vr_num){print type,value,x_loc,y_loc,z_loc,val_qc,ob_num,temp,temp,temp,day,second,oerror,temp,temp,temp,value2}} {if(increment==12 && type==vr_num){radar_xloc=$1;radar_yloc=$2;radar_zloc=$3}} {if(increment==14 && type==vr_num){radar_xdir=$1;radar_ydir=$2;radar_zdir=$3}}  {if(increment==18 && type==vr_num){print type,value,x_loc,y_loc,z_loc,val_qc,ob_num,radar_xloc,radar_yloc,radar_zloc,day,second,oerror,radar_xdir,radar_ydir,radar_zdir,value2}} NR>begin{last=$1} last=="OBS"{ob_num=$2}' ${1} > ${2}
