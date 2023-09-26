#!/bin/bash
######
#output a summary space delimited file containing type,value, x,y,z locs ,QC val,and ob number
#currently only works for radial velocity obs_seq files
#####
begin=`awk '$1=="first:"{begin_row=NR} END{print begin_row}' ${1}`
vr_num=`awk '$0 ~ "DOPPLER_RADIAL_VELOCITY" {begin_row=$1} END{print begin_row}' ${1}`
fields=`awk '$1=="num_obs:"{begin_row=NR+1} $1=="first:"{end_row=NR} END{print end_row-begin_row}' ${1}`

echo $fields
echo $begin
# awk -v begin="$begin" 'BEGIN{OFS=",";temp=-888888;val_qc=-888888}  NR==0{flag=0;last=NA;increment=0} increment!=0{increment+=1} {flag=0} last=="OBS"{flag=1;increment=1} last=="kind"{flag=3} {if(last=="loc3d" && increment==6){flag=2}} increment==16{time=$1} flag==1{value=$1} flag==2{x_loc=$1; y_loc=$2;z_loc=$3} flag==3{type=$1} {if(flag==3 && type!="164"){print type,value,x_loc,y_loc,z_loc,val_qc,ob_num}} {if(increment==13 && type=="164"){nan=$1}} {if(increment==17 && type=="164"){print type,value,x_loc,y_loc,z_loc,val_qc,ob_num}} NR>begin{last=$1} last=="OBS"{ob_num=$2}' ${1} > ${2}

awk -v begin="$begin" -v vr_num="$vr_num" -v fields="$fields" 'BEGIN{OFS=",";temp=-888888;val_qc=-888888}  NR==0{flag=0;last=NA;increment=0} increment!=0{increment+=1} {flag=0} last=="OBS"{increment=1;value=$1} last=="kind"{type=$1}  {if(last=="loc3d" && increment==(4+fields)){x_loc=$1; y_loc=$2;z_loc=$3}} {if (increment==(8+fields) && type!=vr_num){oerror=$1}} {if (increment==(15+fields) && type==vr_num){oerror=$1}} {if (increment==(7+fields) && type!=vr_num){day=$1; second=$2}} {if (increment==(14+fields) && type==vr_num){day=$1; second=$2}} {if(increment==(8+fields) && type!=vr_num){print type,value,x_loc,y_loc,z_loc,val_qc,ob_num,temp,temp,temp,day,second,oerror,temp,temp,temp}} {if(increment==(9+fields) && type==vr_num){radar_xloc=$1;radar_yloc=$2;radar_zloc=$3}} {if(increment==(11+fields) && type==vr_num){radar_xdir=$1;radar_ydir=$2;radar_zdir=$3}} {if(increment==(15+fields) && type==vr_num){print type,value,x_loc,y_loc,z_loc,val_qc,ob_num,radar_xloc,radar_yloc,radar_zloc,day,second,oerror,radar_xdir,radar_ydir,radar_zdir}} NR>begin{last=$1} last=="OBS"{ob_num=$2}' ${1} > ${2}

#awk -v begin="$begin" 'BEGIN{OFS=",";temp=-888888;val_qc=-888888}  NR==0{flag=0;last=NA;increment=0} increment!=0{increment+=1} {flag=0} last=="OBS"{increment=1;value=$1} last=="kind"{type=$1}  {if(last=="loc3d" && increment==6){x_loc=$1; y_loc=$2;z_loc=$3}} {if (increment==10 && type!="36"){oerror=$1}} {if (increment==17 && type=="36"){oerror=$1}} {if (increment==9 && type!="36"){day=$1; second=$2}} {if (increment==16 && type=="36"){day=$1; second=$2}} {if(increment==10 && type!="36"){print type,value,x_loc,y_loc,z_loc,val_qc,ob_num,temp,temp,temp,day,second,oerror,temp,temp,temp}} {if(increment==11 && type=="36"){radar_xloc=$1;radar_yloc=$2;radar_zloc=$3}} {if(increment==13 && type=="36"){radar_xdir=$1;radar_ydir=$2;radar_zdir=$3}} {if(increment==17 && type=="36"){print type,value,x_loc,y_loc,z_loc,val_qc,ob_num,radar_xloc,radar_yloc,radar_zloc,day,second,oerror,radar_zdir,radar_ydir,radar_zdir}} NR>begin{last=$1} last=="OBS"{ob_num=$2}' ${1} > ${2}
