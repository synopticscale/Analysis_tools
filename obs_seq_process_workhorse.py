import matplotlib.pyplot as plt
from pylab import *
import numpy as np 
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import FormatStrFormatter
import sys
from os import path 
from os import system
import os
import subprocess
from scipy.interpolate import interp1d
import math
from netCDF4 import Dataset as netcdf_dataset


#wrf-python libs
import wrf
from wrf import to_np, vertcross, CoordPair
from matplotlib.cm import get_cmap

#########
#written by JMC
#plot obs from obs_sequence file 
#for diagnostic purposes
########

vr_num=36 #change vr_num here

class obs_seq:
   def __init__(self,timestamp,obs_seq,custom='no',perfect_obs='no'):
        self.data = load_obs_seq(timestamp,obs_seq,custom,perfect_obs)
        self.varnames = dict([('obs_type',0),('value',1),('vals',1),('X_loc',2),('Y_loc',3),('Z_loc',4),('QC',5),('obs_num',6),('radar_xloc',7),('radar_yloc',8),('second',10),('day',11),('oerror',12),('value2',13)])
        self.radar_locs = dict([('KVNX',[4.57053374392750,0.6412447157008800]),('KOAX',[4.601266986004480,0.7211764997633801]),('KEAX',[4.637959886345540,0.6773666343042200]),('KDVN',[4.702251524537250,0.7262606074423900])])
   def filter_data_type(self,obs_type):
      data_int = self.data.astype(int)
      self.data = self.data[(data_int[:,0]==obs_type)]
      return self 
   def filter_data_type_list(self,obs_type_list):
      data_int = self.data.astype(int)
      self.data = self.data[np.where(np.isin(data_int[:,0],obs_type_list))]
      return self 
   def return_data_type(self,obs_type):
      data_int = self.data.astype(int)
      new = self.data[(data_int[:,0]==obs_type)]
      return new
   def return_data_type_list(self,obs_type_list):
      data_int = self.data.astype(int)
      new = self.data[np.where(np.isin(data_int[:,0],obs_type_list))]
      return new 
   def filter_outliers(self,z,low_thresh,high_thresh):
      self.data[:,self.varnames[str(z)]][self.data[:,self.varnames[str(z)]]< low_thresh] = 'NaN'  
      self.data[:,self.varnames[str(z)]][self.data[:,self.varnames[str(z)]]> high_thresh] = 'NaN' 
      self.data = self.data[~np.isnan(self.data[:,self.varnames[str(z)]])]
      return self 
   def filter_radar_location(self,radar):
      self.data[:,self.varnames['radar_xloc']][self.data[:,self.varnames['radar_xloc']]!=self.radar_locs[radar][0]] = 'NaN'  
      self.data[:,self.varnames['radar_yloc']][self.data[:,self.varnames['radar_yloc']]!=self.radar_locs[radar][1]] = 'NaN'  
      self.data = self.data[~np.isnan(self.data[:,self.varnames['radar_xloc']])]
      self.data = self.data[~np.isnan(self.data[:,self.varnames['radar_yloc']])]
      return self 

class obs_seq_final:
   def __init__(self,timestamp,obs_seq,rundir='',forecast='no',outputmem='',member='no',custom='no',output='prior'):

    if ( forecast=='no'):
     self.data = load_obs_final(rundir,timestamp,obs_seq,outputmem,custom,output)

    self.varnames = dict([('obs_type',0),('value',1),('vals',1),('X_loc',2),('Y_loc',3),('Z_loc',4),('QC',5),('obs_num',6),('prior',7),('truth',7),('post',8),('prior_sp',9),('post_sp',10),('obs_err',11),('priorA',7),('priorB',8),('priorC',9),('priorD',10),('priorE',11),('priorF',12),('priorG',13),('priorH',14),('priorI',15),('priorJ',16),('priorK',17),('priorL',18),('priorM',19),('priorN',20),('priorO',21),('priorP',22),('priorQ',23),('priorR',24),('priorS',25),('priorT',26)])
   def filter_data_QC(self,QC):
      data_int = self.data.astype(int)
      self.data = self.data[(data_int[:,5]==QC)]
      return self
   def filter_data_QC_list(self,QC):
      data_int = self.data.astype(int)
      self.data = self.data[np.where(np.isin(data_int[:,5],QC))]
      return self 
   def filter_data_type(self,obs_type):
      data_int = self.data.astype(int)
      self.data = self.data[(data_int[:,0]==obs_type)]
      return self 
   def filter_outliers(self,z,low_thresh,high_thresh):
      self.data[:,self.varnames[str(z)]][self.data[:,self.varnames[str(z)]]< low_thresh] = 'NaN'  
      self.data[:,self.varnames[str(z)]][self.data[:,self.varnames[str(z)]]> high_thresh] = 'NaN'  
      self.data = self.data[~np.isnan(self.data[:,self.varnames[str(z)]])]
      return self 
   def filter_outliers_range(self,z1,z2,low_thresh,high_thresh):
      for z in range(z1,z2):
          self.data[:,z][self.data[:,z]< low_thresh] = 'NaN'  
          self.data[:,z][self.data[:,z]> high_thresh] = 'NaN'  
          self.data = self.data[~np.isnan(self.data[:,z])]
      return self 
   def filter_obs_num_list(self,nums):
      data_int = self.data.astype(int)
      self.data = self.data[np.where(np.isin(data_int[:,6],nums))]
      return self 
   def return_data_QC(self,QC):
      data_int = self.data.astype(int)
      new = self.data[(data_int[:,5]==QC)]
      return new
   def return_data_type(self,obs_type):
      data_int = self.data.astype(int)
      new = self.data[(data_int[:,0]==obs_type)]
      return new 

def load_obs_final(rundir,timestamp,obs_seq,outputmem,custom='no',output='prior'):
#eventually update so that 'custom' is default setting, remove rundir, timestamp, obs_seq
   if (custom=='no'):
    procname = '{}/obs_seq.processed.{}.{}'.format(os.getcwd(),timestamp,obs_seq)
    path='/glade/scratch/jmccurry/WOF/realtime/{}/{}/obs_seq.final.{}.{}.{}'.format(rundir,timestamp,timestamp,obs_seq,timestamp)
   else:
    procname = '{}'.format(custom[0])
    path ='{}'.format(custom[1])        
   subprocess.call(['{}/process_obs_final.sh'.format(os.getcwd()),'{}'.format(path),'{}'.format(procname),'{}'.format(output)])
   data = np.loadtxt(procname,delimiter=',')
   return data

def load_obs_seq(timestamp,obs_seq,custom='no',perfect_obs='no'):

   if (custom=='no'):
    procname = '{}/obs_seq.processed.{}.{}'.format(os.getcwd(),obs_seq,timestamp)
    path='/glade/scratch/jmccurry/WOF/realtime/OBSGEN/OBS_SEQ_OSSE/obs_seq.{}.{}'.format(obs_seq,timestamp)
   else:
    procname = '{}'.format(custom[0])
    path ='{}'.format(custom[1])
   if (perfect_obs=='no'):
    subprocess.call(['{}/process_obs_seq.sh'.format(os.getcwd()),'{}'.format(path),'{}'.format(procname),str(vr_num)])
    if('platform' in open(procname).read()):
        subprocess.call(['{}/process_obs_seq.sh'.format(os.getcwd()),'{}'.format(path),'{}'.format(procname),str(164)])
   else:
    subprocess.call(['{}/process_obs_seq_truth.sh'.format(os.getcwd()),'{}'.format(path),'{}'.format(procname),str(vr_num)])
    if('platform' in open(procname).read()):

        subprocess.call(['{}/process_obs_seq_truth.sh'.format(os.getcwd()),'{}'.format(path),'{}'.format(procname),str(164)])

   
   data = np.loadtxt(procname,delimiter=',')
   return data
def main():
   pass
if __name__ == "__main__":
    main()
