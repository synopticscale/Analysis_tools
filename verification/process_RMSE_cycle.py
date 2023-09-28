from datetime import datetime
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import metpy.calc as mpcalc
from metpy.interpolate import log_interpolate_1d
from metpy.units import units
import numpy as np
from scipy.ndimage import gaussian_filter
from scipy.interpolate import griddata
import xarray as xr
from numpy import dtype
import os
import netCDF4 as nc
from mrms_util import readwdssii, haversine
import pandas as pd
from scipy.spatial.distance import cdist
from netCDF4 import Dataset as netcdf_dataset
import wrf
from wrf import to_np
import wrfpy_OSD 
from datetime import datetime, timedelta
import sys
from scipy import stats
import simple_utils
from scipy.interpolate import RegularGridInterpolator


###########################
#global variables and function definitions
scratch=os.getenv('scratch')
campaign=os.getenv('campaign_nature_store')
verification_boxes=['6km_maryland']

exp_dict={20200413:'EXP1',20200812:'EXP9',20200903:'EXP3',20210717:'EXP4',20220608:'EXP5',20220702:'EXP6',20220716:'EXP7'}
config_dict={'default':'config_0','case1':'config_1','case2':'config_2','case3':'config_3','case4':'config_4'}
def get_root(infile_d01):
  rootgroup = netcdf_dataset(infile_d01, 'r')
  return rootgroup
def regrid(data, out_x, out_y):
    m = max(data.shape[0], data.shape[1])
    y = np.linspace(0, 1.0/m, data.shape[0])
    x = np.linspace(0, 1.0/m, data.shape[1])
    interpolating_function = RegularGridInterpolator((y, x), data)

    yv, xv = np.meshgrid(np.linspace(0, 1.0/m, out_y), np.linspace(0, 1.0/m, out_x))

    return interpolating_function((xv, yv))
###########################
#main program

if __name__ == "__main__":
   
    #######################
    #KENT
    exp = exp_dict[int(sys.argv[1])]
    case = str(sys.argv[2])
    config = config_dict[case]
    
    #######################
    #JMC
    #exp = str(sys.argv[1])
    #config = str(sys.argv[2])
    #######################
    templist = os.listdir('{}/{}/CYCLE/{}'.format(scratch,exp,config)) #KENT

    reallist = [x[-12:] for x in templist if "202" in x]
    reallist.sort()
    reallist=reallist[4:-1]
    variables=str(sys.argv[3]).split(",")
    #begin loop
    case_data=np.empty((len(reallist),3,len(variables)))
    for t,cycle_date in enumerate(reallist):
            inputDir = '{}/{}/CYCLE/{}/{}/filter_out'.format(scratch,exp,config,cycle_date) #CYCLE
            natureDir = '{}/{}.naturerun/WRFOUT/{}'.format(campaign,str(sys.argv[1])[0:8],cycle_date)

            for n,verif in enumerate(verification_boxes):
                nature_root=get_root('{}/nature_{}'.format(natureDir,verif))
                case_root=get_root('{}/preassim_mean.nc'.format(inputDir))

                for n_var,variable in enumerate(variables):
                 print('working on variable: {} init: {}'.format(variable,cycle_date))

                 if(variable[-4:]=='_500'):
                    #load arrays for naturerun and cases
                    nature_np=np.asarray(nature_root[variable[0]])[0,18,:,:].squeeze()
                    nature_np[nature_np<-1000]=np.nan
                    
                    case_np=np.asarray(case_root[variable[0]])[0,18,:,:].squeeze()
                    case_np=regrid(case_np,149,149)
                    case_np[case_np<-1000]=np.nan

                    case_data[t,n+2,n_var]=np.sqrt(np.nanmean((nature_np-case_np)**2))
                    case_data[t,0,n_var]= cycle_date
                    case_data[t,1,n_var]= t*15
                 elif(variable[-4:]=='_850'):
                    #load arrays for naturerun and cases
                    nature_np=np.asarray(nature_root[variable[0]])[0,9,:,:].squeeze()
                    nature_np[nature_np<-1000]=np.nan
                    
                    case_np=np.asarray(case_root[variable[0]])[0,9,:,:].squeeze()
                    case_np=regrid(case_np,149,149)
                    case_np[case_np<-1000]=np.nan

                    case_data[t,n+2,n_var]=np.sqrt(np.nanmean((nature_np-case_np)**2))
                    case_data[t,0,n_var]= cycle_date
                    case_data[t,1,n_var]= t*15               
                 else:
                    #load arrays for naturerun and cases
                    nature_np=np.asarray(nature_root[variable]).squeeze()
                    nature_np[nature_np<-1000]=np.nan
                    
                    case_np=np.asarray(case_root[variable]).squeeze()
                    case_np=regrid(case_np,149,149)
                    case_np[case_np<-1000]=np.nan

                    case_data[t,n+2,n_var]=np.sqrt(np.nanmean((nature_np-case_np)**2))
                    case_data[t,0,n_var]= cycle_date
                    case_data[t,1,n_var]= t*15
                    
           

        #load into pandas dataframes
    for n_var,variable in enumerate(variables):
          case_df = pd.DataFrame(case_data[:,:,n_var].squeeze(),columns=['init','time','6km_maryland'])
          os.makedirs('/glade/work/jmccurry/WOF/VERIF/RMSE/{}.{}'.format(exp,config),exist_ok=True)
          case_df.to_pickle('/glade/work/jmccurry/WOF/VERIF/RMSE/{}.{}/df_CYCLE_{}'.format(exp,config,variable))
    exit()
