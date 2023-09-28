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

###########################
#global variables and function definitions
scratch=os.getenv('scratch')
campaign=os.getenv('campaign_nature_store')
verification_boxes=['6km','6km_midatlantic','6km_maryland','6km_pop']
exp_dict={20200413:'EXP21',20200812:'EXP17',20200903:'EXP23',20210717:'EXP24',20220608:'EXP25',20220702:'EXP26',20220716:'EXP27'}
config_dict={'default':'config_0','case1':'config_1','case2':'config_2','case3':'config_3','case4':'config_4'}
def get_root(infile_d01):
  rootgroup = netcdf_dataset(infile_d01, 'r')
  return rootgroup

###########################
#main program

if __name__ == "__main__":
   
    #######################
    #KENT
    exp = exp_dict[int(sys.argv[1])]
    case = str(sys.argv[5])
    config = config_dict[case]
    
    #######################
    #JMC
    #exp = str(sys.argv[1])
    #config = str(sys.argv[5])
    #######################
    templist = os.listdir('{}/{}/CYCLE/{}/forecasts'.format(scratch,exp,config)) #KENT
    #templist = os.listdir('{}/{}.{}'.format(scratch,exp,config)) #JMC

    reallist = [x[-12:] for x in templist if "WRFOUTS_FCST20" in x]
    reallist.sort()
    forecast_init = reallist[int(sys.argv[2])]
    start_date = datetime.strptime(forecast_init,"%Y%m%d%H%M")
    end_date = start_date + timedelta(minutes=180)
    forecast_dates = simple_utils.cycle_timer(start_date,end_date,15)  
    variables=str(sys.argv[4]).split(",")
    checklist=0
    for n_var,variable in enumerate(variables):
      fname='/glade/work/jmccurry/WOF/VERIF/RMSE/{}.{}/df_{}_{}'.format(exp,config,forecast_init,variable)  
 
      if(not os.path.isfile(fname)):
       print('missing /glade/work/jmccurry/WOF/VERIF/RMSE/{}.{}/df_{}_{}'.format(exp,config,forecast_init,variable))
       checklist+=1
      else:
       pass
       #variables.remove(variable)  
    if(checklist<1):
      pass
      #print('nothing to do')
      #exit()

    #begin loop
    case_data=np.empty((13,6,len(variables)))
    for t,forecast_date in enumerate(forecast_dates):
            inputDir = '{}/{}/CYCLE/{}/forecasts/WRFOUTS_FCST{}/{}'.format(scratch,exp,config,forecast_init,forecast_date) #KENT
            #inputDir = '{}/{}.{}/WRFOUTS_FCST{}/{}'.format(scratch,exp,config,forecast_init,forecast_date) #JMC
            natureDir = '{}/{}.naturerun/WRFOUT/{}'.format(campaign,str(sys.argv[1])[0:8],forecast_date)

            for n,verif in enumerate(verification_boxes):
                nature_root=get_root('{}/nature_{}'.format(natureDir,verif))
                try:
                   case_root=get_root('{}/forecast_{}'.format(inputDir,verif))
                except:
                   print('file {} is corrupt -- removing'.format('{}/forecast_{}'.format(inputDir,verif)))
                   if os.path.exists('{}/forecast_{}'.format(inputDir,verif)):
                       os.remove('{}/forecast_{}'.format(inputDir,verif))

                for n_var,variable in enumerate(variables):
                 print('working on variable: {} init: {}'.format(variable,forecast_init))

                 if(variable[-4:]=='_500'):
                    #load arrays for naturerun and cases
                    nature_np=np.asarray(nature_root[variable[0]])[0,18,:,:].squeeze()
                    nature_np[nature_np<-1000]=np.nan
                    
                    case_np=np.asarray(case_root[variable[0]])[0,18,:,:].squeeze()
                    case_np[case_np<-1000]=np.nan

                    case_data[t,n+2,n_var]=np.sqrt(np.nanmean((nature_np-case_np)**2))
                    case_data[t,0,n_var]= forecast_init
                    case_data[t,1,n_var]= t*15
                 elif(variable[-4:]=='_850'):
                    #load arrays for naturerun and cases
                    nature_np=np.asarray(nature_root[variable[0]])[0,9,:,:].squeeze()
                    nature_np[nature_np<-1000]=np.nan
                    
                    case_np=np.asarray(case_root[variable[0]])[0,9,:,:].squeeze()
                    case_np[case_np<-1000]=np.nan

                    case_data[t,n+2,n_var]=np.sqrt(np.nanmean((nature_np-case_np)**2))
                    case_data[t,0,n_var]= forecast_init
                    case_data[t,1,n_var]= t*15               
                 else:
                    #load arrays for naturerun and cases
                    nature_np=np.asarray(nature_root[variable]).squeeze()
                    nature_np[nature_np<-1000]=np.nan
                    
                    case_np=np.asarray(case_root[variable]).squeeze()
                    case_np[case_np<-1000]=np.nan

                    case_data[t,n+2,n_var]=np.sqrt(np.nanmean((nature_np-case_np)**2))
                    case_data[t,0,n_var]= forecast_init
                    case_data[t,1,n_var]= t*15
                    
           

        #load into pandas dataframes
    for n_var,variable in enumerate(variables):
          case_df = pd.DataFrame(case_data[:,:,n_var].squeeze(),columns=['init','time','6km','6km_midatlantic','6km_maryland','6km_pop'])
          os.makedirs('/glade/work/jmccurry/WOF/VERIF/RMSE/{}.{}'.format(exp,config),exist_ok=True)
          case_df.to_pickle('/glade/work/jmccurry/WOF/VERIF/RMSE/{}.{}/df_{}_{}'.format(exp,config,forecast_init,variable))
    exit()
