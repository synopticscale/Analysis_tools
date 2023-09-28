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
##########################
#global variables
##########################
#scratch='/glade/scratch/jmccurry/WOF/realtime'
scratch=str(sys.argv[4])
exp_dict={20200413:'EXP21',20200812:'EXP17',20200903:'EXP23',20210717:'EXP24',20220608:'EXP25',20220702:'EXP26',20220716:'EXP27'}
config_dict={'default':'config_0','case1':'config_1','case2':'config_2','case3':'config_3','case4':'config_4'}
obsDir='/glade/scratch/jmccurry/WOF/GRIDDED_OBS/{}'.format(sys.argv[2])
#obs_type='cropped_osse'
obs_type=str(sys.argv[5])
prescribed_thresh=int(sys.argv[6])
FSS_MASK='/glade/work/jmccurry/OSSE/MARYLAND_NP_MASK_NEW.npy'
template='/glade/scratch/jmccurry/WOF/ensemble_init/20200812/output/2020081211/wrfinput_d01_153260_39600_mean'
##########################
#grab necessary variables
##########################
def get_vars(infile_d01):
  rootgroup = netcdf_dataset(infile_d01, 'r')
  QRAIN = wrf.getvar(rootgroup,'QRAIN')[:,:,:]
  QGRAUP = wrf.getvar(rootgroup,'QGRAUP')[:,:,:]
  QVAPOR = wrf.getvar(rootgroup,'QVAPOR')[:,:,:]
  T_SURF = wrf.getvar(rootgroup,'T')[0,:,:]
  return QRAIN,QGRAUP,QVAPOR,T_SURF
  

def mask_obs(obspath,wrfinput_d01=template):
    airlist=[]
    with open('NEXRAD_loc5.txt',mode='r') as f:
      airlist.append([c.rstrip().split()[1:] for c in f])
      [airlist] = airlist
    airlist = np.array(airlist).astype(float)
    rootgroup = netcdf_dataset(wrfinput_d01, 'r')
    var = wrf.getvar(rootgroup, 'T')
    xy = to_np(wrf.ll_to_xy(rootgroup,airlist[:,0],airlist[:,1]))
    lats = xy[0][(xy[0]>=-50) & (xy[0]<349) & (xy[1]>=-50) & (xy[1]<349)]
    lons = xy[1][(xy[0]>=-50) & (xy[0]<349) & (xy[1]>=-50) & (xy[1]<349)]
    airp = np.array([lons,lats]).T

    #obs_test = np.load('/glade/campaign/univ/umcp0011/GRIDDED_OBS/type_165_{}.npy'.format(timestr))
    obs_test = np.load(obspath) #NEED TO FIX TO ACCOUNT FOR OBS OTHER THAN COMPOSITE REFL

    obs_filtered = np.zeros((299,299))
    for i in range(0,np.shape(obs_test)[0]):
        for j in range(0,np.shape(obs_test)[1]):
            dist_array = cdist(np.array([[i,j]]), airp)*3
            dist_array = dist_array[dist_array>10]
            dist = dist_array.min()

            if(dist > 150):
                obs_filtered[i,j] = np.nan
            elif((i<7) or (i > 291) or (j < 7) or (j > 291)):
                obs_filtered[i,j] = np.nan
            else:
                obs_filtered[i,j] = obs_test[i,j]
    return obs_filtered

def ensemble_probabilities(exp,config,forecast_init,timestamp,thresh,mems=20):
    #global scratch
    print('Processing probabilities for date '+str(timestamp))
    
    obspath = obsDir+'/{}.{}.npy'.format(obs_type,timestamp)
    #ncin_o = mask_obs(obspath) #obs data
    ncin_o = np.load(obspath)

    inputDir = '{}/{}/CYCLE/{}/forecasts/WRFOUTS_FCST{}/{}'.format(scratch,exp,config,forecast_init,timestamp)
    thres = int(thresh)
    for i in range(mems):
        fpath = inputDir+'/wrfout_d01_forecast_'+ timestamp + '_{}'.format(i+1)
        ncin_f = nc.Dataset(fpath, 'r', format='NETCDF4') #fcst data
        # Extract variables
        reflh_preprocess = ncin_f.variables['REFL_10CM'][:,:,:]
        reflh = np.amax(reflh_preprocess,axis=1).squeeze()
        reflh[np.isnan(ncin_o)] = np.nan       
        # create field of ensemble probabilities
        if (i==0):
            model_probabilities = np.zeros_like(reflh)
            obs_probablities = np.zeros_like(reflh)
        model_probabilities += np.where(reflh>=thres, 1, 0)
        
        try:
             model_probabilities[np.isnan(reflh)] = np.nan
        except:
             pass
    obs_probabilities_f = np.where(ncin_o>=thres, 1, 0)
    try:
        obs_probabilities_f[np.isnan(ncin_o)] = np.nan
    except:
        pass
    model_probabilities_f = model_probabilities/mems
    print(np.nansum(model_probabilities))
    return model_probabilities_f, obs_probabilities_f
    
        
    # create field of obs probabilities
def evaluate_fss(model_probabilities, obs_probabilities):
    roi = [1,5,45]  # About 3km, 6 km, 15 km, 45 km and 135km respectively
    FN = [] #false negatives
    FP = [] #false positives
    TN = [] #true negatives
    TP = [] #true positives
    O = [] #observed fractions
    M = [] #modeled fractions
    possible = np.zeros((len(roi)+1,np.shape(model_probabilities)[0],np.shape(model_probabilities)[1])) #field showing total non-NAN gridpoints for scaling 
    if FSS_MASK!='none':
        MASK=np.load(FSS_MASK)
    for i in range(35,264): #change to allow variability in domain size
        for j in range(35,264):   
                if FSS_MASK!='none':
                    if(MASK[j,i]<0.5):
                       continue
                # Loop through roi's
                obcount = []
                modcount = []
                fneg = []
                fpos = []
                tneg = []
                tpos = []
                
                ri = 0
                for r in roi:    
                    # Take subset of points arround current one
                    a_in = obs_probabilities[max(j-r,0):min(j+r+1,299),max(i-r,0):min(i+r+1,299)]
                    b_in = model_probabilities[max(j-r,0):min(j+r+1,299),max(i-r,0):min(i+r+1,299)]
                    # calculate number of valid elements 
                    
                    possible[ri,j,i] = float(np.count_nonzero(~np.isnan(a_in)))
                    #if > 50% of FSS window is within boundaries and not outside radar FOI or otherwise NAN, record fractions
                    if(possible[ri,j,i]>=0.5*(((2*r)+1)**2)):
                        #print('{} valid neighboring points at i:{} j:{}'.format(possible,i,j))
                        # Number of obs field points where threshold is exceeded
                 
                        obcount.append(float(np.nansum(a_in))/possible[ri,j,i])
                        
                        # Number of model points where threshold is exceeded

                        modcount.append(float(np.nansum(b_in))/possible[ri,j,i])

                        # Store obs fraction calculations for each roi 
                        #refl3[ri,j,i] = float(((np.nansum(a_in)/possible[ri,j,i])-(np.nansum(b_in)/possible[ri,j,i])))
                        
                        fneg.append(np.nansum(np.where((a_in==1)&(b_in<1),1-b_in,0))) 
                        fpos.append(np.nansum(np.where((a_in==0)&(b_in>0),b_in,0)))
                        tneg.append(np.nansum(np.where((a_in==0)&(b_in<1),1-b_in,0)))
                        tpos.append(np.nansum(np.where((a_in==1)&(b_in>0),b_in,0)))
                    # if < 50% of FSS window is valid, then record NAN's instead of fractions 
                    else:
                   
                        
                        obcount = [np.nan for r in roi]

                        # Number of model points where threshold is exceeded

                        modcount = [np.nan for r in roi]

                        # Store obs fraction calculations for each roi 
                        #refl3[:,j,i] = np.nan
                          
                        #detection theory stuff
                        fneg = [np.nan for r in roi]
                        fpos = [np.nan for r in roi]
                        tneg = [np.nan for r in roi]
                        tpos = [np.nan for r in roi]
                        break 
                    ri = ri + 1                       

                O.append(obcount)
                M.append(modcount)
                FN.append(fneg)
                FP.append(fpos)
                TN.append(tneg)
                TP.append(tpos)
                    # Loop through roi's
   

    FN = np.array(FN)
    FP = np.array(FP)
    TN = np.array(TN)
    TP = np.array(TP)
    O = np.array(O)
    M = np.array(M)
    fss = []
    bias = []
    pod =[]
    false_alarm = []
    for i in range(np.size(roi)):
        pod.append(1 - np.nanmean(FN[:,i])/np.nanmean((FN[:,i]+TP[:,i])))
        false_alarm.append(np.nanmean(FP[:,i])/np.nanmean((FP[:,i]+TN[:,i])))
        MSE = np.nansum( (O[:,i] - M[:,i])**2 ) / np.count_nonzero(~np.isnan(O[:,i]))
        MSE_ref = ( np.nansum( O[:,i]**2 ) + np.nansum( M[:,i]**2 ) ) / np.count_nonzero(~np.isnan(O[:,i]))
        fss.append(1.0 - (MSE / MSE_ref))
        bias.append(np.nansum( (O[:,i] - M[:,i]) ) / np.count_nonzero(~np.isnan(O[:,i])))
    return fss,bias,false_alarm,pod                 
if __name__ == "__main__":
    mems = int(sys.argv[1])
    exp = exp_dict[int(sys.argv[2])]
    case = str(sys.argv[7])
    config = config_dict[case]
    templist = os.listdir('{}/{}/CYCLE/{}/forecasts'.format(scratch,exp,config))
    print('{}/{}/CYCLE/{}/forecasts'.format(scratch,exp,config))
    reallist = [x[-12:] for x in templist if "WRFOUTS_FCST20" in x]
    reallist.sort()
    forecast_init = reallist[int(sys.argv[3])]

    fname_final='/glade/work/jmccurry/WOF/VERIF/FSS/{}.{}/df_{}_thresh{}_MARYLAND'.format(exp,config,prescribed_thresh,forecast_init)
    start_date = datetime.strptime(forecast_init,"%Y%m%d%H%M")
    end_date = start_date + timedelta(minutes=180)
    forecast_dates = simple_utils.cycle_timer(start_date,end_date,15)  
    for n,forecast_date in enumerate(forecast_dates):
      print('working on time: {} init: {}'.format(forecast_date,forecast_init))
      indie = n*15
      model_probs,obs_probs = ensemble_probabilities(exp,config,forecast_init,forecast_date,prescribed_thresh,mems=mems)
      fss,bias,falarm,pod = evaluate_fss(model_probs,obs_probs)
      if(n==0):
             df = pd.DataFrame(np.array([[start_date,indie,fss[0],fss[1],fss[2],bias[0],bias[1],bias[2],pod[0],pod[1],pod[2],falarm[0],falarm[1],falarm[2]]]),columns=['init','time','3km_FSS','15km_FSS','135km_FSS','3km_BIAS','15km_BIAS','135km_BIAS','3km_POD','15km_POD','135km_POD','3km_FALARM','15km_FALARM','135km_FALARM']).set_index('time')        
      else:
             df = df.append(pd.DataFrame(np.array([[start_date,indie,fss[0],fss[1],fss[2],bias[0],bias[1],bias[2],pod[0],pod[1],pod[2],falarm[0],falarm[1],falarm[2]]]),columns=['init','time','3km_FSS','15km_FSS','135km_FSS','3km_BIAS','15km_BIAS','135km_BIAS','3km_POD','15km_POD','135km_POD','3km_FALARM','15km_FALARM','135km_FALARM']).set_index('time'))        

    df= df.reset_index()
    os.makedirs('/glade/work/jmccurry/WOF/VERIF/FSS/{}.{}'.format(exp,config),exist_ok=True)
    df.to_pickle('/glade/work/jmccurry/WOF/VERIF/FSS/{}.{}/df_{}_thresh{}_MARYLAND'.format(exp,config,forecast_init,prescribed_thresh))
                 
                 
