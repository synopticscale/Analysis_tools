from datetime import datetime, timedelta
def cycle_timer(startdate,enddate,inc):
    start_date = startdate
    end_date = enddate
    delta = timedelta(minutes = inc)
    datestrings = []
    #this is a stopgap solution. Need to find number of obs and increment ob_buffer by that
    #will break if > 1million obs per timestamp
    ob_buffer = 1000000 
    while start_date <= end_date:
      datestrings.append(start_date.strftime("%Y%m%d%H%M"))
      start_date += delta
    return datestrings
def forecast_timer(startdate,enddate,init_inc,forecast_inc,forecast_len):
    delta = timedelta(minutes = init_inc)
    delta_forecast = timedelta(minutes = forecast_inc)
    forecast_len = timedelta(minutes = forecast_len)
    datestrings = []
    datestrings_all = []
    while startdate <= enddate:
         datestrings.append(startdate.strftime("%Y%m%d%H%M"))
         startdate += delta
    for timestr_init in enumerate(datestrings):
        datetime_init = datetime.strptime(timestr_init[1],'%Y%m%d%H%M')
        start_date = datetime.strptime(timestr_init[1],'%Y%m%d%H%M')
        end_date = datetime.strptime(timestr_init[1],'%Y%m%d%H%M')+forecast_len
        datestrings_fcst = []
        while start_date <= end_date:
            datestrings_fcst.append(start_date.strftime("%Y%m%d%H%M"))
            start_date += delta_forecast
        datestrings_all.append(datestrings_fcst)
    return datestrings_all
def forecast_init_timer(startdate,enddate,init_inc):
    delta = timedelta(minutes = init_inc)
    datestrings = []
    while startdate <= enddate:
         datestrings.append(startdate.strftime("%Y%m%d%H%M"))
         startdate += delta
    return datestrings