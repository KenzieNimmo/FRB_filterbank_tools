###requires astropy 2.0 and python3. Will not work with numpy v1.14 due to a numpy function bug###
###Kenzie Nimmo 2020

import pandas as pd
import numpy as np
import astropy
from astropy import time, coordinates as coord, units as u

def get_bary(toas,source,location):
    times = time.Time([toas],format='mjd',scale='utc',location=location)
    ltt_bary = times.light_travel_time(source) #time offset between barycentre and Earth
    toas_bary = times.tdb + ltt_bary
    return toas_bary.value[0]

def barycorr(obs_start, burst_time, f_ref, dm, FRB='R1', telescope='Eff'):
    """
    obs_start is the start time of the scan in MJD.
    burst_time is the time of the burst relative to the obs_start in seconds.
    f_ref is the reference frequency to correct to infinite frequency.
    dm is the dispersion measure.
    FRB_RA is the right ascension of the FRB in the format 01:58:00.7502
    FRB_dec is the declination of the FRB in the format +65:43:00.3152
    telescope  is the telescope you used for your observation (Eff, CHIME, DSS43)
    """
    burst_time_MJD=burst_time/(24.*3600.)

    #obs_start from the filterbank header
    dm_shift = (4150.*(1./(f_ref)**2)*dm)/(24.*3600.) #top of the top frequency channel (readfile reports mid of channel) #in MJD
    FRB=str(FRB)
    if FRB == 'R3':
        FRB_coord=coord.SkyCoord("01:58:00.7502", "+65:43:00.3152",unit=(u.hourangle, u.deg), frame='icrs') #R3 obs pos
    if FRB == 'R1':
        FRB_coord=coord.SkyCoord("05:31:58.600", "+33:08:49.60",unit=(u.hourangle, u.deg), frame='icrs')
    

    telescope=str(telescope)
    if telescope=='Eff':
        telescope_coord = coord.EarthLocation.from_geodetic(lon = (06. +53./60. + 00.99425/3600.)*u.deg,lat=(50. + 31./60.+29.39459/3600.)*u.deg,height=369.082*u.m) #effelsberg geodetic coords in deg
    if telescope=='CHIME':
        telescope_coord = coord.EarthLocation.from_geodetic(lon = (-119. +37./60. + 26./3600.)*u.deg,lat=(49. + 19./60.+16./3600.)*u.deg,height=545.0*u.m) #CHIME geodetic coords in deg
    if telescope=='DSS43':
        telescope_coord = coord.EarthLocation.from_geodetic(lon = (148. +58./60. +  52.55394/3600.)*u.deg,lat=(35. +24./60. + 8.74388/3600.)*u.deg,height=689.608*u.m) #DSS-43 geodetic coords in deg 

    start=obs_start
    TOA = burst_time
    TOA_correctDM = (TOA-dm_shift)+start
    TOA_bary = get_bary(TOA_correctDM, source=FRB_coord, location=telescope_coord)
    return TOA_bary
