#!/raid1/zl382/Libs/anaconda3/bin/python3

import obspy
from obspy import read_events, read_inventory
from obspy.clients.fdsn import Client
from obspy import UTCDateTime
import matplotlib.pyplot as plt
import os.path
import time
from obspy.geodetics import kilometers2degrees
# from geographiclib.geodesic import Geodesic 
from obspy.geodetics.base import gps2dist_azimuth

import numpy as np
import shutil
import requests
from bs4 import BeautifulSoup
import sys
from multiprocessing import Pool
from obspy.core import Stream
from obspy import read

from obspy.taup import TauPyModel

nproc = 4

EventListName = 'ChinaArray'

# SearchCenterLat = -15.9650    # Easter Search Lat -30.192 Lon -110.589; Hawaii Lat 15.4 Lon -172.3;
# SearchCenterLon = -5.7089  # Macdonald Search Lat -29 Lon -140.3;
SearchCenterLat, SearchCenterLon = 30.5928, 114.3052

MinRadius = 60  # 30-70 deg affect on the Sdiff path, 40-60 effects are dominant
MaxRadius = 180  # 30-70 deg affect on the Sdiff path, 40-60 effects are dominant
MinDepth = 0  # avoid shallow events for clean stf
MaxDepth = 1000
MinMagnitude = 6.0
MaxMagnitude = 8.0
magnitudetype = 'Mww'
OrderBy = 'time-asc'
StartTime = UTCDateTime("2010-01-01T00:00:00")
EndTime = UTCDateTime.utcnow()

# define the station geometry
StaDistMin = 90
StaDistMax = 140
AziRange = 40 # half azimuth range
phase_list = ['Sdiff','S'] # calculate ref phase time
model = TauPyModel(model='prem')

# define the waveform parameter
bandcode = "BH*"
second_before_phase = 200.
second_after_phase = 200.
# define a filter band to prevent amplifying noise during the deconvolution
# not being used, as responses are not being removed
fl1 = 0.005
fl2 = 0.01
fl3 = 5.
fl4 = 10.

# load IRIS client
irisclient = Client("IRIS") 
cat = irisclient.get_events(starttime=StartTime, endtime=EndTime, latitude=SearchCenterLat, longitude=SearchCenterLon,
                            minradius=MinRadius, maxradius=MaxRadius, mindepth=MinDepth, maxdepth=MaxDepth,
                            minmagnitude=MinMagnitude, maxmagnitude=MaxMagnitude, orderby=OrderBy, magnitudetype=magnitudetype) 
# cat = irisclient.get_events(starttime=StartTime, endtime=EndTime, mindepth=MinDepth, maxdepth=MaxDepth,
#                             minmagnitude=MinMagnitude, maxmagnitude=MaxMagnitude, orderby=OrderBy, magnitudetype=magnitudetype) 


print(len(cat), "events found!")
# Make Request List
with open('../Log/REQUEST_Catalog_%s.txt' %EventListName, 'w+') as f:
    f.write(cat.__str__(print_all=True))          

# Make Saved List
OUTPUTFILE = '../Log/SAVED_Catalog_%s.txt' %EventListName 
with open(OUTPUTFILE, 'w+') as f:
    line = '%%%%%%%%%%%%%%%%%% This is the Event Output Result %%%%%%%%%%%%%%%%%%%%%%\n'
    f.write(line)


def DownloadEvent(EVENT):
# for EVENT in cat:
    TimeString = str(EVENT.origins[0].time)
    # EventName = TimeString[0:4] + TimeString[5:7] + TimeString[8:10]
    # Extend to second to avoid repeated EQs
    EventName = TimeString[0:4] + TimeString[5:7] + TimeString[8:10] + TimeString[11:13] + TimeString[14:16] + TimeString[17:19]
    EventMag = EVENT.magnitudes[0]['mag']
    EvtLat = EVENT.origins[0]['latitude']
    EvtLon = EVENT.origins[0]['longitude']
    EvtDep = EVENT.origins[0]['depth']/1.e3 # convert to km from m
    EvtStartTime = EVENT.origins[0].time
    EvtEndTime = EvtStartTime + 2000.

    SearchAzimuth =  gps2dist_azimuth(EvtLat, EvtLon, SearchCenterLat, SearchCenterLon)[1]
    SearchAzimuth = (SearchAzimuth + 360)%360

    EventDir='../Data/%s/' %(EventListName)
    FilePathName = '%sOriginal/%s.PICKLE' %(EventDir, EventName)
    # Make directories for data

    os.makedirs(EventDir, exist_ok=True)
    # os.makedirs(EventDir+'Original/', exist_ok=True)

    # else:
    #     print("%s already exists" %EventDir)

    # Make CMTSOLUTION file
    timeA = EvtStartTime
    year, month, day = timeA.year, timeA.month, timeA.day
    max_dt = 1*60.
    minlat, maxlat = max(EvtLat-5, -90), min(EvtLat+5, 90)        #latitude range [-90,90]
    minlon, maxlon = max(EvtLon-5, -180), min(EvtLon+5, 180)      #longitude range [-180,180]
    mindepth, maxdepth = max(EvtDep-100,0), EvtDep+100    #depth range [0,1000]
    # can also add parameter limits for magnitude but IRIS and GLOBALCMT often conflict        
    url = "https://www.globalcmt.org/cgi-bin/globalcmt-cgi-bin/CMT5/form?itype=ymd&yr="+str(year)+"&mo="+str(month)+"&day="+str(day)+"&oyr=1976&omo=1&oday=1&jyr=1976&jday=1&ojyr=1976&ojday=1&otype=nd&nday=1&lmw="+str(EventMag-1)+"&umw="+str(EventMag+1)+"&lms=0&ums=10&lmb=0&umb=10&llat="+str(minlat)+"&ulat="+str(maxlat)+"&llon="+str(minlon)+"&ulon="+str(maxlon)+"&lhd="+str(mindepth)+"&uhd="+str(maxdepth)+"&lts=-9999&uts=9999&lpe1=0&upe1=90&lpe2=0&upe2=90&list=4"
    res = requests.get(url)
    html_page = res.content
    soup = BeautifulSoup(html_page, 'html.parser')
    try:
        text = soup.find_all("pre")[1]
    except:
        print('CMT Webpage not working!')
        with open(OUTPUTFILE, 'a+') as f:
            line = 'CMT Webpage not working!'
            f.write(line)
        return

    if len(text.contents)==4:
        parsed = text.contents[3].split("\n\n")[0:-1]
        parsed[0] = parsed[0][1:]
    else:
        parsed = text.contents[0].split("\n\n")[1:-1]

    hr, mn, sc, lt, ln, dp = [], [], [], [], [], []
    output = []
    
    #print(len(parsed))
    for i in range(len(parsed)):
        line = parsed[i].split("\n")
        hr.append(int(line[0][16:18]))
        mn.append(int(line[0][19:21]))
        sc.append(int(float(line[0][22:27])))
        lt.append(float(line[0][28:36]))
        ln.append(float(line[0][37:46]))
        dp.append(float(line[0][47:52]))
        timeB = UTCDateTime(year,month,day,hr[-1],mn[-1],min(sc[-1],59))
        delta = abs(timeA-timeB)  
    
        # check for an event under the allowed time difference
        if delta < max_dt: output.append(parsed[i])
        #print('==>',i,delta,hr[-1],mn[-1],sc[-1],lt[-1],ln[-1],dp[-1])
    
    if len(output) == 0:
        print("No matching CMTSOLUTION entries found. Check URL:%s \n" %url)
        with open(OUTPUTFILE, 'a+') as f:
            line = "No matching CMTSOLUTION entries found. Check URL:%s \n" %url
            f.write(line)
        return
    elif len(output) > 1:
        print("Multiple CMTSOLUTION entries found.")
        with open(OUTPUTFILE, 'a+') as f:
            line = "Multiple CMTSOLUTION entries found. \n"
            f.write(line)
        return
    else:
        # add space between PDEW#### and correct to PDEW from PDE if necessary
        output = list(output)
        if output[0][4] == " ": output[0] = output[0][:4]+"W"+output[0][5:]
        output[0] = output[0][:5]+" "+output[0][5:]
        print(output[0])
        ofn = open(EventDir+'/%s.CMTSOLUTION' %EventName,'w')
        ofn.write(output[0])
        ofn.close()
    res.close()






with Pool(nproc) as p:
    p.map(DownloadEvent,cat)  # Multiprocessing DownloadEvent