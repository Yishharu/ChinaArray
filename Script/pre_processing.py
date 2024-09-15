from glob import glob
import subprocess
import numpy as np
import os
import warnings
warnings.filterwarnings("ignore")

from obspy.core import Stream
from obspy import read, read_inventory, read_events

# from geographiclib.geodesic import Geodesic
from obspy.geodetics.base import gps2dist_azimuth
from obspy.geodetics import kilometers2degrees
# import obspy
from obspy.signal.rotate import rotate_ne_rt
from obspy.taup import TauPyModel
#Calculate the takeoff angle of the ray
model = TauPyModel('prem')

import time
from multiprocessing import Pool

nproc = 4
pre_filt = [0.005, 0.01, 5, 10]
MseedDataDirectory = "../Data/fetch_fdsn_lizhi"
ResponseDirectory = "../Data/fetch_fdsn_sxd/response"

StationCatalog = dict()
with open("../Data/fetch_fdsn_sxd/station.meta", 'r') as f:
    for line in f:
        StationName = line.split('|')[0] + '.' + line.split('|')[1]
        StationCatalog[StationName] = dict()
        StationCatalog[StationName]['NetworkCode'] = line.split('|')[0]
        StationCatalog[StationName]['StationCode'] = line.split('|')[1]
        StationCatalog[StationName]['Channel'] = line.split('|')[3]
        StationCatalog[StationName]['Latitude'] = float(line.split('|')[4])
        StationCatalog[StationName]['Longitude'] = float(line.split('|')[5])

# for loop through all mseed file
# for MseedFilePath in glob(MseedDataDirectory+'/*.mseed'):
def ProcessMseed(MseedFilePath):
    EVENTNAME = MseedFilePath.split('/')[-1].split('.')[0]
    if os.path.exists(f"../Data/ProcessedSecondRequest/{EVENTNAME}.PICKLE"):
        print(f"{EVENTNAME} PICKLE exists, skipped!!")
        return
    DataStream = read(MseedFilePath,format='MSEED')
    print(EVENTNAME)

    CatFind = glob(f"../Data/CMTSOLUTION/{EVENTNAME[0:12]}*.CMTSOLUTION")
    if len(CatFind)==1:
        cat = read_events(CatFind[0])
    elif len(CatFind)>1: 
        print("Multiple CMTSOLUTION (>1) found!!!")
        return        
    else:
        print("CMTSOLUTION not found!!!")
        return
    SourceLat = cat[0].origins[0].latitude
    SourceLon = cat[0].origins[0].longitude
    SourceDepth = cat[0].origins[0].depth/1.0e3


    # trim DataStream for same length
    DataStream.trim(starttime=DataStream[0].stats.starttime, endtime=DataStream[0].stats.starttime + 2400)
    # calculate baz for rotating
    for itrace, trace in enumerate(DataStream.select(component="N")):
        StationName = DataStream[itrace].stats.network + "." + DataStream[itrace].stats.station
        distm, azimuth, backazimuth = gps2dist_azimuth(SourceLat,SourceLon,StationCatalog[StationName]['Latitude'],StationCatalog[StationName]['Longitude'])
        distance_in_degrees = kilometers2degrees(distm/1.0e3)
        
        trace.stats.distance = distm
        trace.stats.distance_in_degrees = distance_in_degrees
        trace.stats.azimuth = azimuth
        trace.stats.backazimuth = backazimuth

        # if not hasattr(trace.stats,'traveltimes'): # add traveltime
        trace.stats.traveltimes=dict()
        arrivals = model.get_travel_times(source_depth_in_km = SourceDepth, distance_in_degree = distance_in_degrees,
                                            phase_list = ['Sdiff','S'], receiver_depth_in_km = 0.)
        for arrival in arrivals:
            trace.stats.traveltimes[arrival.name]=arrival.time

    RTZDataStream = Stream()

    for itrace, trace in enumerate(DataStream.select(component="N")):
        TChannelName = trace.stats.network + "." + trace.stats.station + "." \
                    + trace.stats.location + "." "BHT"
        seisCheck = RTZDataStream.select(id=TChannelName)
        if len(seisCheck)>0:
            print('processed trace exists!!!')
            continue

        # print(trace)
        EChannelName = trace.stats.network + "." + trace.stats.station + "." \
                    + trace.stats.location + "." "BHE"
        seisE = DataStream.select(id=EChannelName)
        
        NChannelName = trace.stats.network + "." + trace.stats.station + "." \
                    + trace.stats.location + "." "BHN"
        seisN = DataStream.select(id=NChannelName)

        ZChannelName = trace.stats.network + "." + trace.stats.station + "." \
                    + trace.stats.location + "." "BHZ"
        seisZ = DataStream.select(id=ZChannelName)

        if len(seisE) != 1 or len(seisN) != 1:
            # Merge segaments into one single trace
            seisE.merge(method=1,fill_value='interpolate',interpolation_samples=2)
            seisN.merge(method=1,fill_value='interpolate',interpolation_samples=2)
            seisZ.merge(method=1,fill_value='interpolate',interpolation_samples=2)
        

        # Remove response
        try:
            invE = read_inventory(ResponseDirectory+"/RESP."+EChannelName)
            seisE.remove_response(inventory=invE, pre_filt=pre_filt, output="DISP") 
            invN = read_inventory(ResponseDirectory+"/RESP."+NChannelName)
            seisN.remove_response(inventory=invN, pre_filt=pre_filt, output="DISP")
            invZ = read_inventory(ResponseDirectory+"/RESP."+ZChannelName)
            seisZ.remove_response(inventory=invZ, pre_filt=pre_filt, output="DISP")
        except:
            print(trace, "remove response failed")
            continue

        if len(seisN[0].data) != len(seisE[0].data): # check length
            seisN.resample(100)
            seisE.resample(100)
            trimstart = max(seisN[0].stats.starttime, seisE[0].stats.starttime)
            trimend = min(seisN[0].stats.endtime, seisE[0].stats.endtime)
            seisN.trim(trimstart, trimend)
            seisE.trim(trimstart, trimend)
        [seisRtmp,seisTtmp] = rotate_ne_rt(seisN[0].data, seisE[0].data, seisN[0].stats.backazimuth)
        
        seisR=seisN[0].copy()
        seisR.stats['channel']='BHR'
        seisR.data=seisRtmp
        seisT=seisN[0].copy()
        seisT.stats['channel']='BHT'
        seisT.data=seisTtmp

        seisZ[0].stats.traveltimes = seisT.stats.traveltimes
        seisZ[0].stats.distance = seisT.stats.distance
        seisZ[0].stats.distance_in_degrees = seisT.stats.distance_in_degrees
        seisZ[0].stats.azimuth = seisT.stats.azimuth
        seisZ[0].stats.backazimuth = seisT.stats.backazimuth
        
        RTZDataStream += seisR
        RTZDataStream += seisT
        RTZDataStream += seisZ[0]

    RTZDataStream.resample(10)
    RTZDataStream.write(f"../Data/ProcessedSecondRequest/{EVENTNAME}.PICKLE",format='PICKLE')
    print(f"../Data/ProcessedSecondRequest/{EVENTNAME}.PICKLE Saved!!!")


with Pool(nproc) as p:
    p.map(ProcessMseed,glob(MseedDataDirectory+'/*.mseed'))  # Multiprocessing DownloadEvent