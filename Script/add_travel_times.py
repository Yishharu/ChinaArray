#!/usr/bin python3
#Usage python3 add_travel_times.py [event] S P Sdiff Pdiff SKS SKKS SKKKS SKKKKS ...  (edit dir below)
import obspy
from obspy import read, read_events
from obspy.core import Stream
from obspy.core import event
#from obspy.taup.taup import getTravelTimes
from obspy import UTCDateTime
import obspy.signal
import matplotlib.pyplot as plt
import os.path
import time
from glob import glob
import shutil
import numpy as np
import scipy
from subprocess import call
import subprocess
#from obspy.taup.taup import getTravelTimes
import sys
from multiprocessing import Pool

nproc = 4
phase=['S','Sdiff']
# phase=[]
# for i in range(2,len(sys.argv)):
#     phase.append(sys.argv[i])

# This code assumes the data is in a Data directory and in PICKLE format. Change here if different. 
PICKLEDataDirectory = "../Data/ProcessedSecondRequest"
# for PICKLEFilePath in glob.glob(PICKLEDataDirectory+'/*.PICKLE')[0:1]:
def processed_traveltime(PICKLEFilePath):
    EVENTNAME = os.path.basename(PICKLEFilePath).split('.PICKLE')[0]
    print(EVENTNAME)
    DataStream = read(PICKLEFilePath,format='PICKLE')
    TraveltimeDataStream = Stream()


    CatFind = glob(f"../Data/CMTSOLUTION/{EVENTNAME[0:12]}*.CMTSOLUTION")
    cat = read_events(CatFind[0])

    SourceLat = cat[0].origins[0].latitude
    SourceLon = cat[0].origins[0].longitude
    SourceDepth = cat[0].origins[0].depth/1.0e3


    for itrace, trace in enumerate(DataStream.select(component="R")):
        StationName = trace.stats.network + '.' + trace.stats.station

        if not hasattr(trace.stats,'traveltimes'):
            trace.stats.traveltimes=dict()
        
        for ph in range(len(phase)):
        #Start travetime dictionary
            test=['taup_time -mod prem -deg '+str(trace.stats.distance_in_degrees)+' -h '+ str(SourceDepth) +' -ph ' + phase[ph]]
            out=subprocess.check_output(test,shell=True,universal_newlines=True) 
            t=out.split()
            print(t)
            l=[x for x in range(len(t)) if t[x]==phase[ph]]
            try:
                time= float(t[l[0]+1])
                trace.stats.traveltimes[phase[ph]]=time
                print(phase[ph],time)
            except:
                time=None
                trace.stats.traveltimes[phase[ph]]=time
        
        # traceR = DataStream.select(id="%s*R" %StationName)[0]
        # traceR.stats.traveltimes = trace.stats.traveltimes

        traceT = DataStream.select(id="%s*T" %StationName)[0]
        traceT.stats.traveltimes = trace.stats.traveltimes
        traceZ = DataStream.select(id="%s*Z" %StationName)[0]
        traceZ.stats.traveltimes = trace.stats.traveltimes


        TraveltimeDataStream += trace
        TraveltimeDataStream += traceT
        TraveltimeDataStream += traceZ


    # Write out seismogram again
    TraveltimeDataStream.write(f"../Data/ProcessedSecondRequest/{EVENTNAME}.PICKLE",format='PICKLE')
    print(f"traveltime saved in ../Data/ProcessedSecondRequest/{EVENTNAME}.PICKLE")

with Pool(nproc) as p:
    p.map(processed_traveltime,glob(PICKLEDataDirectory+'/*.PICKLE'))  # Multiprocessing DownloadEvent