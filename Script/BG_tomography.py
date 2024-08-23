import sys

sys.path.append('/raid1/sc845/Tomographic_models/LMClust/')
import LMClust_g6_ryb

sys.path.append('/raid1/sc845/Tomographic_models/')
import basemap_circle


import obspy
from obspy import read
from obspy.core import Stream
import obspy.signal
import matplotlib.pyplot as plt
import os.path
import glob
import numpy as np
import scipy
import mpl_toolkits
import mpl_toolkits.basemap
print(mpl_toolkits.basemap.__path__)
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import addcyclic
import matplotlib.image as mpimg
import matplotlib.cm as cm
import subprocess
from matplotlib.colors import LinearSegmentedColormap
from geographiclib.geodesic import Geodesic
from obspy.imaging.beachball import beach


dpi = 300
fig = plt.figure(dpi=dpi,figsize=(10, 10))
ax = fig.add_subplot(111)


LMC=LMClust_g6_ryb.LMClust()
# LMC.read('/raid1/sc845/Tomographic_models/LMClust/','clustgr.txt')
# RGB_light =LinearSegmentedColormap.from_list('rgbmap', LMC.colors/2.+0.5,N=len(LMC.colors))
jet=cm.get_cmap('jet',12)
jet_vals=jet(np.arange(12))+0.4
for i in range(np.shape(jet_vals)[0]):
    for j in range(np.shape(jet_vals)[1]):
        if jet_vals[i,j]>1.:
            jet_vals[i,j]=1.

m = Basemap(projection='ortho',lat_0=19,lon_0=-166,resolution='i')

#m = Basemap(llcrnrlon=lonmin,llcrnrlat=latmin,urcrnrlon=lonmax,urcrnrlat=latmax#,
#                resolution='i',projection='merc',lon_0=27.,lat_0=46.7)
clip_path = m.drawmapboundary()


#m.contour(x,y,layer,cmap=RGB_light,linewidth=0,rasterized=True)   

tmp = np.loadtxt('/raid1/sc845/Tomographic_models/SEMUCB_WM1/UCB_a3d_dist.SEMUCB-WM1.r20151019/SEMUCB_WM1_2800km.dat')
lon = tmp[:,1].reshape((181,361))
lat = tmp[:,2].reshape((181,361))
dvs = tmp[:,3].reshape((181,361))  
s = m.transform_scalar(dvs, lon[0,:], lat[:,0], 1000, 500)
im = m.imshow(s, cmap=plt.cm.seismic_r, clip_path=clip_path, vmin=-10, vmax=10)
# cb = plt.colorbar()
# cb.set_label('dlnVs (%)')
# # cb.set_clim(-10,10)
# cb.remove()
x,y=m(lon,lat)


#m.shadedrelief()
m.drawcoastlines()
m.drawcountries()

plt.savefig('BG_tomography.png',dpi=dpi)