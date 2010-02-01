#!/usr/bin/env mypython

"""
script to calculate array impulse response for a
given geometry
"""

import os
import sys
import glob
import obspy.sac
from obspy.sac import *
from pylab import *
from obspy.signal import rotate

dirn = '/home/data/dev/beamforming/geopsy/dataII'
fl = glob.glob(os.path.join(dirn,'*.SAC'))
slon = []
slat = []
for _f in fl:
    x = ReadSac(_f,headonly=True)
    slon.append(x.GetHvalue('stlo'))
    slat.append(x.GetHvalue('stla'))

meanlat = mean(slat)
meanlon = mean(slon)
#plot(slon,slat, 'bo')
#plot([meanlon],[meanlat],'g*')
#xlabel('Longitude (degrees)')
#ylabel('Latitude (degrees)')
#### read in NZ coordinates
#a = loadtxt('NZ_coastline.dat')
#plot(a[:,0],a[:,1])

theta= arange(0,362,2)
nsources=len(theta)

sta_origin_dist = array([])
sta_origin_bearing = array([])

for lat,lon in zip(slat,slon):
    dist, az, baz = rotate.gps2DistAzimuth(meanlat,meanlon,lat,lon)
    sta_origin_dist = append(sta_origin_dist,dist)
    sta_origin_bearing = append(sta_origin_bearing,az)

sta_origin_x = sta_origin_dist*cos(sta_origin_bearing*pi/180.)
sta_origin_y = sta_origin_dist*sin(sta_origin_bearing*pi/180.)

zeta_x = -cos(theta*pi/180.)
zeta_y = -sin(theta*pi/180.)
zeta_x.resize(len(zeta_x),1)
sta_origin_x.resize(1,len(sta_origin_x))
zeta_y.resize(len(zeta_y),1)
sta_origin_y.resize(1,len(sta_origin_y))
zeta = dot(zeta_x,sta_origin_x) + dot(zeta_y,sta_origin_y)
print zeta.shape
slowness=arange(0.03,0.5,0.005)  ###slowness in s/km
beam = zeros((slowness.size,len(theta)))
R=ones((28,28))
for _s in slowness:
    omega = 2*pi*0.05
    velocity = 1./_s*1000
    e_steer = exp(-1j*zeta*omega/velocity)
    beam[where(slowness==_s),:] = diag(abs(dot(dot(e_steer,R),conj(e_steer).T))**2)

tre = 10*log10(beam)

#polar(theta*pi/180.,beam[0,:])
fig = figure(figsize=(6, 6))
project='polar'
#project='rectilinear'
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], projection=project, axisbg='white')
ax.contourf(theta*pi/180.,slowness,tre,100)
ax.set_rmax(0.5)
grid(True)
#figure()
#pcolor(zetax)


#figure()
#for x,y in zip(sta_origin_x,sta_origin_y):
#    plot([x],[y],'ro')
    
#show()


#f = open('stats.coord').readlines()
#stats = {}
#for _l in f:
#    a = _l.split()
#    if not a[0] in stats.keys():
#        stats[a[0]] = {'x': float(a[1]), 'y': float(a[2]), 'z': float(a[3])}

        
