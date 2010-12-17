#!/usr/bin/env mypython
"""
make a polar plot from beamformer output
"""

from pylab import *
import scipy.io as sio
from matplotlib import cm
from obspy.signal.rotate import gps2DistAzimuth
from collections import defaultdict
from matplotlib import rcParams
import sys
import os
import glob
sys.path.append(os.path.join(os.environ['PROC_SRC'],'misc'))
import delaz
from monthdict import monthdict
rcParams = {'backend':'Agg'}
import datetime

#fig = figure(figsize=(6, 6))
#ax = fig.add_subplot(1,1,1,projection='polar')
cmap = cm.get_cmap('jet')
if 0:
    ### array center of gravity
    meanlon = 175.5851
    meanlat = -39.1686
    sites = loadtxt('SiteLocations.dat')
    data1 = loadtxt('22February2001.SWH.dat')
    #data1 = loadtxt('03March2001.SWH.dat')
    swh = []
    for _l in xrange(sites.shape[0]):
        dist,az,baz = gps2DistAzimuth(meanlat,meanlon,sites[_l,2],sites[_l,1])
        swh.append([az,data1[_l,1]])
    swh = array(swh)
    theta= arange(0,362,2)
    swh_grid = zeros(theta.size)
    idx = theta.searchsorted(swh[:,0])
    cnt = 0
    for _i in idx:
        if swh[cnt,1] > swh_grid[_i]:
            swh_grid[_i] = swh[cnt,1]
        #swh_grid[_i] += swh[cnt,1]
        cnt += 1
    swh_grid = swh_grid.reshape(swh_grid.size,1)
    sl = array([0.5,0.6])
    ax.contourf((theta[::-1]+90.)*pi/180.,sl,hstack((swh_grid,swh_grid)).T,100,cmap=cmap)
    

if 1:
    fig = figure(figsize=(6, 6))
    ax = fig.add_subplot(1,1,1,projection='polar')
    matfile = '6s_average_beam.mat'
    matfile = 'all_beam_8s.mat'
    #matfile = '/home/data/dev/beamforming/laura/all_stations/8s/2001_6_6_0_0_0/beam88.mat'
    #_d = '/home/data/dev/beamforming/laura/all_stations/8s/2001_6_6_0_0_0'
    #year = int(os.path.basename(_d).split('_')[0])
    #month = int(os.path.basename(_d).split('_')[1])
    #mday = int(os.path.basename(_d).split('_')[2])
    #odir = os.path.join('/data/wanakaII/yannik/start/sacfiles',
    #                        str(year),monthdict[month],os.path.basename(_d))
    #fl = glob.glob(os.path.join(odir,'ft_grid*.HHZ.SAC'))
    #tre = sio.loadmat(matfile)['beam']
    a = sio.loadmat(matfile)['beam']
    tre = a.mean(axis=0)
    #tre = 10*log10(tre)
    #tre -= tre.max()
    #beam = sio.loadmat('beam62.mat')['beam']
    #tre = squeeze(beam[:,13,:,:])
    #tre = tre.mean(axis=2)
    #tre = 10*log10(tre)
    #tre = tre-tre.max()
    
    theta= arange(0,362.,2.)
    slowness=arange(0.03,0.505,0.005)  ###slowness in s/km
    ax.contourf((theta[::-1]+90.)*pi/180.,slowness[31:],tre[:,31:].T,100,cmap=cmap,antialiased=True,
                linewidths=0.1,linstyles='dotted')
    ax.contour((theta[::-1]+90.)*pi/180.,slowness[31:],tre[:,31:].T,100,cmap=cmap)
    #ax.set_thetagrids([0,45.,90.,135.,180.,225.,270.,315.],
    #                  labels=['90','45','0','315','270','225','180','135'])
    ### add tvz wedge
    lon0 = 175.52
    lat0 = -38.90
    lon1 = 176.10
    lat1 = -38.87
    ### mean start coordinates:
    mlat = -39.182
    mlon = 175.60
    delta1,az1,baz1 = delaz.delaz(mlat,mlon,lat0,lon0,0)
    delta1,az2,baz2 = delaz.delaz(mlat,mlon,lat1,lon1,0)
    ax.plot([-(az1-90)*pi/180.,-(az1-90)*pi/180.],[0.,0.5],'k')
    ax.plot([-(az2-90)*pi/180.,-(az2-90)*pi/180.],[0.,0.5],'k')

    ax.set_rmax(0.5)
    ax.set_thetagrids([0,45.,90.,135.,180.,225.,270.,315.],
                      labels=[])
    ax.set_rgrids([0.2,0.3,0.4,0.5],labels=['0.2','0.3','0.4','s/km'],color='r')
    ax.set_title("8 s")
    ax.grid(True)
    show()
    savefig('average_beam_8s_start.pdf')

if 0:
    ### evaluate beamformer output 
    dl = glob.glob('/home/data/dev/beamforming/laura/all_stations/10s/2001*')
    dates = array([])
    azims = array([])
    PLOT = False
    for _d in dl:
        year = int(os.path.basename(_d).split('_')[0])
        month = int(os.path.basename(_d).split('_')[1])
        mday = int(os.path.basename(_d).split('_')[2])
        odir = os.path.join('/data/wanakaII/yannik/start/sacfiles',
                            str(year),monthdict[month],os.path.basename(_d))
        fl = glob.glob(os.path.join(odir,'ft_grid*.HHZ.SAC'))
        if len(fl) < 10: continue
        matfile = os.path.join(_d,'beam88.mat')
        if not os.path.isfile(matfile):
            print matfile, ' does not exist'
            continue
        beam = sio.loadmat(matfile)['beam']
        tre = squeeze(beam[:,10,:,:])
        tre = tre.mean(axis=2)
        tre = tre-tre.max()
        theta= arange(0,362.,2.)
        slowness=arange(0.03,0.505,0.005)  ###slowness in s/km
        idx = where((slowness>0.2) & (slowness<0.4))
        new_tre = tre[:,idx[0]]
        id1,id2 = unravel_index(new_tre.argmax(),new_tre.shape)
        if PLOT:
            fig = figure(figsize=(6, 6))
            ax = fig.add_subplot(1,1,1,projection='polar')
            cmap = cm.get_cmap('jet')
            ax.contourf((theta[::-1]+90.)*pi/180.,slowness[14:],tre[:,14:].T,
                        100,cmap=cmap,antialiased=True,
                        linewidths=0.1,linstyles='dotted')
            ax.set_thetagrids([0,45.,90.,135.,180.,225.,270.,315.],
                              labels=['90','45','0','315','270','225','180','135'])
            ax.set_rgrids([0.1,0.2,0.3,0.4,0.5],labels=['0.1','0.2','0.3','0.4','0.5'],color='r')
            ax.grid(True)
            ax.set_title(os.path.basename(_d))
            ax.plot(-(theta[id1]-90)*pi/180.,slowness[idx][id2],'ko')
            ax.set_rmax(0.5)
        print _d,theta[id1],slowness[idx][id2]
        dates = append(dates,datetime.date(year,month,mday))
        azims = append(azims,theta[id1])
        show()

    sidx = argsort(dates)
    plot(dates[sidx],azims[sidx])
    plot(dates[sidx],azims[sidx],'ro')
    ax = gca()
    labels = ax.get_xticklabels()
    setp(labels, 'rotation', 45, fontsize=10)
    months   = MonthLocator()
    monFmt   = DateFormatter(' %b %y')
    ax.xaxis.set_major_locator(months)
    ax.xaxis.set_major_formatter(monFmt)
    #savefig('beam_max_vs_azimuth_10s.pdf')
    #savetxt('beam_max_vs_azimuth_10s.txt',vstack((dates[sidx],azims[sidx])).T,fmt="%s %f")


if 0:
    from obspy.core import UTCDateTime
    def mycmp(fn1,fn2):
        year = int(os.path.basename(fn1).split('_')[0])
        month = int(os.path.basename(fn1).split('_')[1])
        mday = int(os.path.basename(fn1).split('_')[2])
        a = UTCDateTime(year,month,mday)
        year = int(os.path.basename(fn2).split('_')[0])
        month = int(os.path.basename(fn2).split('_')[1])
        mday = int(os.path.basename(fn2).split('_')[2])
        b = UTCDateTime(year,month,mday)
        if a<b:
            return -1
        if a>b:
            return 1
        if a==b:
            return 0
        
    ### compress beamformer output into one big matrix
    dl = glob.glob('/home/data/dev/beamforming/laura/all_stations/8s/2001*')
    dl_new = sorted(dl,cmp=mycmp)
    all_beam = zeros((len(dl),181,95))
    dates = zeros(len(dl))
    cnt = 0
    for _d in dl_new:
        matfile = os.path.join(_d,'beam88.mat')
        if not os.path.isfile(matfile):
            print matfile," does not exist"
            continue
        print matfile
        beam = sio.loadmat(matfile)['beam']
        tre = squeeze(beam[:,13,:,:])
        tre = tre.mean(axis=2)
        tre = tre-tre.max()
        all_beam[cnt,:,:] = tre
        year = int(os.path.basename(_d).split('_')[0])
        month = int(os.path.basename(_d).split('_')[1])
        mday = int(os.path.basename(_d).split('_')[2])
        a = UTCDateTime(year,month,mday)
        dates[cnt] = a.julday
        cnt += 1
    sio.savemat('all_beam_8s.mat',{'beam':all_beam,'dates':dates})
    
