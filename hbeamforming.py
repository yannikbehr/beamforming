#!/usr/bin/env mypython
###/usr/local/python2/bin/python
"""
next try to rewrite laura's beamformer
"""

import os
import sys
import glob
import obspy.sac
from obspy.sac import *
from obspy.core import read
from pylab import *
import obspy.signal
import scipy.io as sio
from matplotlib import cm, rcParams
import ctypes as C
import pickle
sys.path.append(os.path.join(os.environ['PROC_SRC'],'misc'))
rcParams = {'backend':'Agg'}

DEBUG = True
def prep_beam_h(files,matfile,nhours=1,fmax=10.,fact=10,new=True):
    if new:
        ntimes = int(round(24/nhours))
        step = nhours*3600*fmax/fact

        slons = array([])
        slats = array([])
        nfiles = len(files)
        seisbandn = zeros((nfiles,ntimes,step))
        seisbande = zeros((nfiles,ntimes,step))
        sigmas = []
        for i,_ftup in enumerate(files):
            trn = read(_ftup[0])[0]
            tre = read(_ftup[1])[0]
            slons = append(slons,trn.stats.sac.stlo)
            slats = append(slats,trn.stats.sac.stla)
            trn.downsample(decimation_factor=fact, strict_length=True)
            tre.downsample(decimation_factor=fact, strict_length=True)
            nptsn = trn.stats.npts
            nptse = tre.stats.npts
            df = trn.stats.sampling_rate
            dt = trn.stats.delta
            seis0n = zeros(24*3600*int(df))
            seis0e = zeros(24*3600*int(df))
            seis0n[0:nptsn] = trn.data
            seis0e[0:nptse] = tre.data
            seis0n -= seis0n.mean()
            seis0e -= seis0e.mean()
            for j in xrange(ntimes):
                ilow = j*step
                iup = (j+1)*step
                seisbandn[i,j,:] = seis0n[ilow:iup]
                seisbandn[i,j,:] -= seisbandn[i,j,:].mean()
                seisbande[i,j,:] = seis0e[ilow:iup]
                seisbande[i,j,:] -= seisbande[i,j,:].mean()

        if 1:
            fftpower = 7
            ismall = 2**fftpower
            ipick = arange(ismall)
            n=nhours*3600*df
            nsub = int(np.floor(n/ismall)) # Number of time pieces -20 mins long each
            #seissmall = zeros((len(ipick),ntimes,nsub,nfiles))
            seissmalln = zeros((nfiles,ntimes,nsub,len(ipick)))
            seissmalle = zeros((nfiles,ntimes,nsub,len(ipick)))
            for ii in xrange(nfiles):
                for jj in xrange(ntimes):
                    for kk in xrange(nsub):
                        #seissmall[:,jj,kk,ii] = seisband[kk*ismall+ipick,jj,ii]
                        seissmalln[ii,jj,kk,:] = seisbandn[ii,jj,kk*ismall+ipick]
                        seissmalle[ii,jj,kk,:] = seisbande[ii,jj,kk*ismall+ipick]

        LonLref= 165
        LonUref= 179.9
        LatLref= -48
        LatUref= -34
        stacoord=vstack((slons,slats))
        ##Find the stations which belong to this grid
        idx = where((stacoord[0] >= LonLref) & (stacoord[0] <= LonUref) & \
                    (stacoord[1] >= LatLref) & (stacoord[1] <= LatUref))
        #nb. this simple 'mean' calc only works if we don't cross lat=0 or lon=180
        meanlat = slats.mean()
        meanlon = slons.mean()

        sio.savemat(matfile,{'seissmalln':seissmalln,'seissmalle':seissmalle,'slats':slats,'slons':slons,'dt':dt})
        return seissmalln, seissmalle, meanlat, meanlon, slats, slons, dt
    else:
        a = sio.loadmat(matfile)
        seissmalln = a['seissmalln']
        seissmalle = a['seissmalle']
        slats = a['slats']
        slons = a['slons']
        dt = a['dt'][0][0]
        slats = slats.reshape(slats.shape[0],)
        slons = slons.reshape(slons.shape[0],)
        meanlat = slats.mean()
        meanlon = slons.mean()
        return seissmalln, seissmalle, meanlat, meanlon, slats, slons, dt


def calc_steer(slats,slons):
    theta= arange(0,362,2)
    theta = theta.reshape((theta.size,1))
    sta_origin_dist = array([])
    sta_origin_bearing = array([])
    for lat,lon in zip(slats,slons):
        dist, az, baz = obspy.signal.rotate.gps2DistAzimuth(meanlat,meanlon,lat,lon)
        sta_origin_dist = append(sta_origin_dist,dist)
        sta_origin_bearing = append(sta_origin_bearing,az)
    sta_origin_x = sta_origin_dist*cos(sta_origin_bearing*pi/180.)
    sta_origin_y = sta_origin_dist*sin(sta_origin_bearing*pi/180.)
    #zeta_x = sta_origin_dist*cos(sta_origin_bearing*pi/180.)
    #zeta_y = sta_origin_dist*sin(sta_origin_bearing*pi/180.)
    zeta_x = -cos(theta*pi/180.)
    zeta_y = -sin(theta*pi/180.)
    #dot product betwen zeta and x
    zetax = zeta_x*sta_origin_x+zeta_y*sta_origin_y
    #slowness in s/km
    slowness = arange(0.03,0.505,0.005)
    slowness = slowness.reshape((1,slowness.size))
    return zetax,theta,slowness,sta_origin_x,sta_origin_y


def mkfilelist(filesN, filesE):
    """
    Order the file lists for the north and east component.
    """
    newlist = []
    for _fN in filesN:
        for _fE in filesE:
            a = os.path.basename(_fN).split('_')
            stN = a[2].split('.')[0]
            a = os.path.basename(_fE).split('_')
            stE = a[2].split('.')[0]
            if stN == stE:
                newlist.append((_fN,_fE))
                break
            
    return newlist

def rotate_fft(seisbandn,seisbande,az):
    """
    Rotate traces into radial direction and calculate fft of radial trace.
    """
    az *= pi/180.
    r = cos(az)*seisbandn + sin(az)*seisbande
    t = -sin(az)*seisbandn + cos(az)*seisbande
    R = fft(r,axis=2)
    T= fft(t,axis=2)
    return R,T

def beamforming(seisn,seise,slowness,zetax,theta,dt,new=True,matfile=None,freq_int=(0.02,0.4)):
    if new:
        _p = 6.
        _f = 1./_p
        nstat, ntimes, nsub, nfft = seisn.shape
        nsources = theta.size
        freq = fftfreq(nfft,dt)
        I = np.where((freq>freq_int[0]) & (freq<freq_int[1]))
        beam = zeros((nsources,slowness.size,ntimes))
        df = dt/nfft
        ind = int(_f/df)
        print ind, freq[ind]
        N = fft(seisn,n=nfft,axis=3)
        E = fft(seise,n=nfft,axis=3)
        for i,az in enumerate(theta):
            daz = az*pi/180.
            R = cos(daz)*N + sin(daz)*E
            T = -sin(daz)*N + cos(daz)*E
            dist = zetax[i,:]
            for ww in [ind]:
                FF = freq[ww]
                omega = 2*pi*FF
                for cc in xrange(slowness.shape[1]):
                    velocity = 1./slowness[0][cc]*1000
                    e = exp(-1j*dist*omega/velocity)
                    eT = e.T.copy()
                    #for tt in xrange(ntimes):
                    for tt in [0]:
                        #for TT in xrange(nsub):
                        for TT in [0]:
                            Y = asmatrix(squeeze(T[:,tt,TT,ww]))
                            #Y = squeeze(R[:,tt,TT,ww])
                            YT = Y.T.copy()
                            cov = dot(YT,conjugate(Y))
                            #import ipdb
                            #ipdb.set_trace()
                            #r,t = rotate_fft(seisn,seise,az)
                            beam[i,cc,tt] = abs(asarray(dot(conjugate(eT),dot(cov,e).T)))**2
                            #beam[i,cc,tt] = abs(dot(Y,conjugate(e)))**2

        #sio.savemat(matfile,{'beam':beam})
    else:
        beam = sio.loadmat(matfile)['beam']
    return beam


def polar_plot(beam,theta,slowness):
    theta = theta[:,0]
    slowness = slowness[0,:]
    tre = squeeze(beam[:,:,0])
    tre = tre-tre.max()
    fig = figure(figsize=(6,6))
    ax = fig.add_subplot(1,1,1,projection='polar')
    cmap = cm.get_cmap('jet')
    ax.contourf((theta[::-1]+90.)*pi/180.,slowness[14:],tre[:,14:].T,
                100,cmap=cmap,antialiased=True,
                linewidths=0.1,linstyles='dotted')
    ax.contour((theta[::-1]+90.)*pi/180.,slowness[14:],tre[:,14:].T,
                100,cmap=cmap)
    ax.set_thetagrids([0,45.,90.,135.,180.,225.,270.,315.],
                      labels=['90','45','0','315','270','225','180','135'])
    ax.set_rgrids([0.1,0.2,0.3,0.4,0.5],labels=['0.1','0.2','0.3','0.4','0.5'],color='r')
    ax.grid(True)
    #ax.set_title(os.path.basename(_d))
    ax.set_rmax(0.5)
    show()


if __name__ == '__main__':
    datdir = '/Volumes/Wanaka_01/yannik/start/sacfiles/10Hz/2001/Feb/2001_2_22_0_0_0/'
    datdir = '/Volumes/Wanaka_01/yannik/start/sacfiles/10Hz/2001/Mar/2001_3_3_0_0_0/'
    filesN = glob.glob(os.path.join(datdir,'ft_grid*.HHN.SAC'))
    filesE = glob.glob(os.path.join(datdir,'ft_grid*.HHE.SAC'))
    matfile = 'prep_beam_h_2001_3_3.mat'
    nlist = mkfilelist(filesN, filesE)
    seisn, seise, meanlat, meanlon, slats, slons, dt = prep_beam_h(nlist,matfile,nhours=1,fmax=10.,fact=10,new=False)
    zetax,theta,slowness,sta_origin_x,sta_origin_y = calc_steer(slats,slons)
    matfile = 'beam_h_2001_3_3.mat'
    beam = beamforming(seisn,seise,slowness,zetax,theta,dt,new=True,freq_int=(0.01,0.4),
                       matfile=matfile)
    polar_plot(beam,theta,slowness)
    
