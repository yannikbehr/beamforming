#!/usr/bin/env mypython
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
rcParams = {'backend':'Agg'}

DEBUG = True
def prep_beam(nhours=1,fmax=10.,threshold_std=0.5,fftpower=7,freq_int=(0.02,0.4)):
    datdir = '/Volumes/Wanaka_01/yannik/start/sacfiles/10Hz/2001/Feb/2001_2_22_0_0_0/'
    datdir = '/Volumes/Wanaka_01/yannik/start/sacfiles/10Hz/2001/Mar/2001_3_3_0_0_0/'
    files = glob.glob(os.path.join(datdir,'ft_grid*.HHZ.SAC'))
    ntimes = int(round(24/nhours))
    times = arange(0,24+nhours,nhours)
    step = nhours*3600*fmax
    fs = fmax/step
    freq = fmax/2.*arange(0,2**(fftpower-1))/2**(fftpower-1)

    ### initialise station no.s
    ista = 0
    Ista = 0
    stations = []
    slons = array([])
    slats = array([])
    nfiles = len(files)
    seisband = zeros((nfiles,ntimes,step))
    for i,_f in enumerate(files):
        tr = read(_f)[0]
        staname = tr.stats.station
        kcomp = tr.stats.channel
        slat = tr.stats.sac.stla
        slon = tr.stats.sac.stlo
        slons = append(slons,tr.stats.sac.stlo)
        slats = append(slats,tr.stats.sac.stla)
        df = tr.stats.sampling_rate
        dt = tr.stats.delta
        if staname not in stations:
            stations.append(staname)
        npts = tr.stats.npts
        seis0 = zeros(24*3600*int(df))
        seis0[0:npts] = tr.data
        seis0 -= seis0.mean()
        for j in xrange(ntimes):
            ilow = j*step
            iup = (j+1)*step
            seisband[i,j,:] = seis0[ilow:iup]
            seisband[i,j,:] -= seisband[i,:,j].mean()
            seisband[i,j,:] = sign(seisband[i,j,:])
        fseis = fft(seisband,axis=2)
        if np.isnan(fseis).any():
            print "NaN found"
            return
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

    return fseis, meanlat, meanlon, slats, slons, ntimes, dt


def read_matfiles():
    year=2001
    sacpath = './START_DATA_TEST/'
    matpath = './Matfiles_start_TEST/'
    JulianDay=88
    JD='88'
    beamdict = sio.loadmat('BeamformInputData_test.mat')
    infom = beamdict['infom']
    I = beamdict['I']
    matpath = beamdict['matpath']
    matpath = [u'./Matfiles_test']
    Nsub = beamdict['Nsub']
    Ntimes = beamdict['Ntimes']
    freq = beamdict['freq']
    slons = array([])
    slats = array([])
    for i in infom[0]:
        slons = append(slons,i.slon)
        slats = append(slats,i.slat)

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

    Nfreq = I.size;
    seis1 = zeros((Nfreq,Nsub,Ntimes,idx[0].size),'complex128')
    ic = 0
    for ista in xrange(idx[0].size):
        sta1 = infom[0][ista].staname
        filename = os.path.join(matpath[0],sta1[0]+JD+'.mat')
        if not os.path.isfile(filename):
            print filename, ' did not exist'
            continue
        fseis = sio.loadmat(filename)['fseis']
        seis1[:,:,:,ic] = fseis
        ic += 1
    #frequencies within our range of interest
    freqs = freq[0][I-1]
    return Ntimes, freqs, slons, slats, seis1,meanlat,meanlon


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



def beamforming(seis1,slowness,zetax,nsources,Ntimes,dt,new=True,matfile=None,freq_int=(0.02,0.4)):
    if new:
        _p = 6.
        _f = 1./_p
        freq = fftfreq(seis1.shape[2],dt)
        ind = searchsorted(freq[0:int(seis1.shape[2]/2)],_f)
        I = np.where((freq>freq_int[0]) & (freq<freq_int[1]))
        beam = zeros((nsources,slowness.size,Ntimes))
        print ind, freq[ind]
        for ww in [ind]:
            FF = freq[ww]
            for cc in xrange(slowness.shape[1]):
                omega = 2*pi*FF
                velocity = 1./slowness[0][cc]*1000
                e_steer=exp(-1j*zetax*omega/velocity).T
                beamtemp = empty((len(theta),1))
                beamtemp = None
                for tt in xrange(Ntimes):
                    #for TT in xrange(seis1.shape[1]):
                    Y = asmatrix(squeeze(seis1[:,tt,ww],))
                    R = dot(Y.T,conjugate(Y))
                    if beamtemp is None:
                        beamtemp = atleast_2d(sum(abs(asarray(dot(conjugate(e_steer.T),dot(R,e_steer))))**2,axis=1))
                    else:
                        beamtemp = vstack((beamtemp,atleast_2d(sum(abs(asarray(dot(conjugate(e_steer.T),dot(R,e_steer))))**2,axis=1))))

                    beam[:,cc,tt] = transpose(beamtemp).mean(axis=1)
        sio.savemat(matfile,{'beam':beam})
    else:
        beam = sio.loadmat(matfile)['beam']
    return beam

def arr_resp(freqs,slowness,zetax,theta,sta_origin_x,sta_origin_y,
             new=True,matfile=None,src=False,fout=None,pplot=True):
    """
    calculate array response
    """
    if new:
        theta1 = 90
        zeta_x = -cos(theta1*pi/180.)
        zeta_y = -sin(theta1*pi/180.)
        zeta_src = zeta_x*sta_origin_x + zeta_y*sta_origin_y
        c1 = 3000
        beam = zeros((zetax.shape[0],freqs.size,slowness.size))
        for ww in xrange(freqs.shape[1]):
            FF = freqs[0][ww]
            for cc in xrange(slowness.shape[1]):
                omega = 2*pi*FF
                velocity = 1./slowness[0][cc]*1000
                e_steer=exp(-1j*zetax*omega/velocity).T
                if src:
                    e_src = exp(-1j*zeta_src*omega/c1).T
                    Y = multiply(ones((zetax.shape[1],1),'complex128'),atleast_2d(e_src).T)
                    Y = ones((zetax.shape[1],1),'complex128')
                    R = dot(Y,conjugate(Y).T)
                else:
                    R = ones((zetax.shape[1],zetax.shape[1]))
                beam[:,ww,cc] = atleast_2d(sum(abs(asarray(dot(conjugate(e_steer.T),dot(R,e_steer))))**2,axis=1))
        sio.savemat(matfile,{'beam':beam})
    else:
        beam = sio.loadmat(matfile)['beam']
    theta = theta[:,0]
    slowness = slowness[0,:]
    fig = figure(figsize=(10,10))
    ff = [5,7,10,18]
    if pplot:
        for _r in xrange(4):
            tre = squeeze(beam[:,ff[_r],:])
            tre = tre-tre.max()
            ax = fig.add_subplot(2,2,_r+1,projection='polar')
            cmap = cm.get_cmap('jet')
            ax.contourf((theta[::-1]+90.)*pi/180.,slowness,tre.T,
                        100,cmap=cmap,antialiased=True,
                        linewidths=0.1,linstyles='dotted')
            ax.contour((theta[::-1]+90.)*pi/180.,slowness,tre.T,
                       100,cmap=cmap)
            #ax.set_thetagrids([0,45.,90.,135.,180.,225.,270.,315.],
            #                  labels=['90','45','0','315','270','225','180','135'])
            ax.set_thetagrids([0,45.,90.,135.,180.,225.,270.,315.],
                              labels=[])
            ax.set_rgrids([0.1,0.2,0.3,0.4,0.5],labels=['0.1','0.2','0.3','0.4','s/km'],color='r')
            ax.grid(True)
            ax.set_title(os.path.basename("%.1f s"%(round(1./freqs[0][ff[_r]],1))))
            ax.set_rmax(0.5)
        savefig(fout)
        show()
        

def polar_plot(beam,theta,slowness,resp=False):
    theta = theta[:,0]
    slowness = slowness[0,:]
    if resp:
        tre = squeeze(beam[:,18,:])
    else:
        #tre = squeeze(beam[:,18,:,:])
        tre = beam
        tre = tre.mean(axis=2)
    tre = tre-tre.max()
    idx = where((slowness>0.2) & (slowness<0.4))
    new_tre = tre[:,idx[0]]
    id1,id2 = unravel_index(new_tre.argmax(),new_tre.shape)
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
    if 1:
        if DEBUG:
            print 'preparing raw data'
        fseis, meanlat, meanlon, slats, slons, ntimes, dt = prep_beam()
        if DEBUG:
            print 'calculating steering vector'
        zetax,theta,slowness,sta_origin_x,sta_origin_y = calc_steer(slats,slons)
        if DEBUG:
            print 'beamforming'
        beam = beamforming(fseis,slowness,zetax,theta.size,ntimes,dt,new=True,matfile='6s_average_beam.mat')
        polar_plot(beam,theta,slowness,resp=False)
    if 0:
        Ntimes, freqs, slons, slats, seis,meanlat,meanlon = read_matfiles()
        zetax,theta,slowness,sta_origin_x,sta_origin_y = calc_steer(slats,slons)
        beam = beamforming(seis,freqs,slowness,zetax,theta.size,Ntimes,new=False,matfile='6s_average_beam.mat')
        polar_plot(beam,theta,slowness,resp=True)
    if 0:
        arr_resp(freqs,slowness,zetax,theta,sta_origin_x,sta_origin_y,
                 matfile='array_response_start.mat',new=False,src=True,
                 fout='array_response_start.pdf')
