#!/usr/local/python2/bin/python
####/usr/bin/env mypython
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
from ctypes import *
rcParams = {'backend':'Agg'}

DEBUG = True
def prep_beam(files,matfile,nhours=1,fmax=10.,threshold_std=0.5,onebit=True,
              tempfilter=False,fact=10,new=True):
    if new:
        ntimes = int(round(24/nhours))
        step = nhours*3600*fmax/fact

        stations = []
        slons = array([])
        slats = array([])
        nfiles = len(files)
        seisband = zeros((nfiles,ntimes,step))
        sigmas = []
        for i,_f in enumerate(files):
            tr = read(_f)[0]
            staname = tr.stats.station
            kcomp = tr.stats.channel
            slat = tr.stats.sac.stla
            slon = tr.stats.sac.stlo
            slons = append(slons,tr.stats.sac.stlo)
            slats = append(slats,tr.stats.sac.stla)
            if staname not in stations:
                stations.append(staname)
            tr.downsample(decimation_factor=fact, strict_length=True)
            npts = tr.stats.npts
            df = tr.stats.sampling_rate
            dt = tr.stats.delta
            seis0 = zeros(24*3600*int(df))
            seis0[0:npts] = tr.data
            seis0 -= seis0.mean()
            for j in xrange(ntimes):
                ilow = j*step
                iup = (j+1)*step
                seisband[i,j,:] = seis0[ilow:iup]
                seisband[i,j,:] -= seisband[i,j,:].mean()
                sigmas.append(seisband[i,j,:].std())
                if onebit:
                    seisband[i,j,:] = sign(seisband[i,j,:])
        if tempfilter:
            sgm = ma.masked_equal(array(sigmas),0.).compressed()
            sigma = sqrt(sum(sgm**2)/sgm.size)
            threshold = threshold_std*sigma
            seisband = where(abs(seisband) > threshold,threshold*sign(seisband),seisband)
            seisband = apply_along_axis(lambda e: e-e.mean(),2,seisband)

        if 1:
            fftpower = 7
            ismall = 2**fftpower
            ipick = arange(ismall)
            n=nhours*3600*df
            nsub = int(np.floor(n/ismall)) # Number of time pieces -20 mins long each
            seissmall = zeros((nfiles,ntimes,nsub,len(ipick)))
            for ii in xrange(nfiles):
                for jj in xrange(ntimes):
                    for kk in xrange(nsub):
                        seissmall[ii,jj,kk,:] = seisband[ii,jj,kk*ismall+ipick]

        fseis = fft(seissmall,n=2**fftpower,axis=3)
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

        sio.savemat(matfile,{'fseis':fseis,'seissmall':seissmall,'slats':slats,'slons':slons,'dt':dt})
        return fseis, meanlat, meanlon, slats, slons, dt, seissmall
    else:
        a = sio.loadmat(matfile)
        fseis = a['fseis']
        seissmall = a['seissmall']
        slats = a['slats']
        slons = a['slons']
        dt = a['dt'][0][0]
        slats = slats.reshape(slats.shape[0],)
        slons = slons.reshape(slons.shape[0],)
        meanlat = slats.mean()
        meanlon = slons.mean()
        return fseis, meanlat, meanlon, slats, slons, dt, seissmall
        

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



def beamforming(seis1,slowness,zetax,nsources,dt,new=True,matfile=None,freq_int=(0.02,0.4)):
    if new:
        _p = 6.
        _f = 1./_p
        freq = fftfreq(seis1.shape[3],dt)
        ind = searchsorted(freq[0:int(seis1.shape[3]/2)],_f,side='left')
        I = np.where((freq>freq_int[0]) & (freq<freq_int[1]))
        beam = zeros((nsources,slowness.size,seis1.shape[1]))
        ind = 21
        print ind, freq[ind]
        for ww in [ind]:
            FF = freq[ww]
            for cc in xrange(slowness.shape[1]):
                omega = 2*pi*FF
                velocity = 1./slowness[0][cc]*1000
                e_steer=exp(-1j*zetax*omega/velocity).T
                beamtemp = empty((len(theta),1))
                beamtemp = None
                for tt in xrange(seis1.shape[1]):
                    for TT in xrange(seis1.shape[2]):
                        Y = asmatrix(squeeze(seis1[:,tt,TT,ww],))
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


def beamforming_c(seis1,seissmall,slowness,zetax,nsources,dt,new=True,
                  matfile=None,freq_int=(0.02,0.4)):
    #int beam(double *traces, int stride1, int stride2, int stride3, int stride4, 
    #	 int nfft, int digfreq, double flow, double fhigh){
    lib = cdll.LoadLibrary('./beam_c.so')
    st1, st2, st3, st4 = seissmall.strides
    digfreq = int(round(1./dt))
    lib.beam.argtypes = [ \
        POINTER(c_double),
        c_int,
        c_int,
        c_int,
        c_int,
        c_int,
        c_int,
        c_int,
        c_int,
        c_int,
        c_double,
        c_double
        ]
    nstat,ntimes,nsub,nfft = seissmall.shape
    print zetax.shape
    ffttest = zeros(nfft)
    errcode = lib.beam(seissmall.ctypes.data_as(POINTER(c_double)),st1,st2,
                       st3,st4,nfft,nstat,ntimes,nsub,digfreq,freq_int[0],freq_int[1])
    #power = []
    #for i in xrange(nfft/2):
    #    re = ffttest[2*i]
    #    im = ffttest[2*i+1]
    #    power.append(sqrt(re*re+im*im))
    #freq = fftfreq(2*len(power),1.)
    #plot(freq[0:len(power)],power)
    #plot(freq[0:len(power)],abs(fft(seissmall[0,0,0,:],n=nfft)[0:len(power)]))
    #show()
    
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
        datdir = '/Volumes/Wanaka_01/yannik/start/sacfiles/10Hz/2001/Feb/2001_2_22_0_0_0/'
        #datdir = '/Volumes/Wanaka_01/yannik/start/sacfiles/10Hz/2001/Mar/2001_3_3_0_0_0/'
        files = glob.glob(os.path.join(datdir,'ft_grid*.HHZ.SAC'))
        matfile = 'prep_beam_2001_2_22.mat'
        fseis, meanlat, meanlon, slats, slons, dt, seissmall = prep_beam(files,matfile,onebit=False,tempfilter=True,new=False)
        if DEBUG:
            print 'calculating steering vector'
        zetax,theta,slowness,sta_origin_x,sta_origin_y = calc_steer(slats,slons)
        if DEBUG:
            print 'beamforming'
        beamforming_c(fseis,seissmall,slowness,zetax,theta.size,dt,new=True,matfile='6s_average_beam.mat')
        #beam = beamforming(fseis,slowness,zetax,theta.size,dt,new=True,matfile='6s_average_beam.mat')
        #polar_plot(beam,theta,slowness,resp=False)
    if 0:
        Ntimes, freqs, slons, slats, seis,meanlat,meanlon = read_matfiles()
        zetax,theta,slowness,sta_origin_x,sta_origin_y = calc_steer(slats,slons)
        beam = beamforming(seis,freqs,slowness,zetax,theta.size,Ntimes,new=False,matfile='6s_average_beam.mat')
        polar_plot(beam,theta,slowness,resp=True)
    if 0:
        arr_resp(freqs,slowness,zetax,theta,sta_origin_x,sta_origin_y,
                 matfile='array_response_start.mat',new=False,src=True,
                 fout='array_response_start.pdf')
