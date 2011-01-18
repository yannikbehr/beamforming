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
import ctypes as C
import pickle
#import ipdb
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
        #seisband = zeros((step,ntimes,nfiles))
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
                #seisband[:,j,i] = seis0[ilow:iup]
                #seisband[:,j,i] -= seisband[:,j,i].mean()
                #sigmas.append(seisband[:,j,i].std())
                seisband[i,j,:] = seis0[ilow:iup]
                seisband[i,j,:] -= seisband[i,j,:].mean()
                sigmas.append(seisband[i,j,:].std())
                if onebit:
                    #seisband[:,j,i] = sign(seisband[:,j,i])
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
            #seissmall = zeros((len(ipick),ntimes,nsub,nfiles))
            seissmall = zeros((nfiles,ntimes,nsub,len(ipick)))
            for ii in xrange(nfiles):
                for jj in xrange(ntimes):
                    for kk in xrange(nsub):
                        #seissmall[:,jj,kk,ii] = seisband[kk*ismall+ipick,jj,ii]
                        seissmall[ii,jj,kk,:] = seisband[ii,jj,kk*ismall+ipick]
            
        #temp = fft(seissmall,n=2**fftpower,axis=0)
        #fseis = temp.copy()
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



def beamforming(seis,slowness,zetax,nsources,dt,new=True,matfile=None,freq_int=(0.02,0.4)):
    if new:
        _p = 6.
        _f = 1./_p
        #nfft,ntimes,nsub,nstat = seis.shape
        nstat,ntimes,nsub,nfft = seis.shape
        freq = fftfreq(nfft,dt)
        I = np.where((freq>freq_int[0]) & (freq<freq_int[1]))
        beam = zeros((nsources,slowness.size,ntimes,nfft))
        df = dt/nfft
        ind = int(_f/df)
        print ind, freq[ind]
        for ww in [ind]:
            FF = freq[ww]
            for cc in xrange(slowness.shape[1]):
                omega = 2*pi*FF
                velocity = 1./slowness[0][cc]*1000
                e_steer=exp(-1j*zetax*omega/velocity).T
                beamtemp = empty((len(theta),1))
                beamtemp = None
                #for tt in xrange(ntimes):
                for tt in [0]:
                    #for TT in xrange(nsub):
                    for TT in [0]:
                        #Y = asmatrix(squeeze(seis[ww,tt,TT,:],))
                        Y = asmatrix(squeeze(seis[:,tt,TT,ww],))
                        R = dot(Y.T,conjugate(Y))
                        if beamtemp is None:
                            beamtemp = atleast_2d(sum(abs(asarray(dot(conjugate(e_steer.T),dot(R,e_steer))))**2,axis=1))
                        else:
                            beamtemp = vstack((beamtemp,atleast_2d(sum(abs(asarray(dot(conjugate(e_steer.T),dot(R,e_steer))))**2,axis=1))))

                    beam[:,cc,tt,ww] = transpose(beamtemp).mean(axis=1)
        sio.savemat(matfile,{'beam':beam})
    else:
        beam = sio.loadmat(matfile)['beam']
    return beam


def ndarray2ptr2D(ndarray):
    """
    Construct ** pointer for ctypes from numpy.ndarray
    """
    ptr = C.c_void_p
    dim1, dim2 = ndarray.shape
    voids = []
    return (ptr * dim1)(*[row.ctypes.data_as(ptr) for row in ndarray])

def ndarray2ptr3D(ndarray):
    """
    Construct *** pointer for ctypes from numpy.ndarray
    """
    ptr = C.c_void_p
    dim1, dim2, _dim3 = ndarray.shape
    voids = []
    for i in xrange(dim1):
        row = ndarray[i]
        p = (ptr * dim2)(*[col.ctypes.data_as(ptr) for col in row])
        voids.append(C.cast(p, C.c_void_p))
    return (ptr * dim1)(*voids)

def ndarray2ptr4D(ndarray):
    """
    Construct **** pointer for ctypes from numpy.ndarray
    """
    dim1,dim2,dim3,dim4 = ndarray.shape
    ptr = C.c_void_p
    voids1 = []
    for i in xrange(dim1):
        voids2 = []
        for j in xrange(dim2):
            row = ndarray[i][j]
            p = (ptr * dim3)(*[col.ctypes.data_as(ptr) for col in row])
            voids2.append(C.cast(p, C.c_void_p))
        p2 = (ptr*dim2)(*voids2)
        voids1.append(C.cast(p2,C.c_void_p))
    return (ptr*dim1)(*voids1)

def beamforming_c(seis1,seissmall,slowness,zetax,nsources,dt,theta,new=True,
                  matfile=None,freq_int=[0.02,0.4]):
    lib = C.cdll.LoadLibrary('./beam_c.so')
    digfreq = int(round(1./dt))
    lib.beam.argtypes = [ \
        C.c_void_p,
        C.c_void_p,
        C.c_void_p,
        np.ctypeslib.ndpointer(dtype='float64', ndim=1, flags='C_CONTIGUOUS'),
        C.c_int,
        C.c_int,
        C.c_int,
        C.c_int,
        C.c_int,
        C.c_int,
        C.c_int,
        C.c_double,
        C.c_double
        ]
    #nfft,ntimes,nsub,nstat = seissmall.shape
    nstat,ntimes,nsub,nfft = seissmall.shape
    df = digfreq/float(nfft)
    wlow = int(freq_int[0]/df+0.5)
    whigh = int(freq_int[1]/df+0.5)
    beam = zeros((nfft,nsources,slowness.size))
    beam_p = ndarray2ptr3D(beam)
    zetax_p = ndarray2ptr2D(zetax)
    nslow = slowness.size
    if 0:
        seissmall_p = ndarray2ptr4D(seissmall)
        errcode = lib.beam_fft(C.byref(seissmall_p),C.byref(beam_p),C.byref(zetax_p),
                           slowness[0],nslow,nsources,nfft,nstat,ntimes,nsub,digfreq,
                           freq_int[0],freq_int[1])
        print 'python',seissmall[3,4,6,100]
        f = open('beam_c.dump','w')
        pickle.dump(beam,f)
        f.close()
    else:
        f = open('beam_c.dump')
        beam = pickle.load(f)
        
    if new:
        beam = zeros((nsources,slowness.size,ntimes,nfft))
        beam_p = ndarray2ptr4D(beam)
        ntimes = 1
        nsub = 1
        fnc = lib.beam
        fnc.argtypes = [\
            np.ctypeslib.ndpointer(dtype='complex128',ndim=4,flags='aligned, contiguous'),
            C.POINTER(C.c_long),
            C.POINTER(C.c_long),
            C.c_void_p,
            C.c_void_p,
            np.ctypeslib.ndpointer(dtype='float64', ndim=1, flags='C_CONTIGUOUS'),
            C.c_int,
            C.c_int,
            C.c_int,
            C.c_int,
            C.c_int,
            C.c_int,
            C.c_int,
            np.ctypeslib.ndpointer(dtype='float64', ndim=1, flags='C_CONTIGUOUS')
            ]
        errcode = fnc(seis1,seis1.ctypes.strides,seis1.ctypes.shape,
                      C.byref(beam_p),C.byref(zetax_p),slowness[0],
                      nslow,nsources,nfft,nstat,ntimes,nsub,digfreq,array(freq_int))

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
        tre = squeeze(beam[:,:,:,21])
        tre = tre.mean(axis=2)
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
    if 1:
        if DEBUG:
            print 'preparing raw data'
        datdir = '/Volumes/Wanaka_01/yannik/start/sacfiles/10Hz/2001/Feb/2001_2_22_0_0_0/'
        #datdir = '/Volumes/Wanaka_01/yannik/start/sacfiles/10Hz/2001/Mar/2001_3_3_0_0_0/'
        files = glob.glob(os.path.join(datdir,'ft_grid*.HHZ.SAC'))
        matfile = 'prep_beam_2001_2_22.mat'
        fseis, meanlat, meanlon, slats, slons, dt, seissmall = prep_beam(files,matfile,onebit=False,tempfilter=True,new=True)
        if DEBUG:
            print 'calculating steering vector'
        zetax,theta,slowness,sta_origin_x,sta_origin_y = calc_steer(slats,slons)
        if DEBUG:
            print 'beamforming'
        beam = beamforming_c(fseis,seissmall,slowness,zetax,theta.size,dt,theta,new=True,matfile='6s_short_beam_c.mat')
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
