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
import scipy.io as sio

def get_sac_list():
    sacdir = '/data/wanakaII/yannik/cnipse/sacfiles/2001/Apr/2001_4_30_0_0_0/'
    sacdirII = '/data/wanakaII/yannik/start/sacfiles/2001/Apr/2001_4_30_0_0_0/'
    a = sio.loadmat('/home/data/dev/proc-scripts_git/beamforming/BeamformInputData_start_cnipse.mat',struct_as_record=True)
    fl = glob.glob(os.path.join(a['matpath'][0],'*.mat'))
    newlist = []
    for _f in fl:
        stat = os.path.basename(_f).split('_')[0]
        sacf = glob.glob(os.path.join(sacdir,'*'+stat+'*HZ.SAC'))
        if len(sacf) < 1:
            sacf = glob.glob(os.path.join(sacdirII,'*'+stat+'*HZ.SAC'))
        newlist.append(sacf[0])
    return newlist

def arr_geom(filelist):
    """
    calculate array geometry from given
    list of sac-files
    """
    slon = []
    slat = []
    sta_origin_dist = array([])
    sta_origin_bearing = array([])
    for _f in filelist:
        x = ReadSac(_f,headonly=True)
        slon.append(x.GetHvalue('stlo'))
        slat.append(x.GetHvalue('stla'))

    meanlat = mean(slat)
    meanlon = mean(slon)
    for lat,lon in zip(slat,slon):
        dist, az, baz = rotate.gps2DistAzimuth(meanlat,meanlon,lat,lon)
        sta_origin_dist = append(sta_origin_dist,dist)
        sta_origin_bearing = append(sta_origin_bearing,az)

    sta_origin_x = sta_origin_dist*sin(sta_origin_bearing*pi/180.)
    sta_origin_y = sta_origin_dist*cos(sta_origin_bearing*pi/180.)
    return sta_origin_x, sta_origin_y


def plot_arr(x,y):
    """
    plot array geometry
    """
    fig = figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot(x,y,'b*')
    ax.set_xlabel('X-dist (km)')
    ax.set_ylabel('Y-dist (km)')
    return fig


def get_delay(theta,x,y):
    """
    calculate projection of location vector onto
    steering vector
    """
    zeta_x = -cos(theta*pi/180.)
    zeta_y = -sin(theta*pi/180.)
    zeta_x.resize(len(zeta_x),1)
    zeta_y.resize(len(zeta_y),1)
    x.resize(1,len(x))
    y.resize(1,len(y))
    zeta = dot(zeta_x,x) + dot(zeta_y,y)
    return zeta


def arr_resp_src(zeta,slowness,freq,R,x,y):
    """
    calculate response for source at theta1 with slowness s1
    """
    theta1 = 270
    zeta_x = -cos(theta1*pi/180.)
    zeta_y = -sin(theta1*pi/180.)
    zeta_src = zeta_x*x + zeta_y*y
    s1 = 0.3
    
    beam = zeros((len(slowness),len(freq),len(theta)))
    for _s in slowness:
        for _f in freq:
            omega = 2*pi*_f
            velocity = 1./_s*1000
            e_steer = exp(-1j*zeta*omega/velocity)
            c1 = 1./s1*1000
            e_src = exp(1j*zeta_src*omega/c1)
            e_steer *= e_src
            beam[where(slowness==_s),where(freq==_f),:] = diag(abs(dot(dot(e_steer,R),conj(e_steer).T))**2)
    return beam

def arr_resp(zeta,slowness,freq,R):
    """
    calculate array response
    """
    beam = zeros((len(slowness),len(freq),len(theta)))
    for _s in slowness:
        for _f in freq:
            omega = 2*pi*_f
            velocity = 1./_s*1000
            e_steer = exp(-1j*zeta*omega/velocity)
            beam[where(slowness==_s),where(freq==_f),:] = diag(abs(dot(dot(e_steer,R),conj(e_steer).T))**2)
    return beam

def plot_beam(ax,theta,slowness,beam):
    """
    Make a polar plot for the beamformer output.
    """
    project='polar'
    #project='rectilinear'
    #ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], projection=project, axisbg='white')
    ax.contourf(theta*pi/180.,slowness,beam,100)
    ax.set_rmax(0.5)
    ax.grid(True)
    #y = mat(slowness)
    #x = mat(theta*pi/180.).T
    #print x.shape, y.shape
    #Y1 = dot(ones((len(x),1)),y)
    #X1 = dot(x,ones((1,len(y))))
    return fig

def load_trace():
    """
    Load .mat-files with preprocessed traces
    """
    a = sio.loadmat('/home/data/dev/proc-scripts_git/beamforming/BeamformInputData.mat')
    fl = glob.glob(os.path.join(a['matpath'][0],'*.mat'))
    nst = len(fl)
    Nsub = int(a['Nsub'])
    Ntimes = int(a['Ntimes'])
    Nfreq = len(a['I'][0])
    seis = zeros((Nfreq,Nsub,Ntimes,nst),dtype='complex128')
    _nst = 0
    for _f in fl:
        b = sio.loadmat(_f)
        fseis = b['fseis']
        seis[:,:,:,_nst] = fseis.copy()
        _nst += 1
    return  seis

def get_freqs():
    """
    Load .mat-file and get frequency array
    """
    a = sio.loadmat('/home/data/dev/proc-scripts_git/beamforming/BeamformInputData.mat')
    return a['freq'][a['I']]
    

def run_beam(zeta,slowness,freqs,freq_ind,traces):
    """
    run beamforming
    """
    Nsub = traces.shape[1]
    Ntimes = traces.shape[2]
    beam = zeros((len(slowness),len(freqs[0]),Ntimes,len(theta)))
    for _f in freqs[0][freq_ind]:
        ff = where(freqs[0]==_f)
        print "beamforming for frequency:",_f
        omega = 2*pi*_f
        for _s in slowness:
            velocity = 1./_s*1000
            e_steer = exp(-1j*zeta*omega/velocity)
            for tt in range(Ntimes):
                beamtemp = empty((len(theta),1))
                for TT in range(Nsub):
                    Y = asmatrix(squeeze(traces[ff,TT,tt,:]))
                    R = dot(Y.T,Y)
                    beamtemp = append(beamtemp,asmatrix(diag(abs(dot(dot(e_steer,R),e_steer.T))**2)).T,axis=1)
                    #beamtemp = append(beamtemp,sum(abs(dot(dot(e_steer,R),e_steer.T))**2,axis=1),axis=1)
                beam[where(slowness==_s),ff,tt,:] = squeeze(beamtemp[:,1::].mean(axis=1))
    return beam

if __name__ == '__main__':
    sta_origin_x, sta_origin_y = arr_geom(get_sac_list())
    fig = plot_arr(sta_origin_x, sta_origin_y)
    if True:
        theta= arange(0,362,2)
        zeta = get_delay(theta, sta_origin_x, sta_origin_y)
        slowness=arange(0.03,0.505,0.005)  ###slowness in s/km
        R=ones((46,46))
        indx = array([0,10,30,48])
        #beam = run_beam(zeta,slowness,get_freqs(),indx,load_trace())
        #sio.savemat('beam.mat',{'beam':beam})
        #beam = sio.loadmat('beam.mat')['beam']
        cnt = 1
        fig = figure(figsize=(6, 6))
        beam = arr_resp(zeta,slowness,get_freqs()[0],R)
        #beam = arr_resp_src(zeta,slowness,get_freqs()[0],R,sta_origin_x, sta_origin_y)
        for _f in get_freqs()[0][indx]:
            ff = where(get_freqs()[0]==_f)
            #tre = squeeze(beam.mean(axis=2)[:,ff,:])
            tre = squeeze(beam[:,ff,:])
            #tre = 10*log10(tre)
            ax = fig.add_subplot(2,2,cnt, projection='polar')
            plot_beam(ax,theta,slowness,tre)
            cnt += 1
            title(str(_f))
    show()


