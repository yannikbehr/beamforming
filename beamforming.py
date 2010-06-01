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
from matplotlib import cm

def get_sac_list(sacdirs,matfile):
    a = sio.loadmat(matfile,struct_as_record=True)
    fl = glob.glob(os.path.join(a['matpath'][0],'*.mat'))
    newlist = []
    for _f in fl:
        #stat = os.path.basename(_f).split('_')[0]
        stat = os.path.basename(_f).split('62')[0]
        for _sd in sacdirs:
            sacf = glob.glob(os.path.join(_sd,'*'+stat+'*HZ.SAC'))
            if len(sacf) > 1:
                newlist.append(sacf[0])
                break
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

    sta_origin_x = sta_origin_dist*cos(sta_origin_bearing*pi/180.)
    sta_origin_y = sta_origin_dist*sin(sta_origin_bearing*pi/180.)
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
    try:
        beam = zeros((len(slowness),len(freq),len(theta)))
        for _s in slowness:
            for _f in freq:
                omega = 2*pi*_f
                velocity = 1./_s*1000
                e_steer = exp(-1j*zeta*omega/velocity)
                beam[where(slowness==_s),where(freq==_f),:] = diag(abs(dot(dot(e_steer,R),conj(e_steer).T))**2)
    except:
        import pdb; pdb.set_trace()
    return beam

def plot_beam(ax,theta,slowness,beam):
    """
    Make a polar plot for the beamformer output.
    """
    cmap = cm.get_cmap('jet')
    ax.contourf((theta[::-1]-90.)*pi/180.,slowness,beam,100,cmap=cmap,antialiased=True,
                linstyles='dotted')
    ax.set_rmax(0.5)
    ax.set_thetagrids([0,45.,90.,135.,180.,225.,270.,315.],
                      labels=['90','45','0','315','270','225','180','135'])
    ax.set_rgrids([0.1,0.2,0.3,0.4,0.5],labels=['0.1','0.2','0.3','0.4','0.5'],color='r')
    ax.grid(True)
    return 

def load_trace(matfile):
    """
    Load .mat-files with preprocessed traces
    """
    a = sio.loadmat(matfile)
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

def get_freqs(matfile):
    """
    Load .mat-file and get frequency array
    """
    a = sio.loadmat(matfile)
    return a['freq'][a['I']]
    

def run_beam(zeta,slowness,freqs,freq_ind,traces,new=False,fout='beam.mat'):
    """
    run beamforming
    """
    if new:
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
                        #beamtemp = append(beamtemp,asmatrix(diag(abs(dot(dot(e_steer,R),e_steer.T))**2)).T,axis=1)
                        beamtemp = append(beamtemp,sum(abs(dot(dot(e_steer,R),e_steer.T))**2,axis=1),axis=1)
                    beam[where(slowness==_s),ff,tt,:] = squeeze(beamtemp[:,1::].mean(axis=1))
        sio.savemat(fout,{'beam':beam})
    if not new:
        beam = sio.loadmat(fout)['beam']
    return beam

if __name__ == '__main__':
    theta= arange(0,362,2)
    slowness=arange(0.03,0.505,0.005)  ###slowness in s/km
    indx = array([18])
    if True:
        matfile = '/home/data/dev/proc-scripts_git/beamforming/BeamformInputData_start.mat'
        sacdirs = ['/data/wanakaII/yannik/start/sacfiles/2001/Mar/2001_3_3_0_0_0/']
        sta_origin_x, sta_origin_y = arr_geom(get_sac_list(sacdirs,matfile))
        fig = plot_arr(sta_origin_x, sta_origin_y)
    if False:
        zeta = get_delay(theta, sta_origin_x, sta_origin_y)
        #beam = arr_resp(zeta,slowness,get_freqs()[0],R)
        #beam = arr_resp_src(zeta,slowness,get_freqs()[0],R,sta_origin_x, sta_origin_y)
        #R=ones((sta_origin_x.size,sta_origin_x.size))
        #indx = array([0,10,30,48])
        beam = run_beam(zeta,slowness,get_freqs(matfile),indx,load_trace(matfile),new=True,fout='beam_test.mat')
    if False:
        rawdatf = 'BeamformInputData53.mat'
        rawdat = sio.loadmat(rawdatf)
        
    if True:
        cnt = 1
        fig = figure(figsize=(6, 6))
        beam1 = sio.loadmat('beam62.mat')['beam']
        for _f in get_freqs(matfile)[0][indx]:
            ff = where(get_freqs(matfile)[0]==_f)
            tre = squeeze(beam[:,:,ff,:])
            tre = squeeze(beam1[:,ff,:,:])
            #tre = tre.mean(axis=1)
            tre = tre.mean(axis=2)
            tre = tre-tre.max()
            #tre = 10*log10(tre)
            print _f
            ax = fig.add_subplot(1,1,cnt, projection='polar')
            plot_beam(ax,theta,slowness,tre.T)
            cnt += 1
            title(str(_f))
        show()


