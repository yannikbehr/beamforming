#!/usr/bin/env mypython
"""
Plot output from beamforming.
"""
from pylab import *
import scipy.io as sio
import sys
import glob
import os
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import Normalize
from matplotlib import rcParams

def mycmp(beama, beamb):
    tmp = os.path.basename(beama).split('_')
    year = int(tmp[-6])
    month = int(tmp[-5])
    day = int(tmp[-4])
    datea = UTCDateTime(year,month,day,0,0,0,0).getTimeStamp()
    tmp = os.path.basename(beamb).split('_')
    year = int(tmp[-6])
    month = int(tmp[-5])
    day = int(tmp[-4])
    dateb = UTCDateTime(year,month,day,0,0,0,0).getTimeStamp()
    if datea > dateb: return 1
    if datea < dateb: return -1
    if datea == dateb: return 0

def getmonth(beam):
    tmp = os.path.basename(beam).split('_')
    year = int(tmp[-6])
    month = int(tmp[-5])
    day = int(tmp[-4])
    return month

def polar_plot(beam,theta,slowness,dt,nfft,wtype,fout=None):
    df = dt/nfft
    periods = [8.]
    #periods = [4.,5.,6.,7.,8.,9.,10.]
    idx = [int(1./(p*df)) for p in periods]
    for ind in idx:
        tre = squeeze(beam[:,:,:,ind])
        tre = tre.mean(axis=2)
        #tre = tre-tre.max()
        #tre = log10(abs(tre))
        fig = figure(figsize=(6,6))
        cax = fig.add_axes([0.85, 0.2, 0.05, 0.5])
        ax  = fig.add_axes([0.10, 0.1, 0.70, 0.7], polar=True)
        #ax = fig.add_subplot(1,1,1,projection='polar')
        cmap = cm.get_cmap('jet')
        ax.contourf((theta[::-1]+90.)*pi/180.,slowness,tre.T,
                    100,cmap=cmap,antialiased=True,
                    linstyles='dotted')
        ax.contour((theta[::-1]+90.)*pi/180.,slowness,tre.T,
                   100,cmap=cmap)
        ax.set_thetagrids([0,45.,90.,135.,180.,225.,270.,315.],
                          labels=['90','45','0','315','270','225','180','135'])
        ax.set_rgrids([0.1,0.2,0.3,0.4,0.5],labels=['0.1','0.2','0.3','0.4','0.5'],color='r')
        ax.grid(True)
        ax.set_title("%s %ds period"%(wtype,1./(ind*df)))
        ax.set_rmax(0.5)
        ColorbarBase(cax, cmap=cmap,
                     norm=Normalize(vmin=tre.min(), vmax=tre.max()))
    if fout is not None:
        savefig(fout)

def polar_plot_panel(beam,theta,slowness,dt,nfft,wtype,fout=None):
    rcParams['figure.subplot.left'] = 0.05
    rcParams['figure.subplot.right'] = 0.95
    rcParams['figure.subplot.top'] = 0.95
    rcParams['figure.subplot.bottom'] = 0.05
    rcParams['figure.subplot.wspace'] = 0.29
    df = dt/nfft
    periods = [6.]
    periods = [5.,6.,7.,8.,9.,10.,12.,15.,18.]
    idx = [int(1./(p*df)) for p in periods]
    cnt = 1
    cmap = cm.get_cmap('jet')
    fig = figure(figsize=(10,10))
    for ind in idx:
        tre = squeeze(beam[:,:,:,ind])
        tre = tre.mean(axis=2)
        #tre = tre-tre.max()
        #tre = log10(abs(tre))
        ax = fig.add_subplot(3,3,cnt,projection='polar')
        cnt += 1
        inds = where(slowness > 0.2)
        ax.contourf((theta[::-1]+90.)*pi/180.,slowness[inds],tre[:,inds[0]].T,
                    100,cmap=cmap,antialiased=True,
                    linstyles='dotted')
        ax.contour((theta[::-1]+90.)*pi/180.,slowness[inds],tre[:,inds[0]].T,
                   100,cmap=cmap)
        ax.set_thetagrids([0,45.,90.,135.,180.,225.,270.,315.],
                          labels=['90','45','0','315','270','225','180','135'])
        ax.set_rgrids([0.1,0.2,0.3,0.4,0.5],labels=['0.1','0.2','0.3','0.4','0.5'],color='r')
        ax.grid(True)
        ax.set_title("%s %ds period"%(wtype,1./(ind*df)))
        ax.set_rmax(0.5)
    if fout is not None:
        savefig(fout)

def monthly_average(fl,months,comp,new=True):
    fl.sort(cmp=mycmp)
    if new:
        for month in months.keys():
            avbeam = None
            cnt = 0
            for _f in fl:
                if getmonth(_f) == months[month]:
                    print _f
                    if comp == 'transverse':
                        beam = sio.loadmat(_f)['beamt']
                    if comp == 'vertical':
                        beam = sio.loadmat(_f)['beam']
                    else:
                        print 'comp has to be transverse or vertical'
                        return
                    
                    if avbeam is None:
                        avbeam = beam
                    else:
                        avbeam += beam
                    cnt += 1
            avbeam /= cnt
            sio.savemat(os.path.join(dirn,'average_beam_%s_%s.mat'%(comp,month)),{'avbeam':avbeam})
            dt = 1.0
            theta= arange(0,365,5)
            nfft = 128
            slowness = arange(0.125,0.51,0.01)
            fout = os.path.join(dirn,'average_beam_%s_%s.png'%(comp,month))
            polar_plot_panel(avbeam,theta,slowness,dt,nfft,'%s %s'%(comp,month),fout=fout)
    else:
        for month in months.keys():
            matfile = os.path.join(dirn,'average_beam_%s_%s.mat'%(comp,month))
            print matfile
            avbeam = sio.loadmat(matfile)['avbeam']
            dt = 1.0
            theta= arange(0,365,5)
            nfft = 128
            slowness = arange(0.125,0.51,0.01)
            fout = os.path.join(dirn,'average_beam_%s_%s.png'%(comp,month))
            polar_plot_panel(avbeam,theta,slowness,dt,nfft,'%s %s'%(comp,month),fout=fout)

def average(fl,comp,new=True):
    if new:
        avbeam = None
        for _f in fl:
            print _f
            if comp == 'transverse':
                beam = sio.loadmat(_f)['beamt']
            if comp == 'vertical':
                beam = sio.loadmat(_f)['beam']
            else:
                print 'comp has to be transverse or vertical'
                return
            if avbeam is None:
                avbeam = beam
            else:
                avbeam += beam
        avbeam /= len(fl)
        sio.savemat(os.path.join(dirn,'average_beam_%s.mat'%(comp)),{'avbeam':avbeam})
    else:
        avbeam = sio.loadmat(os.path.join(dirn,'average_beam_%s.mat'%(comp)))['avbeam']
    dt = 1.0
    theta= arange(0,365,5)
    nfft = 128
    slowness = arange(0.125,0.51,0.01)
    #polar_plot(avbeam,theta,slowness,dt,nfft,'rayleigh')
    polar_plot_panel(avbeam,theta,slowness,dt,nfft,'%s '%(comp))

if __name__ == '__main__':
    taranaki = False
    start = True
    comp = 'transverse'
    if taranaki:
        dirn = '/Volumes/GeoPhysics_05/users-data/yannik78/taranaki/beamforming/beams'
        months = {'march':3,'april':4,'may':5,'june':6,'july':7,'august':8,'september':9}
        if comp == 'transverse':
            fl = glob.glob(os.path.join(dirn,'beam_h*.mat'))
        if comp == 'vertical':
            fl = glob.glob(os.path.join(dirn,'beam_200*.mat'))
    if start:
        dirn = '/Volumes/Wanaka_01/yannik/start/beamforming'
        months = {'march':3,'april':4,'may':5,'june':6,'july':7,'august':8,'september':9}
        if comp == 'transverse':
            fl = glob.glob(os.path.join(dirn,'beam_h*.mat'))
        if comp == 'vertical':
            fl = glob.glob(os.path.join(dirn,'beam_200*.mat'))
    if 0:
        monthly_average(fl,months,comp,new=False)
    if 0:
        average(fl,comp,new=False)

### Bash lines to print the output of monthly_average
# for i in average*transverse*png;do convert $i -resize 75% `echo $i|cut -d. -f1`_small.png;done
# for i in average*transverse*small.png; do lpr -PCO505 $i;done
