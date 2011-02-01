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

def polar_plot(beam,theta,slowness,dt,nfft,wtype,fout=None):
    df = dt/nfft
    periods = [6.]
    idx = [int(1./(p*df)) for p in periods]
    for ind in idx:
        tre = squeeze(beam[:,:,:,ind])
        tre = tre.mean(axis=2)
        tre = tre-tre.max()
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


if __name__ == '__main__':
    dirn = '/data/wanakaII/yannik/taranaki/beamforming/beams/'
    fl = glob.glob(os.path.join(dirn,'beam_2002_*.mat'))
    if 1:
        avbeam = None
        for _f in fl:
            print _f
            beam = sio.loadmat(_f)['beam']
            if avbeam is None:
                avbeam = beam
            else:
                avbeam += beam
        avbeam /= len(fl)
        sio.savemat(os.path.join(dirn,'average_beam_2002_4_vertical.mat'),{'avbeam':avbeam})
    else:
        avbeam = sio.loadmat(os.path.join(dirn,'average_beam_2002_4_vertical.mat'))['avbeam']
    dt = 1.0
    theta= arange(0,365,5)
    nfft = 128
    slowness = arange(0.125,0.51,0.01)
    polar_plot(avbeam,theta,slowness,dt,nfft,'rayleigh')
    show()
