#!/usr/bin/env mypython
"""
Plot output from beamforming.
"""
from pylab import *
import scipy.io as sio
import sys

def polar_plot(beam,theta,slowness,dt,nfft,wtype,fout):
    df = dt/nfft
    periods = [6.]
    idx = [int(1./(p*df)) for p in periods]
    for ind in idx:
        tre = squeeze(beam[:,:,:,ind])
        tre = tre.mean(axis=2)
        tre = tre-tre.max()
        fig = figure(figsize=(6,6))
        ax = fig.add_subplot(1,1,1,projection='polar')
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
    savefig(fout)


if __name__ == '__main__':
    matfile = sys.argv[1]
    beam = sio.loadmat(matfile)['beam']
    dt = 1.0
    theta= arange(0,365,5)
    nfft = 128
    slowness = arange(0.125,0.51,0.01)
    fout = matfile.replace('.mat','.png')
    polar_plot(beam,theta,slowness,dt,nfft,'rayleigh',fout)
