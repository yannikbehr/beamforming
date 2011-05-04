#!/usr/bin/env mypython
"""
Measure dispersion curves from beamformer results.
"""

import os
import scipy.io as sio
from pylab import *

dirn = '/Volumes/GeoPhysics_05/users-data/yannik78/taranaki/beamforming/beams'
beamf = os.path.join(dirn,'average_beam_vertical_august.mat')
beam = sio.loadmat(beamf)['avbeam']
periods = [5.,6.,7.,8.,9.,10.,12.,15.,18.]
slowness = arange(0.125,0.51,0.01)
theta= arange(0,365,5)
dt = 1.0
nfft = 128
df = dt/nfft
idx = [int(1./(p*df)) for p in periods]
maxc = []
cnt = 1
cmap = cm.get_cmap('jet')
if 0:
    fig = figure(figsize=(10,10))
for ind in idx[0:1]:
    tre = squeeze(beam[:,:,:,ind])
    tre = tre.mean(axis=2)
    if 0:
        ax = fig.add_subplot(3,3,cnt,projection='polar')
        cnt += 1
        inds = where(slowness > 0.2)
        ax.contourf((theta[::-1]+90.)*pi/180.,slowness[inds],tre[:,inds[0]].T,
                    100,cmap=cmap,antialiased=True,
                    linstyles='dotted')
        ax.contour((theta[::-1]+90.)*pi/180.,slowness[inds],tre[:,inds[0]].T,
                   100,cmap=cmap)
        id1, id2 = unravel_index(tre.argmax(),tre.shape)
        ax.plot((theta[::-1][id1]+90.)*pi/180.,slowness[id2],'ko')
        ax.set_thetagrids([0,45.,90.,135.,180.,225.,270.,315.],
                          labels=['90','45','0','315','270','225','180','135'])
        ax.set_rgrids([0.1,0.2,0.3,0.4,0.5],labels=['0.1','0.2','0.3','0.4','0.5'],color='r')
        ax.grid(True)
        ax.set_rmax(0.5)
    if 1:
        hist,edges = histogram(tre.max(axis=1),normed=True,bins=10)
        for i,val in enumerate(hist):
            newval = val
            ax.bar(left=f-df/2,bottom=edges[i],height=edges[i+1]-edges[i],
                   width=df,color=cmap(newval),linewidth=0)
