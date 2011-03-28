#!/usr/bin/env mypython
"""
Project beamformer output onto the coastline.
"""

from gmtpy import GMT
import os
from collections import defaultdict
from pylab import *
import scipy.io as sio
from obspy.signal.rotate import gps2DistAzimuth
import glob
import tempfile
import cStringIO

os.environ['GMTHOME'] = '/usr/local/gmt/'
gmt = GMT()
rng = '160./181/-48/-37.'
scl = 'M10c'
gmt.pscoast(R=rng,J=scl,B='a0.5wsne',D='l',A='400/1/1',W='thinnest',m=True)
a = gmt.output.getvalue().split('\n')
d = defaultdict(list)
cnt = -1
for _l in a:
    if _l.startswith('#'):continue
    if _l.startswith('>'):
        cnt += 1
        continue
    try:
        d[cnt].append(map(float,_l.split('\t')))
    except ValueError:
        print _l

ni_cl = vstack((array(d[0]),array(d[4])))
si_cl = vstack((array(d[1]),array(d[3])))
if 0:
    plot(ni_cl[:,0],ni_cl[:,1],'b.')
    plot(si_cl[:,0],si_cl[:,1],'g.')

mean_lat_tara = -39.34574
mean_lon_tara = 174.18534
dirn = '/Volumes/GeoPhysics_05/users-data/yannik78/taranaki/beamforming/beams'
beams = glob.glob(os.path.join(dirn,'beam_h_2002_*.mat'))
beams = glob.glob(os.path.join(dirn,'average_beam_transverse.mat'))
ntimes = 24
nfft = 128
dt = 1.0
df = dt/nfft
periods = [5.,6.,7.,8.,9.,10.,12.,15.,18.]
#periods = [8.]
indices = [int(1./(p*df)) for p in periods]
slowness = np.arange(0.125,0.51,0.01)
theta= np.arange(0,365,5)
nz_cl = vstack((ni_cl,si_cl))
nval = zeros((nz_cl.shape[0],len(indices),3))
comp = 't'
for beamf in beams:
        a = sio.loadmat(beamf)
        if comp == 'v':
            beam = a['beam']
        elif comp == 't':
            beam = a['avbeam']
        elif comp == 'r':
            beam = a['beamr']
        else:
            print "argument 'comp' as to be either 'v', 't' or 'r'"
        nsources, nvels, ntimes, nfft = beam.shape
        for k,freqind in enumerate(indices):
            for cnt,coord in enumerate(nz_cl):
                lon, lat = coord
                dist, az, baz = gps2DistAzimuth(lat, lon, mean_lat_tara, mean_lon_tara)
                ind = int(baz/5.)
                nb = np.squeeze(beam[:,:,:,freqind])
                nb = nb[:,:,0:22].mean(axis=2)
                bmax = sum(nb[ind+1,:])*0.01
                bmin = sum(nb[ind,:])*0.01
                grd = (bmax - bmin)/5.
                nval[cnt,k,0] = lon
                nval[cnt,k,1] = lat
                nval[cnt,k,2] = (bmin + grd*(baz-theta[ind]))*np.sqrt(dist/1000.)

ind = 6
fstr = cStringIO.StringIO()
for _i in range(nval.shape[0])[0:ni_cl.shape[0]-1]:
    lon0,lat0,val0 = nval[_i,ind,:]
    lon1,lat1,val1 = nval[_i+1,ind,:]
    fstr.write("> -Z%s\n%.4f\t%.4f\n%.4f\t%.4f\n"%((val0+val1)/2.,lon0,lat0,lon1,lat1))
for _i in range(nval.shape[0])[ni_cl.shape[0]+1:-1]:
    lon0,lat0,val0 = nval[_i,ind,:]
    lon1,lat1,val1 = nval[_i+1,ind,:]
    fstr.write("> -Z%s\n%.4f\t%.4f\n%.4f\t%.4f\n"%((val0+val1)/2.,lon0,lat0,lon1,lat1))

gmt = GMT()
cptfile = gmt.tempfilename()
rng = '160./179/-48/-32.'
scl = 'M10ic'
anot = 'a5f2'
scalebar = '3.c/-1.6c/6c/.2ch'
gmt.makecpt(C='seis',I=True,T='%f/%f/1.0'%(nval[:,ind,2].min(),nval[:,ind,2].max()),
            D=True,out_filename=cptfile)
gmt.psbasemap(R=rng,J=scl,B=anot,G='white')
gmt.psxy(R=True,J=True,m=True,W='2',C=cptfile,in_string=fstr.getvalue())
gmt.psscale(C=cptfile,D=scalebar,B='%f::/::'%5.)

fout = './.eps'
gmt.save(fout)


if 0:
    ind = 6
    scatter(nval[:,ind,0],nval[:,ind,1],c=nval[:,ind,2],
            norm=Normalize(vmin=nval[:,ind,2].min(),vmax=nval[:,ind,2].max()))
    plot(ni_cl[:,0],ni_cl[:,1],'k-')
    plot(si_cl[:,0],si_cl[:,1],'k-')
    colorbar()

if 0:
    fig = plt.figure(figsize=(6,6))
    cax = fig.add_axes([0.85, 0.2, 0.05, 0.5])
    ax  = fig.add_axes([0.10, 0.1, 0.70, 0.7], polar=True)
    tre = np.squeeze(beam[:,:,:,int(1./(12.0*df))])
    tre = tre.mean(axis=2)
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
    ax.set_rmax(0.6)

