#!/usr/bin/env mypython
"""
Run FTAN on synthetic surface wave trains.
"""

import os, sys, string, glob
sys.path.append(os.environ['AUTO_SRC']+'/src/FTAN/gv')
sys.path.append(os.environ['AUTO_SRC']+'/src/FTAN/pv')
sys.path.append(os.path.join(os.environ['PROC_SRC'],'NA'))
sys.path.append(os.path.join(os.environ['PROC_SRC'],'misc'))
from dinver_run import get_disp
from herman_surface_wave import herman_syn
import ftanc
import ftangv
from obspy.sac import *
import scipy.interpolate as scint

model = ['4\n','3.  6.0 3.5 2.5 100. 200.\n',
         '2.  3.4 2.0 2.3 100. 200.\n',
         '5.  6.5 3.8 2.5 100. 200.\n',
         '0.  8.0 4.7 3.3 500. 900.\n']
rayc,lovc = get_disp(model,0.02,1.0,nmode=1)
rayu,lovu = get_disp(model,0.02,1.0,gv=True)
indlc = where(lovc[:,0] > 0.)
indrc = where(rayc[:,0] > 0.)
indlu = where(lovu[:,0] > 0.)
indru = where(rayu[:,0] > 0.)
if 0:
    plot(1./lovc[indlc,0][0],1./lovc[indlc,1][0],label='Love phase')
    plot(1./rayc[indrc,0][0],1./rayc[indrc,1][0],label='Rayleigh phase')
    plot(1./lovu[indlu,0][0],1./lovu[indlu,1][0],label='Love group')
    plot(1./rayu[indru,0][0],1./rayu[indru,1][0],label='Rayleigh group')
    xlabel('Period [s]')
    ylabel('Velocity [km/s]')
    legend(loc='lower right')

x = 499.89334
df = 0.01
fint = arange(0.02,1.01,0.0001)
wtype = 'rayleigh'
if wtype == 'rayleigh':
    frc = rayc[indrc,0][0]
    fru = rayu[indru,0][0]
    c = 1./rayc[indrc,1][0]
    u = 1./rayu[indru,1][0]
elif wtype == 'love':
    frc = lovc[indlc,0][0]
    fru = lovu[indlu,0][0]
    c = 1./lovc[indlc,1][0]
    u = 1./lovu[indlu,1][0]
else:
    print "incorrect wave type [rayleigh or love]"
dom = 2*pi*df
om0 = 0.005+fint[:-1]*2*pi
om0 = fint*2*pi
t = linspace(0,208,208)
t = linspace(0,256,256)
repc = scint.splrep(2*pi*frc,c)
repu = scint.splrep(2*pi*fru,u)
fsum = 0
 
for _w in om0:
    y = dom/2.*(t-(x/scint.splev(_w,repu)))
    f0 = dom/pi*sin(y)/y*cos(_w*t-_w*x/scint.splev(_w,repc))
    fsum += f0

if 1:
    figure()
    subplot(2,1,1)
    plot(t,fsum/fsum.max(),label='Rayleigh (Aki&Richards)')
    model = ['3. 6.0 3.5 2.5 100.0 200.0 0.0 0.0 1.0 1.0\n',
             '2. 3.4 2.0 2.3 100. 200. 0.0 0.0 1.0 1.0\n',
             '5. 6.5 3.8 2.5 100. 200. 0.0 0.0 1.0 1.0\n',
             '0. 8.0 4.7 3.3 500.0 900.0 0.0 0.0 1.0 1.0\n']
    z,r,t = herman_syn(model,x,npts=256)
    plot(z.data/r.data.max(),label='Rayleigh (Hermann)')
    legend(loc='upper left')
    subplot(2,1,2)
    plot(z.data,label='Rayleigh (Hermann)')
    legend(loc='upper left')

rtr = SacIO()
rtr.fromarray(fsum,distkm=x)
rtr.WriteSacBinary('rayleigh_synthetic.sac')



if 1:
    rfout = 'syn_rayleigh_phase.txt'
    savetxt(rfout,vstack((1./rayc[indrc,0][0][::-1],1./rayc[indrc,1][0][::-1])).T)

if 0:
    rfn = 'rayleigh_synthetic.sac'
    rfn = 'herman_syn.sac'
    t.write(rfn,format='SAC')
    tr = SacIO(rfn)
    cper,caper,gv,pv,gvamp,gvsnr,gvwdth,ampv,amps,refper,refvel = ftanc.myftan(tr,rfout,
                                                                              tmin=2.0,tmax=50.0,
                                                                              vmax=5.0,vmin=1.5,
                                                                              ffact=1.0,
                                                                              extrace=None,
                                                                              level='easy',
                                                                               phm=False)
    cper,aper,gv,gvamp,gvsnr,ampv,amps = ftangv.myftan(tr,tmin=3.0,vmax=5.0,vmin=1.0,ffact=1.0,
                                                       level='easy',phm=False,tmaxmax=50.)
    figure()
    contourf(cper,ampv,amps,100)
    plot(1./lovu[indlu,0][0][::-1],1./lovu[indlu,1][0][::-1],label='Theoretical Love wave group velocity dispersion curve',c='k')
    plot(1./lovc[indlc,0][0][::-1],1./lovc[indlc,1][0][::-1],label='Theoretical Love wave phase velocity dispersion curve',c='g')
    plot(aper,gv,label='measured group velocity',c='k',ls='--')
    plot(caper,pv,label='measured phase velocity',c='g',ls='--')
    ax = gca()
    ax.autoscale_view(tight=True)
    legend(loc='lower right')
    xlim(3.0,50.0)


if 1:
    rfn = 'rayleigh_synthetic.sac'
    rfn = 'herman_syn.sac'
    z.write(rfn,format='SAC')
    tr = SacIO(rfn)
    cper,caper,gv,pv,gvamp,gvsnr,gvwdth,ampv,amps,refper,refvel = ftanc.myftan(tr,rfout,
                                                                              tmin=2.0,tmax=50.0,
                                                                              vmax=5.0,vmin=1.5,
                                                                              ffact=1.0,
                                                                              extrace=None,
                                                                              level='easy',
                                                                               phm=False)
    cper,aper,gv,gvamp,gvsnr,ampv,amps = ftangv.myftan(tr,tmin=3.0,vmax=5.0,vmin=1.0,ffact=1.0,
                                                       level='easy',phm=False,tmaxmax=50.)
    figure()
    contourf(cper,ampv,amps,100)
    plot(1./rayu[indru,0][0][::-1],1./rayu[indru,1][0][::-1],label='Theoretical Rayleigh wave group velocity dispersion curve',c='k')
    plot(1./rayc[indrc,0][0][::-1],1./rayc[indrc,1][0][::-1],label='Theoretical Rayleigh wave phase velocity dispersion curve',c='g')
    plot(aper,gv,label='measured group velocity',c='k',ls='--')
    plot(caper,pv,label='measured phase velocity',c='g',ls='--')
    ax = gca()
    ax.autoscale_view(tight=True)
    legend(loc='lower right')
    xlim(3.0,50.0)
