#!/usr/bin/env mypython
"""
find the days with the greates number of stations 
"""
import sys
import os
sys.path.append(os.path.join(os.environ['AUTO_SRC'],'src/modules'))
import sac_db
from obspy.sac import *
import glob

if 1:
    def get_day_count(sdb,daydict):
        for _ne in xrange(sdb.nev):
            cnt = 0
            for _ns in xrange(sdb.nst):
                fn = os.path.dirname(sdb.rec[_ne][_ns].fname)+\
                     '/ft_grid_'+os.path.basename(sdb.rec[_ne][_ns].fname)
                if os.path.isfile(fn):
                    cnt +=1
            daydict[sdb.ev[_ne].name]=cnt


    sdbf = '/data/wanakaII/yannik/start/tmp/sac_db.out'
    sdb = sac_db.read_db(sdbf)
    daydict = {}
    get_day_count(sdb,daydict)
    cnt = 0
    for _k in daydict.keys():
        if daydict[_k] >= 20:
            print _k
            cnt += 1
    print "number of files: ",cnt

if 0:
    dirn = '/home/data/dev/beamforming/laura/all_stations/START_DATA_82/'
    os.chdir(dirn)
    fl = glob.glob('*')
    for _f in fl:
        x = ReadSac(_f,headonly=True)
        stat = x.kstnm.rstrip()
        os.mkdir(stat)
        os.rename(_f,os.path.join(stat,_f))
                                                    
