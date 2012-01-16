import os

def PhonyTarget(target, action):
    phony = Environment(ENV = os.environ,BUILDERS = { 'phony' : Builder(action = action) })
    AlwaysBuild(phony.phony(target = target, source = 'SConstruct'))

def SemiPhony(target,source,action):
    env = Environment(ENV = os.environ,BUILDERS = { 'cmd' : Builder(action = action) })
    env.cmd(target = target, source = source)

# compute average as done by the matlab script plot_freq_azimuth_day.m written by L. Brooks
dirn = '/Volumes/GeoPhysics_05/users-data/yannik78/taranaki/beamforming/beams/'
fout = os.path.join(dirn,'average_beam_vertical_6s_247kms.mat')
PhonyTarget('avbtara','ls %s |./plot_average_beam.py -l -o %s'%(dirn+'beam_200*.mat',fout))

# compute simple mean over all beamformer outputs
dirn = '/Volumes/GeoPhysics_05/users-data/yannik78/taranaki/beamforming/beams/'
fout = os.path.join(dirn,'average_beam_vertical.mat')
PhonyTarget('avbtaraall','ls %s |./plot_average_beam.py -o %s'%(dirn+'beam_200*.mat',fout))

# compute simple mean over beamformer outputs between April and July
dirn = '/Volumes/GeoPhysics_05/users-data/yannik78/taranaki/beamforming/beams/'
fout = os.path.join(dirn,'average_beam_vertical_4months.mat')
fin = 'vertical_beams_4months_taranaki.txt'
PhonyTarget('avbtara4m','./plot_average_beam.py -o %s < %s'%(fout,fin))

# compute simple mean over beamformer outputs between April and July at 6s period
# and 2.47 km/s only
dirn = '/Volumes/GeoPhysics_05/users-data/yannik78/taranaki/beamforming/beams/'
fout = os.path.join(dirn,'average_beam_vertical_4months_6s_247kms.mat')
fin = 'vertical_beams_4months_taranaki.txt'
PhonyTarget('avbtara4m6s','./plot_average_beam.py -l -o %s < %s'%(fout,fin))

# compute simple mean over all horizontal beamformer outputs
dirn = '/Volumes/GeoPhysics_05/users-data/yannik78/taranaki/beamforming/beams/'
fout = os.path.join(dirn,'average_beam_horizontal.mat')
PhonyTarget('avbtaraallh','ls %s |./plot_average_beam.py --comp=transverse -o %s'%(dirn+'beam_h_200*.mat',fout))

# compute simple mean over all horizontal beamformer outputs from the second beamformer run
dirn = '/Volumes/GeoPhysics_05/users-data/yannik78/taranaki/beamforming/beams_2/'
fout = os.path.join(dirn,'average_beam_horizontal_2nd_run.mat')
PhonyTarget('avbtaraallh2','ls %s |./plot_average_beam.py --comp=transverse -o %s'%(dirn+'beam_h_200*.mat',fout))


#dirn = '/Volumes/Wanaka_01/yannik/start/beamforming'
#            fl = glob.glob(os.path.join(dirn,'beam_h*.mat'))
#            fl = glob.glob(os.path.join(dirn,'beam_200*.mat'))
