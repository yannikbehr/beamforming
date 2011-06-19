#!/bin/sh
#
#$ -S /bin/sh
#### 1-203
#$ -t 1-122
#$ -q all.q
#$ -m n
#$ -N vertical_beamforming
### /Volumes/GeoPhysics_05/users-data/yannik78/taranaki/beamforming/beams
#$ -wd /Volumes/GeoPhysics_05/users-data/yannik78/taranaki/beamforming/laura_beam_test/beams
### /Volumes/GeoPhysics_05/users-data/yannik78/taranaki/beamforming/beams/vertical_beam.out
### /Volumes/GeoPhysics_05/users-data/yannik78/taranaki/beamforming/beams/vertical_beam.err
#$ -o /Volumes/GeoPhysics_05/users-data/yannik78/taranaki/beamforming/laura_beam_test/beams/vertical_beam.out
#$ -e /Volumes/GeoPhysics_05/users-data/yannik78/taranaki/beamforming/laura_beam_test/beams/vertical_beam.err

export PROC_SRC=/Users/home/yannik78/dev/proc-scripts_git/
#LIST=/Volumes/GeoPhysics_05/users-data/yannik78/taranaki/beamforming/beams/dirlist_april.txt
#LIST=/Volumes/GeoPhysics_05/users-data/yannik78/taranaki/beamforming/beams/list_all_months.txt
LIST=/Volumes/GeoPhysics_05/users-data/yannik78/taranaki/beamforming/laura_beam_test/beams/four_month_list.txt
ROOT=/Volumes/GeoPhysics_05/users-data/yannik78/taranaki
DIR=$(cat $LIST | head -n $SGE_TASK_ID | tail -n 1)
INPUT=${ROOT}/${DIR}
/usr/local/python2/bin/python2.7 /Users/home/yannik78/dev/proc-scripts_git/beamforming/beamforming.py  $INPUT -s -d -b --nstat=20


