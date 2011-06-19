#!/bin/sh
#
#$ -S /bin/sh
#$ -t 1-203
#$ -q all.q
#$ -wd /Volumes/GeoPhysics_05/users-data/yannik78/taranaki/beamforming/beams
#$ -m n
#$ -N vertical_beamforming
#$ -o /Volumes/GeoPhysics_05/users-data/yannik78/taranaki/beamforming/beams/vertical_beam.out
#$ -e /Volumes/GeoPhysics_05/users-data/yannik78/taranaki/beamforming/beams/vertical_beam.err

export PROC_SRC=/Users/home/yannik78/dev/proc-scripts_git/
#LIST=/Volumes/GeoPhysics_05/users-data/yannik78/taranaki/beamforming/beams/dirlist_april.txt
LIST=/Volumes/GeoPhysics_05/users-data/yannik78/taranaki/beamforming/beams/list_all_months.txt
ROOT=/Volumes/GeoPhysics_05/users-data/yannik78/taranaki
DIR=$(cat $LIST | head -n $SGE_TASK_ID | tail -n 1)
INPUT=${ROOT}/${DIR}
/usr/local/python2/bin/python2.7 /Users/home/yannik78/dev/proc-scripts_git/beamforming/beamforming.py  $INPUT -s -b -d --nstat=20


