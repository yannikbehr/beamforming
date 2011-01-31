#!/bin/sh
#
#$ -S /bin/sh
#$ -t 1-30
#$ -q all.q
#$ -wd /Volumes/Wanaka_01/yannik/taranaki/beamforming/beams
#$ -m n
#$ -N vertical_beamforming

export PROC_SRC=/Users/home/yannik78/dev/proc-scripts_git/
LIST=/Volumes/Wanaka_01/yannik/taranaki/beamforming/beams/dirlist_april.txt
ROOT=/Volumes/Wanaka_01/yannik/taranaki/sacfiles/5Hz/2002/Apr
DIR=$(cat $LIST | head -n $SGE_TASK_ID | tail -n 1)
INPUT=${ROOT}/${DIR}
#echo "/usr/local/python2/bin/python2.7 /Users/home/yannik78/dev/proc-scripts_git/beamforming.py  ${INPUT}"
#/usr/local/python2/bin/python2.7 /Users/home/yannik78/dev/proc-scripts_git/beamforming/hbeamforming.py  $INPUT -s
/usr/local/python2/bin/python2.7 /Users/home/yannik78/dev/proc-scripts_git/beamforming/beamforming.py  $INPUT -s


