#!/bin/sh
#
#$ -S /bin/sh
#$ -t 1-163
#$ -q all.q
#$ -wd /Volumes/Wanaka_01/yannik/start/beamforming/
#$ -m n
#$ -N vertical_beamforming

export PROC_SRC=/Users/home/yannik78/dev/proc-scripts_git/
LIST=/Volumes/Wanaka_01/yannik/start/beamforming/list_all_months_start.txt
ROOT=/Volumes/Wanaka_01/yannik/start
DIR=$(cat $LIST | head -n $SGE_TASK_ID | tail -n 1)
INPUT=${ROOT}/${DIR}
/usr/local/python2/bin/python2.7 /Users/home/yannik78/dev/proc-scripts_git/beamforming/beamforming.py  $INPUT -s


