#!/bin/sh
#detection-submit.sh
#Torque script to run Matlab program

#Torque directives
#PBS -N psc-detection-$param_ind
#PBS -W group_list=yetistats
#PBS -l nodes=1:ppn=12,walltime=04:00:00,mem=16000mb
#PBS -M shababo@berkeley.edu
#PBS -m abe
#PBS -V

#set output and error directories (SSCC example here)
#PBS -o localhost:/vega/stats/users/bms2156/log
#PBS -e localhost:/vega/stats/users/bms2156/log

#Command to execute Matlab code
matlab -nosplash -nodisplay -nodesktop -r "detect_pscs('$trace_file','$param_file', $param_ind, $noise_type)" > matoutfile.$PBS_ARRAYID

#Command below is to execute Matlab code for Job Array (Example 4) so that each part writes own output
#matlab -nosplash -nodisplay -nodesktop -r "simPoissGLM($LAMBDA)" > matoutfile.$PBS_ARRAYID

#End of script
