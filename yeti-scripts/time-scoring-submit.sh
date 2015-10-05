#!/bin/sh
#detection-submit.sh
#Torque script to run Matlab program

#Torque directives
#PBS -N psc-detection
#PBS -W group_list=yetistats
#PBS -l nodes=1:ppn=12,walltime=00:30:00,mem=8000mb
#PBS -M shababo@berkeley.edu
#PBS -m abe
#PBS -V

#set output and error directories (SSCC example here)
#PBS -o localhost:/vega/stats/users/bms2156/log
#PBS -e localhost:/vega/stats/users/bms2156/log

#Command to execute Matlab code
#echo $trace_file
#echo $param_file
#echo $param_ind
#echo $noise_type
resultsfile_base=${results_file##*/}
resultsfile_base=${resultsfile_base%.*}
echo $resultsfile_base
matlab -nosplash -nodisplay -nodesktop -r "score_detection('$results_file', 20, 0)" > matlab-output/matoutfile-${PBS_JOBID}-${resultsfile_base}

#Command below is to execute Matlab code for Job Array (Example 4) so that each part writes own output
#matlab -nosplash -nodisplay -nodesktop -r "simPoissGLM($LAMBDA)" > matoutfile.$PBS_ARRAYID

#End of script
