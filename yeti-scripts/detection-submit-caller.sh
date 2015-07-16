#!/bin/bash

trace_file = "/vega/stats/users/bms2156/psc-detection/data/psc_traceset_01.mat"
param_file = "/vega/stats/users/bms2156/psc-detection/data/psc_parameters_test01.mat"
noise_type = 3

for i in `seq 1 4`; #486
do
	echo $i
	echo $trace_file
	echo $param_file
	echo $noise_type
	
	response=`qsub -v trace_file=trace_file,param_file=param_file,param_ind=i,noise_type=noise_type /vega/stats/users/bms2156/psc-detection/yeti-scripts/detection-submit.sh`
	echo $response
	sleep 5 
done

