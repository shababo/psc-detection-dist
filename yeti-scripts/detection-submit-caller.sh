#!/bin/bash

trace_file="/vega/stats/users/bms2156/psc-detection/data/psc_traceset_test.mat"
param_file="/vega/stats/users/bms2156/psc-detection/data/psc_parameters_test01.mat"
noise_type=3

for i in `seq 1 486`; #486
do

	response=`qsub -v trace_file=$trace_file,param_file=$param_file,param_ind=$i,noise_type=$noise_type /vega/stats/users/bms2156/psc-detection/yeti-scripts/detection-submit.sh`
	echo $response
	sleep 5 
done

