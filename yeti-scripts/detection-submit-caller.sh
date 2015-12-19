#!/bin/bash

#trace_file="/vega/stats/users/bms2156/psc-detection/data/psc_traceset_test.mat"
#param_file="/vega/stats/users/bms2156/psc-detection/data/psc_parameters_test01.mat"
#noise_type=1


dirname="/vega/stats/users/bms2156/psc-detection/data/param-files/evoked-test-1218"
for param_file in "$dirname"/*; #486
do
#param_file=0
	response=`qsub -v param_file=$param_file /vega/stats/users/bms2156/psc-detection/yeti-scripts/detection-submit.sh`
	echo $param_file
	echo $response
	sleep 5 
#done
