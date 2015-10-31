#!/bin/bash

#trace_file="/vega/stats/users/bms2156/psc-detection/data/psc_traceset_test.mat"
#param_file="/vega/stats/users/bms2156/psc-detection/data/psc_parameters_test01.mat"
#noise_type=1


dirname="/vega/stats/users/bms2156/psc-detection/data"
for results_file in "$dirname"/simulated-epscs-1027*; #486
do

	response=`qsub -v results_file=$results_file /vega/stats/users/bms2156/psc-detection/yeti-scripts/time-scoring-submit.sh`
	echo $results_file
	echo $response
	sleep 5 
done
