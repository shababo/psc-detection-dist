#!/bin/bash

trace_file = "/vega/stats/users/bms2156/psc-detection/data/psc_traceset_01.mat"
param_file = "/vega/stats/users/bms2156/psc-detection/data/psc_parameters_test01.mat"


for i in `seq 1 486`;
do
	echo $i
	response=`qsub -v trace_file=trace_file,param_file=param_file,param_ind=i /vega/stats/users/bms2156/psc-detection/yeti-scripts/detection-submit.sh`
	echo $response
	sleep 5 
done

