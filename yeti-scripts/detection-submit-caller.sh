#!/bin/bash

#trace_file="/vega/stats/users/bms2156/psc-detection/data/psc_traceset_test.mat"
#param_file="/vega/stats/users/bms2156/psc-detection/data/psc_parameters_test01.mat"
#noise_type=1

 for f in "$dir"/*; do echo "$f"; done

dir = "/vega/stats/users/bms2156/psc-detection/data/cluster-param-files"
for param_file in "$dir"/*; #486
do

	response=`qsub -v param_file=$param_file /vega/stats/users/bms2156/psc-detection/yeti-scripts/detection-submit.sh`
	echo $response
	sleep 5 
done
