#!/bin/bash

#trace_file="/vega/stats/users/bms2156/psc-detection/data/psc_traceset_test.mat"
#param_file="/vega/stats/users/bms2156/psc-detection/data/psc_parameters_test01.mat"
#noise_type=1


<<<<<<< HEAD
dirname="/vega/stats/users/bms2156/psc-detection/data/param-files/evoked-test-1218"
for param_file in "$dirname"/*; #486
do

=======
#dirname="/vega/stats/users/bms2156/psc-detection/data/cluster-param-files"
#for param_file in "$dirname"/*-pspike-*.mat; #486
#do
param_file=0
>>>>>>> af8819ce92e8cf52e507a0933e6ce8d8d7bef76d
	response=`qsub -v param_file=$param_file /vega/stats/users/bms2156/psc-detection/yeti-scripts/detection-submit.sh`
#	echo $param_file
	echo $response
#	sleep 5 
#done
