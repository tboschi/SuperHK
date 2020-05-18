#!/bin/bash

for i in {0..200}
do
	echo test point $i
	./submitSens.sh $i
	#./submitSensStatOnly.sh $i
done
