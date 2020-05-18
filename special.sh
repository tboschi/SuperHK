#!/bin/bash

echo "printing hostname"
hostname


echo
echo "running exe"
#/data/tboschi/HKsens/OscAna/SuperHK/bin/momspaghetti $1 $2 $3
/data/tboschi/HKsens/OscAna/SuperHK/bin/fitter $1 $2 $3
