#!/bin/bash

EXEC=~/OscAna/Osc3++/processing/build.contours/BuildContourPlots

while IFS='' read -r file || [[ -n "$file" ]]; do
	echo "Contouring: $file"
	#echo $EXEC $file
	$EXEC $file	&> /dev/null
	file=${file/SpaghettiSens./contour.}
	file=${file/.000000./.}
	echo saving to $file
	mv ChiSquared.root $file
done < "$1"
