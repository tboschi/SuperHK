#!/bin/bash

#############################
#### user defined options ###
#############################


# output name unique to combination of sets
name="one_two_three"

# list of files with CPV scan results
model=("errorstudy/one/NH_NH/contours/CPV_scan.dat"
       "errorstudy/two/NH_NH/contours/CPV_scan.dat"
       "errorstudy/three/NH_NH/contours/CPV_scan.dat")

# names for displaying above results
nickn=("one_cpv"
       "two_cpv"
       "three_cpv")




###################################
#### do not modify script below ###
###################################


Sens=bin/dropsens

lims="plot/limits.gpl"
olim="plot/.limits.gpl"
parms=("M12" "M23" "S12" "S13" "S23" "CP")

echo Creating dir plot/$name
mkdir -p plot/$name
card=plot/$name/$name.card

rm $lims $olim

for i in ${!model[*]}; do 
	m="${model[$i]}"
	oscc=${m%/*}
	oscc=${m/contours/sensitivity}/oscillation.card

	echo "pi = 4.*atan(1.)" > $olim
	for p in "${parms[@]}"
	do
		echo "" >> $lims
		cat "$card" | awk -v par=$p '/^parm/ && $0~par {print "x0_" par " = " $2}' >> $olims
		cat "$card" | awk -v par=$p '/^parm/ && $0~par {print "x1_" par " = " $3}' >> $olims
		cat "$card" | awk -v par=$p '/^parm/ && $0~par {print "nn_" par " = " $4}' >> $olims
	done
	if [ ! -f $lims ] ; then # lims file does not exist
		cp $olim $lims
	elif ! diff $lims $olim ; then
		echo oscillation cards are different! Cannot process "$m" with the rest
		exit 1
	fi

	echo "${nickn[$i]}" \"$m\" >> $card
done

echo Extract exclusion information
$Sens $card


if [ $name ] ; then
	ssfile=($(ls exclusion_*.dat))

	for f in ${ssfile[@]}; do
		echo $f
		mv $f plot/$name/$name'_'$f
	done
fi

option="dataset=\"$name\""

plot=( "sens_all.gpl"
       "sens_diff.gpl" )

cd plot

for p in "${plot[@]}" ; do
	gnuplot -e "$option" $p
done

echo mv -f $name"*.*" $name/
mv -f $name*.* $name/

cd ..

./bin/plot_bundle plot/$name/$name'_'*.tex
