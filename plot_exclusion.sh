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
	oscc=${oscc/contours/sensitivity}/oscillation.card

	echo "pi = 4.*atan(1.)" > $olim
	for p in "${parms[@]}"
	do
		echo "" >> $lims
		cat "$oscc" | awk -v par=$p '/^parm/ && $0~par {gsub('/,/', "", $2); print "x0_" par " = " $2}' >> $olim
		cat "$oscc" | awk -v par=$p '/^parm/ && $0~par {gsub('/,/', "", $3); print "x1_" par " = " $3}' >> $olim
		cat "$oscc" | awk -v par=$p '/^parm/ && $0~par {gsub('/,/', "", $4); print "nn_" par " = " $4}' >> $olim
	done

	cp $olim $lims
	if [ -f $lims ] && ! diff $lims $olim; then # lims file does not exist
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
