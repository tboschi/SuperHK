#! /bin/bash

Filter=/data/tboschi/HKsens/OscAna/SuperHK/bin/filter
Contour=/data/tboschi/HKsens/OscAna/SuperHK/bin/BuildContourPlots
root=/data/tboschi/HKsens/
samples=(T2HK)

pp=false
tr=false
while getopts 'fr:p' flag; do
	case "${flag}" in
		f) tr=true ;;
		r) root="${OPTARG}" ;;
		p) pp=true ;;
		*) exit 1 ;;
	esac
done

shift $((OPTIND-1))

root=${root%/}
sens=$root/sensitivity
cont=$root/contours

list=()
for point in "$@"
do
	list=("${list[@]}" "$sens/$point")
done

#if no array of points is provided, take all points fitted!
if [ ${#list[@]} -eq 0 ] ; then
	list=($(ls -d $sens/point_*/))
fi

name=""
outp=""
if [ "$pp" = true ] ; then
	name="SpaghettiSens_penalised"
	outp="gaussian"
else
	name="SpaghettiSens"
	outp="uniform"
fi

for ff in "${samples[@]}" ; do
	echo Processing $ff set
	for ll in "${list[@]}" ; do
		dir=${ll%/}
		dir=${dir##*/}

		echo file $ll, $dir

		mkdir -p $cont/$dir

		echo find $sens/$dir/ -name $name.$ff.*.root
		nfiles=$(find $sens/$dir/ -name $name.$ff.*.root | wc -l)
		echo $nfiles
		if [ $nfiles -eq '0' ] ; then
			continue
		fi

		#is there an updated version of sens files?
		updt=$(ls -t $sens/$dir/$name.$ff.*.root | head -n1)
		tupdt=$(date -r $updt +%s)

		#check when last contour was done
		echo filtering $tr
		if [ "$tr" = true ] ; then
			last=$cont/$dir/$outp'_filter'.$ff.root
		else
			last=$cont/$dir/$outp.$ff.root
		fi
		tlast=0
		if [ -s $last ] ; then
			tlast=$(date +%s -r $last)
		fi
		echo $last, $updt
		echo $tlast, $tupdt

		if [ $tupdt -gt $tlast ] ; then
			echo "Creating" $outp "contour for" $dir

			#contour of pure files
			hadd -f $cont/$dir/all.$ff.root $sens/$dir/$name.$ff.*.root
			echo $Contour $cont/$dir/all.$ff.root >& /dev/null
			$Contour $cont/$dir/all.$ff.root >& /dev/null
			mv ChiSquared.root $cont/$dir/$outp.$ff.root

			if [ "$tr" = true ] ; then
				$Filter $cont/$dir/all.$ff.root $cont/$dir/all_filter.$ff.root
				echo $Contour $cont/$dir/all_filter.$ff.root >& /dev/null
				$Contour $cont/$dir/all_filter.$ff.root >& /dev/null
				mv ChiSquared.root $cont/$dir/$outp'_filter'.$ff.root
			fi

			rm -f $cont/$dir/all*.$ff.root
		fi

		#nfiles=$(find $root/$dir/ -name "SpaghettiSens_penalised.$ff.*.root" | wc -l)
		#if [ $nfiles -eq '0' ] ; then
		#	continue
		#fi

		#last=$(ls -t $sens/$dir/SpaghettiSens_penalised.$ff.*.root | head -n1)
		#updt=$(ls -t $cont/$dir/gaussian.$ff.root | head -n2)
		#tlast=$(date -r $last +%s)
		#tupdt=$(date -r $last +%s)
		#if [ updt -gt last ] ; then
		#	echo "Creating gaussian contour for" $dir
		#	#contour of penalised files
		#	hadd -f $cont/$dir/all.$ff.root $sens/$dir/SpaghettiSens_penalised.$ff.*.root
		#	$Contour $cont/$dir/all.$ff.root >& /dev/null
		#	mv ChiSquared.root $cont/$dir/gaussian.$ff.root
		#	rm $cont/$dir/all.$ff.root
		#fi
	done
done

echo "DONE"
