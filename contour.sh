#! /bin/bash

Contour=/data/tboschi/HKsens/OscAna/SuperHK/bin/buildcontours
root=/data/tboschi/HKsens/

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

for ll in "${list[@]}" ; do
	dir=${ll%/}
	dir=${dir##*/}

	echo file $ll, $dir

	mkdir -p $cont/$dir

	echo find $sens/$dir/ -name $name.*.root
	nfiles=$(find $sens/$dir/ -name $name.*.root | wc -l)
	echo $nfiles
	if [ $nfiles -eq '0' ] ; then
		continue
	fi

	#is there an updated version of sens files?
	updt=$(ls -t $sens/$dir/$name.*.root | head -n1)
	tupdt=$(date -r $updt +%s)

	#check when last contour was done
	echo filtering $tr
	if [ "$tr" = true ] ; then
		last=$cont/$dir/$outp'_filter.root'
	else
		last=$cont/$dir/$outp'.root'
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
		hadd -f $cont/$dir/all.root $sens/$dir/$name.*.root
		$Contour $cont/$dir/all.root $cont/$dir/$outp.root $card $sens/$dir/this_sensitivity.card
		#mv ChiSquared.root 

		#if [ "$tr" = true ] ; then
		#	$Filter $cont/$dir/all.$ff.root $cont/$dir/all_filter.$ff.root
		#	echo $Contour $cont/$dir/all_filter.$ff.root >& /dev/null
		#	$Contour $cont/$dir/all_filter.$ff.root >& /dev/null
		#	mv ChiSquared.root $cont/$dir/$outp'_filter'.$ff.root
		#fi

		#rm -f $cont/$dir/all.root
	fi
done

echo "DONE"
