#!/bin/bash

Drop=./dropchi2.sh
Plot=./create_plots.sh

#model=("stats/asim/NH_NH" "nuenorm5_corr/asim/NH_NH" "nuenorm5_anti/asim/NH_NH")
#name="validate_small_corr"
#point="point_49773"

#model=("stats/asim/NH_NH" "nuenorm1_anti/asim/NH_NH" "nuenorm2_anti/asim/NH_NH" "nuenorm3_anti/asim/NH_NH" "nuenorm4_anti/asim/NH_NH" "nuenorm5_anti/asim/NH_NH")
#name="validate_nuenorm_anti"
#point="point_49773"

#model=("stats/asim/NH_NH" "new_nuenorm5_corr/asim/NH_NH" "new_nuenorm5_anti/asim/NH_NH")
#name="validate_large_nuenorm"
##point="point_96927"
#point="point_358443"

#model=("stats/asim/NH_NH" "new_nuenorm5_corr/asim/NH_NH" "0/asim/NH_NH")
#model=("stats/asim/NH_NH" "0/asim/NH_NH" "11a/asim/NH_NH" "11b/asim/NH_NH")
#name="validate_scale"
#point="point_118595"

#model=("stats/asim/NH_NH" "new_nuenorm5_corr/asim/NH_NH" "0/asim/NH_NH")
model=("stats_fine/asim/NH_NH" "nominal_fine/asim/NH_NH" "nominal_energyscale_fine/asim/NH_NH" "nominal_energyscale_fine_2/asim/NH_NH" )
name="fine_scale_fix"
point="point_81209"



parms=("M23" "S13" "S23" "CP")
card="errorstudy/${model[0]}/sensitivity/$point/this_sensitivity.card"
lims="plot/limits.gpl"
echo "pi = 4.*atan(1.)" > "$lims"
for p in "${parms[@]}"
do
	echo "" >> $lims
	cat "$card" | awk -v par=$p '/^parm/ && $0~par {print "x0_" par " = " $2}' >> $lims
	cat "$card" | awk -v par=$p '/^parm/ && $0~par {print "x1_" par " = " $3}' >> $lims
	cat "$card" | awk -v par=$p '/^parm/ && $0~par {print "nn_" par " = " $4}' >> $lims
done


echo Extract chi2 information
$Drop -p -f $point -n "$name" -r "errorstudy/" "${model[@]}"

echo Create plot
$Plot "$name"
