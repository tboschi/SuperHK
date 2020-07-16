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

model=("stats/asim/NH_NH" "new_nuenorm5_corr/asim/NH_NH" "0/asim/NH_NH")
name="validate_all"
point="point_118595"

echo Extract chi2 information
$Drop -p -f $point -n "$name" -r "errorstudy/" "${model[@]}"

echo Create plot
$Plot "$name"
