#!/bin/bash

Drop=./dropchi2.sh

#model=("stats/asim/NH_NH" "nuenorm1_corr/asim/NH_NH" "nuenorm2_corr/asim/NH_NH" "nuenorm3_corr/asim/NH_NH" "nuenorm4_corr/asim/NH_NH" "nuenorm5_corr/asim/NH_NH")
#name="validate_nuenorm_corr"
#point="point_49773"

model=("stats/asim/NH_NH" "nuenorm1_anti/asim/NH_NH" "nuenorm2_anti/asim/NH_NH" "nuenorm3_anti/asim/NH_NH" "nuenorm4_anti/asim/NH_NH" "nuenorm5_anti/asim/NH_NH")
name="validate_nuenorm_anti"
point="point_49773"

#model=("stats/asim/NH_NH" "new_nuenorm5_corr/asim/NH_NH" "new_nuenorm5_anti/asim/NH_NH")
#name="validate_new_nuenorm"
#point="point_96927"

echo $Drop -f $point -n $name -r errorstudy/ "${model[@]}"
$Drop -p -f $point -n "$name" -r "errorstudy/" "${model[@]}"
