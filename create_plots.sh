#!/bin/bash

dest=${1%/}

option="dataset=\"$dest\""

echo option $option

cd plot

gnuplot -e "$option" plot_dCP.gpl
gnuplot -e "$option" plot_M23.gpl
gnuplot -e "$option" plot_S23.gpl
gnuplot -e "$option" plot_S13.gpl

gnuplot -e "$option" cont_dCP_M23.gpl
gnuplot -e "$option" cont_S13_dCP.gpl
gnuplot -e "$option" cont_S23_dCP.gpl
gnuplot -e "$option" cont_S13_M23.gpl
gnuplot -e "$option" cont_S23_M23.gpl
gnuplot -e "$option" cont_S13_S23.gpl

gnuplot -e "$option" sens.gpl
gnuplot -e "$option" sens_detail.gpl

echo mv -f $dest"*.*" $dest/
mv -f $dest*.* $dest/

cd ..

./bin/plot_bundle plot/$dest/$dest'_'*.tex
