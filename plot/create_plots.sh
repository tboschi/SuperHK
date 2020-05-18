#!/bin/bash


dest=${1%/}

option="dataset=\"$dest\""

gnuplot -e $option plot_dCP.glp
gnuplot -e $option plot_M23.glp
gnuplot -e $option plot_S23.glp
gnuplot -e $option plot_S13.glp

gnuplot -e $option cont_dCP_M23.glp
gnuplot -e $option cont_S13_dCP.glp
gnuplot -e $option cont_S23_dCP.glp
gnuplot -e $option cont_S13_M23.glp
gnuplot -e $option cont_S23_M23.glp
gnuplot -e $option cont_S13_S23.glp

gnuplot -e $option sens.glp
gnuplot -e $option sens_detail.glp

echo mv -f $dest"*.*" $dest/
mv -f $dest*.* $dest/
