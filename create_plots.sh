#!/bin/bash

dest=${1%/}

option="dataset=\"$dest\""

gnuplot -e $option plot/plot_dCP.glp
gnuplot -e $option plot/plot_M23.glp
gnuplot -e $option plot/plot_S23.glp
gnuplot -e $option plot/plot_S13.glp

gnuplot -e $option plot/cont_dCP_M23.glp
gnuplot -e $option plot/cont_S13_dCP.glp
gnuplot -e $option plot/cont_S23_dCP.glp
gnuplot -e $option plot/cont_S13_M23.glp
gnuplot -e $option plot/cont_S23_M23.glp
gnuplot -e $option plot/cont_S13_S23.glp

gnuplot -e $option plot/sens.glp
gnuplot -e $option plot/sens_detail.glp

echo mv -f $dest"*.*" plot/$dest/
mv -f $dest*.* plot/$dest/
