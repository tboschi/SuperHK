file="file=\"hkdr_prediction/app_FHC.dat\""
outp="outp=\"hkdr_app_fhc.tex\""
gnuplot -e $file -e $outp spectra_app.glp

file="file=\"asim_prediction/app_FHC.dat\""
outp="outp=\"asim_app_fhc.tex\""
gnuplot -e $file -e $outp spectra_app.glp

file="file=\"hkdr_prediction/disapp_FHC.dat\""
outp="outp=\"hkdr_dis_fhc.tex\""
gnuplot -e $file -e $outp spectra_dis.glp

file="file=\"asim_prediction/disapp_FHC.dat\""
outp="outp=\"asim_dis_fhc.tex\""
gnuplot -e $file -e $outp spectra_dis.glp

file="file=\"hkdr_prediction/app_RHC.dat\""
outp="outp=\"hkdr_app_rhc.tex\""
gnuplot -e $file -e $outp spectra_app.glp

file="file=\"asim_prediction/app_RHC.dat\""
outp="outp=\"asim_app_rhc.tex\""
gnuplot -e $file -e $outp spectra_app.glp

file="file=\"hkdr_prediction/disapp_RHC.dat\""
outp="outp=\"hkdr_dis_rhc.tex\""
gnuplot -e $file -e $outp spectra_dis.glp

file="file=\"asim_prediction/disapp_RHC.dat\""
outp="outp=\"asim_dis_rhc.tex\""
gnuplot -e $file -e $outp spectra_dis.glp
