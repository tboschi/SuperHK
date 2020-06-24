app () {
	echo gnuplot -e $1 -e $2 -e $3 -e $4 spectra_app.gpl
	gnuplot -e "$1" -e "$2" -e "$3" -e "$4" spectra_app.gpl
}

dis () {
	echo gnuplot -e $1 -e $2 -e $3 -e $4 spectra_dis.gpl
	gnuplot -e "$1" -e "$2" -e "$3" -e "$4" spectra_dis.gpl
}


file='file="hkdr_prediction/app_FHC.dat"'
outp='outp="hkdr_app_fhc.tex"'
posc='posc="Design"'
mode='mode="\$\\nu\$-mode"'
app "$file" "$outp" "$posc" "$mode"


file='file="asim_prediction/app_FHC.dat"'
outp='outp="asim_app_fhc.tex"'
posc='posc="Asimov A"'
mode='mode="\$\\nu\$-mode"'
app "$file" "$outp" "$posc" "$mode"


file='file="hkdr_prediction/disapp_FHC.dat"'
outp='outp="hkdr_dis_fhc.tex"'
posc='posc="Design"'
mode='mode="\$\\nu\$-mode"'
dis "$file" "$outp" "$posc" "$mode"


file='file="asim_prediction/disapp_FHC.dat"'
outp='outp="asim_dis_fhc.tex"'
posc='posc="Asimov A"'
mode='mode="\$\\nu\$-mode"'
dis "$file" "$outp" "$posc" "$mode"


file='file="hkdr_prediction/app_RHC.dat"'
outp='outp="hkdr_app_rhc.tex"'
posc='posc="Design"'
mode='mode="\$\\cj\{\\nu\}\$-mode"'
app "$file" "$outp" "$posc" "$mode"


file='file="asim_prediction/app_RHC.dat"'
outp='outp="asim_app_rhc.tex"'
posc='posc="Asimov A"'
mode='mode="\$\\cj\{\\nu\}\$-mode"'
app "$file" "$outp" "$posc" "$mode"


file='file="hkdr_prediction/disapp_RHC.dat"'
outp='outp="hkdr_dis_rhc.tex"'
posc='posc="Design"'
mode='mode="\$\\cj\{\\nu\}\$-mode"'
dis "$file" "$outp" "$posc" "$mode"


file='file="asim_prediction/disapp_RHC.dat"'
outp='outp="asim_dis_rhc.tex"'
posc='posc="Asimov A"'
mode='mode="\$\\cj\{\\nu\}\$-mode"'
dis "$file" "$outp" "$posc" "$mode"
