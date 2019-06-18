#! /bin/bash

Oscillator=/home/tboschi/OscAna/SuperHK/bin/oscillator
cardS=/home/tboschi/OscAna/SuperHK/cards

global=/home/tboschi/data/
root=/home/tboschi/data/
MH=NH

#global contains already asim or hkdr

while getopts 'ins:g:' flag; do
	case "${flag}" in
		i) MH=IH ;;
		n) MH=NH ;;
		g) global="${OPTARG}" ;;
		*) exit 1 ;;
	esac
done

#where unoscillated files are
input=${global%/}/$MH/global
output=${global%/}/$MH/oscillated

echo $input
echo $output

#create files	#needs to be done once only

mkdir -p $output


cp $cardS/t2hkk_.card $cardS/t2hkk_oscillated.card

sed -i "s:input.*:input \"$input/*.root\":"	"$cardS/t2hkk_oscillated.card"
sed -i "s:output.*:output \"$output\":"		"$cardS/t2hkk_oscillated.card"

if [[ $MH = "NH" ]]; then
	sed -i "s:hierarchy.*:hierarchy \"normal\":"	"$cardS/t2hkk_oscillated.card"
else
	sed -i "s:hierarchy.*:hierarchy \"inverted\":"	"$cardS/t2hkk_oscillated.card"
fi

#oscillate beam part, updating cards with new systematics
#must need to openmp this one or at least paralellise
echo $Oscillator $cardS/t2hkk_oscillated.card

echo mv $cardS/t2hkk_oscillated.card $output/oscillated.card
