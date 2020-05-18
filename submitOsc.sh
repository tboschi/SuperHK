#!/bin/bash

NCPU=$3

QUEUE_NAME=atmpd
BINDIR=/home/tboschi/OscAna/hk.atm+beam
BINNAME=GlobalOsc

CARDDIR=/home/tboschi/OscAna/hk.atm+beam/Cards
CARD=$1
CARDNAME=${CARD##*/}
OUTDIR=$2

#if [ ! -e ${OUTDIR} ] then
mkdir -p ${OUTDIR}
#fi

SUBRUN=0

FILENAME=${BINNAME}_${CARDNAME}

for ((SUBRUN=0; SUBRUN<${NCPU};SUBRUN++))
do
  SHELLNAME=${OUTDIR}/R${FILENAME}_$SUBRUN.sh

  if [ -e ${SHELLNAME} ]; then
    #echo "find" ${SHELLNAME} "rm"
    rm ${SHELLNAME}
  fi

  OUTPUT=${OUTDIR}/${FILENAME}_$SUBRUN.log
  ERROR=${OUTDIR}/${FILENAME}_$SUBRUN.err

  cat << EOF > ${SHELLNAME}
#!/bin/bash
source /home/tboschi/OscAna/Setup_OscAna.sh
cd ${BINDIR}
date
./bin/${BINNAME} ${SUBRUN} ${NCPU} ${OUTDIR} $CARD
date
EOF

  echo -r job${SUBRUN} -o ${OUTPUT} -e ${ERROR} ${SHELLNAME} >> /home/tboschi/jobManager/${QUEUE_NAME}.list

done
