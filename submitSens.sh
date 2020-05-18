#!/bin/bash

QUEUE_NAME=atmpd
NFILE=$1
NCPU=$1

CARDDIR=/home/tboschi/OscAna/hk.atm+beam/Cards
BINDIR=/home/tboschi/OscAna/hk.atm+beam
BINNAME=GlobalSpaghettiSensitivity

CARD=$2
CARD2=$3
TPOINT=$4
OUTDIR=$5
CARDNAME=${CARD##*/}
#OUTDIR=/home/tboschi/data/asimov_nh_ih.valor/spaghetti_t$TPOINT
#CARD2=scan.asimov_ih
#OUTDIR=/home/tboschi/data/hkdr.large/spaghetti_t$TPOINT
#CARD=scan.hkdr
#CARD2=$CARD

echo $OUTDIR

SUBRUN=0

#if( -e $OUTDIR) then
#mv $OUTDIR ${OUTDIR}_`date +"%y%m%d%H%M"`
#endif
mkdir -p ${OUTDIR}

#if($#argv < 1) then
#echo please input comment
#exit
#endif


#echo $1 >> $OUTDIR/comment

FILENAME=${BINNAME}_${CARDNAME}

#set CPLIST = (0 9 18)	# 3 loop
#set NCP = 13
#set NS23 = 9
#set S23LIST = (0 2 4 6 8)	# 1 loop
#set NM23 = 9
#set M23LIST = (0)
#set NS13 = 9
#set S13LIST = (0)

for ((SUBRUN=0; SUBRUN<${NFILE};SUBRUN++))
  do
#@ IPOINT = ${SUBRUN} / ${NCPU}
#@ CPP = 1 + ( ${IPOINT} / ${#S23LIST} / ${#M23LIST} / ${#S13LIST} ) % ${#CPLIST} 
#@ M23P = 1 +  ( ${IPOINT} / ${#S23LIST} / ${#S13LIST} ) % ${#M23LIST}
#@ S23P = ${IPOINT} % ${#S23LIST} + 1
#@ S13P = ( ${IPOINT} / ${#S23LIST} ) % ${#S13LIST} + 1
#@ TPOINT = $CPLIST[$CPP] * ${NS23} * ${NS13} * ${NM23} + $M23LIST[$M23P] * ${NS23} * ${NS13} + $S13LIST[$S13P] * ${NS23} + $S23LIST[$S23P]

#------------Original---------------------------------
# S23=0.4 ---------------  
#  TPOINT=360 #dCP=0
#  TPOINT=6921 #dCP=90
#  TPOINT=13482 #dCP=180
#  TPOINT=20043 #dCP=270

# S23=0.45 ---------------  
#  TPOINT=362 #dCP=0
#  TPOINT=6923 #dCP=90
#  TPOINT=13484 #dCP=180
#  TPOINT=20045 #dCP=270

# S23=0.5 ---------------  
#  TPOINT=364 #dCP=0
#  TPOINT=6925 #dCP=90
#  TPOINT=13486 #dCP=180
#  TPOINT=20047 #dCP=270

# S23=0.55 ---------------  
#  TPOINT=366 #dCP=0
#  TPOINT=6927 #dCP=90
#  TPOINT=13488 #dCP=180
#  TPOINT=20049 #dCP=270

# S23=0.6 ---------------  
#  TPOINT=368 #dCP=0
#  TPOINT=6929 #dCP=90
#  TPOINT=13490 #dCP=180
#  TPOINT=20051 #dCP=270
#------------Original---------------------------------


#------------Wider---------------------------------
# S23=0.6 ---------------  
#  TPOINT=530 #dCP=0
  #TPOINT=19484 #dCP=180
#------------Wider---------------------------------

  LSUBRUN=`expr ${SUBRUN} % ${NCPU}`

  SHELLNAME=${OUTDIR}/R${FILENAME}_${TPOINT}_${LSUBRUN}.sh
  OUTPUT=${OUTDIR}/${FILENAME}_${TPOINT}_${LSUBRUN}.log
  #OUTPUT=/dev/null
  ERROR=${OUTDIR}/${FILENAME}_${TPOINT}_${LSUBRUN}.err
  #ERROR=/dev/null
  if [ -e ${SHELLNAME} ]; then
      echo "find" ${SHELLNAME} "rm"
      rm ${SHELLNAME}
  fi

  cat > ${SHELLNAME} <<EOF
#!/bin/bash
source /home/tboschi/OscAna/Setup_OscAna.sh
cd ${BINDIR}
./bin/${BINNAME} $LSUBRUN $NCPU ${OUTDIR}/ ${CARD} ${CARD2} ${TPOINT}
EOF

  echo -r job${SUBRUN} -o ${OUTPUT} -e ${ERROR} ${SHELLNAME} >> /home/tboschi/jobManager/${QUEUE_NAME}.list
done

