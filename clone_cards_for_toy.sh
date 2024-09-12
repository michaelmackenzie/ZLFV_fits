#! /bin/bash
# Clone a combine data card collection for use in a toy
# Useage: ./clone_cards_for_toy.sh <input dir> <output dir> <toy file>

DATACARDS=$1
OUTDIR=$2
TOYFILE1=$3
TOYFILE2=$4

if [[ "${DATACARDS}" == "" ]]; then
    echo "No input collection given!"
    exit
fi

if [[ "${OUTDIR}" == "" ]]; then
    echo "No output collection name given!"
    exit
fi

if [[ "${TOYFILE1}" == "" ]]; then
    echo "No input toy file given!"
    exit
fi

if [[ "${TOYFILE2}" == "" ]]; then
    echo "No input toy file given!"
    exit
fi

[ ! -d ${OUTDIR} ] && mkdir -p ${OUTDIR}

if [[ "${DATACARDS}" != *"/" ]]; then
    DATACARDS="${DATACARDS}/"
fi
if [[ "${OUTDIR}" != *"/" ]]; then
    OUTDIR="${OUTDIR}/"
fi

#FIXME: Generalize this linking
if [ ! -e ${OUTDIR}WorkspaceScanTOY ]; then
    cd ${OUTDIR}
    ln -s ../../WorkspaceScanTOY WorkspaceScanTOY
    cd -
fi
if [ ! -e ${OUTDIR}WorkspaceScanSGN ]; then
    cd ${OUTDIR}
    ln -s ../../WorkspaceScanSGN WorkspaceScanSGN
    cd -
fi
if [ ! -e ${OUTDIR}WorkspaceScanBKG ]; then
    cd ${OUTDIR}
    ln -s ../../WorkspaceScanBKG WorkspaceScanBKG
    cd -
fi

cp ${DATACARDS}datacard_zprime_*_mp*.txt ${OUTDIR}

for CARD in `ls -d ${OUTDIR}*_0d3_0d7_*.txt`
do
    # Create binned toy datasets for each input mass point
    echo ${CARD}
    DATAWS1=`cat ${CARD} | grep "shapes background" | awk '{print $4}'` #Retrieve the background workspace name from the datacard
    DATAWS2=`echo ${DATAWS1} | sed "s/0d3_0d7/0d7_1d0/g"`
    MASSPOINT=$(echo $DATAWS1 | tr "_" " " | awk '{print $NF}' | sed 's/.root//')
    TOY_MP_1=`echo ${TOYFILE1} | sed "s/.root/_${MASSPOINT}.root/"`
    TOY_MP_2=`echo ${TOYFILE2} | sed "s/.root/_${MASSPOINT}.root/"`
    python clone_data_binning.py --data ${DATAWS1} --toy ${TOYFILE1} -o ${TOY_MP_1}
    python clone_data_binning.py --data ${DATAWS2} --toy ${TOYFILE2} -o ${TOY_MP_2}


    # Replace the data file names with the toy file names
    cat ${CARD} | awk -v toy1=${TOY_MP_1} -v toy2=${TOY_MP_2} '{
                                                                 if($2 == "data_obs" && $4 ~ /0d3_0d7/){print $1,$2,$3,toy1,$5}
                                                                 else if($2 == "data_obs" && $4 ~ /0d7_1d0/){print $1,$2,$3,toy2,$5}
                                                                 else print $0}' >| tmp.txt
    cp tmp.txt ${CARD}
    rm tmp.txt
    CARD=`echo ${CARD} | sed 's/0d3_0d7/0d7_1d0/g'`
    cat ${CARD} | awk -v toy1=${TOY_MP_1} -v toy2=${TOY_MP_2} '{
                                                                 if($2 == "data_obs" && $4 ~ /0d3_0d7/){print $1,$2,$3,toy1,$5}
                                                                 else if($2 == "data_obs" && $4 ~ /0d7_1d0/){print $1,$2,$3,toy2,$5}
                                                                 else print $0}' >| tmp.txt
    cp tmp.txt ${CARD}
    rm tmp.txt
    CARD=`echo ${CARD} | sed 's/_bdt_0d7_1d0//g'`
    cat ${CARD} | awk -v toy1=${TOY_MP_1} -v toy2=${TOY_MP_2} '{
                                                                 if($2 == "data_obs" && $4 ~ /0d3_0d7/){print $1,$2,$3,toy1,$5}
                                                                 else if($2 == "data_obs" && $4 ~ /0d7_1d0/){print $1,$2,$3,toy2,$5}
                                                                 else print $0}' >| tmp.txt
    cp tmp.txt ${CARD}
    rm tmp.txt

done
