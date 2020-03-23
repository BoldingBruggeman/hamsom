#!/bin/sh

# script to start a single ECOSMO run
# that's why runs can be started directly with this file 

# written by: Cara Nissen (uda081@student.uib.no), ute.daewel@hzg.de

echo -e " \n Preparations and start of new model run. \n"

#define model specifications here 
#runID_old run Id that should be replaced by the new runID

#runSP run specification "ECOSMOCO2" fo simulation with carbon chemistry
#rm *.o
Mrate_old=2
runSP="ECOSMO"
#runSP="ECOSMOCO2"
runFI="NOFISH"
#runFI="FISH"
runID_old=cr1
runID=ttt

#set compiler options
#DSATLAS : start biology from WOA values
#DSASCII : start physics from ascii climatology
#else remember to define file in code (here n011978)
#DECO2 : use carbon module (check inpur files for Alkalinity and DIC)
#DNOFISH : use model without fish and MB 
#DWADDEN : consider specific conditions for wadden sea boundary
#DCARA  : consider CARA parameterisation of sediment remineralizeation under hypoxic conditions
#export OPT="-DCOASTDAT -DCARA"
#export OPT="-DSASCII -DCARA -DNCEP"
#export OPT="-DSATLAS -DSASCII -DCARA -DCOASTDAT"
export OPT="-DSATLAS -DSASCII -DCOASTDAT"
#export OPT="-DCARA -DNOFISH"
#export OPT="-DSATLAS -DSASCII"
#export OPT="-DSASCII"
#export OPT="-DSEDPO"
#export OPT="-DPUMP"

#if [ ${runFI} ==  "NOFISH" ]; then
#cp ECOSMparam_mi0.f ECOSMparam_${runID}.f
#else
#cp ECOSMparam_nsp.f ECOSMparam_${runID}.f
#fi



echo -e "Die Run ID ist $runID \n"

################
# Below this block, nothing needs to be changed
################

start=$(cat /work/uda081/ecosmo/input_$runID | head -n1)
duration=$(cat /work/uda081/ecosmo/input_$runID | tail -n1)
let end_year=$duration-1
let end=$start+$end_year

echo -e "Model run starts in $start and model is run for $duration years. \n"

years=$(seq $start $end)

for i in $years 
do
  if [ -d /work/uda081/ecosmo/north_b/f_out_$i ]; then
    echo -e "Folder "f_out_$i" exists" 
  else
     mkdir /work/uda081/ecosmo/north_b/f_out_$i
     echo -e "Folder "f_out_$i" didn't exist, but was created."
  fi
done
echo -e "All output folders exist. \n"


# create folder with model scripts for runID

# Update sbatchecosmo.pbs
sed -i s/ECOSM${runID_old}/ECOSM${runID}/g  sbatchecosmo.pbs
sed -i s/ECOSMO_${runID_old}/ECOSMO_${runID}/g  sbatchecosmo.pbs
sed -i s/output${runID_old}/output${runID}/g sbatchecosmo.pbs

# Update sbatchecosmo_co2.pbs

# Update Makefile
sed -i s/ECOSMO_${runID_old}/ECOSMO_${runID}/g Makefile
#sed -i s/bio_${runID_old}/bio_${runID}/g Makefile

#cp Makefile Makefile_${runID}

# Update main_new.F
#sed -i 's/ppp='${runID_old}'/ppp='${runID}'"/g' main_new.F
#sed -i s/ECOSMparam_${runID_old}/ECOSMparam_${runID}/g main_new.F
#cp bio_${runID_old}.F bio_${runID}.F
#sed -i s/ECOSMparam_${runID_old}/ECOSMparam_${runID}/g bio_${runID}.F
#sed -i s/${runID_old}/${runID}/g bio_${runID}.F

# Update main_new.F
if [ ${runSP} ==  "ECOSMOCO2" ]; then
sed -i s/"(nbio=nbiox)"/"(nbio=nbiox+5)"/g C_model.f
else 
sed -i s/"(nbio=nbiox+5)"/"(nbio=nbiox)"/g C_model.f
fi

echo -e "sbatchecosmo.pbs, Makefile and main_new.F are updated. \n"

echo -e "All needed changes made. All files and folders created."
echo -e "Model run can be started. \n"
mkdir track${runID}
cp Makefile track${runID}
cp compile_ECOSMO.sh track${runID}

make
#sbatch sbatchecosmo.pbs
runID_old=${runID}

#done
 
cd

