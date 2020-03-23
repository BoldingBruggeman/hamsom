#!/bin/csh

#BASE EUR NAM ASGM IND INTW POWERGEN
set EXEC = sbatchecosmo.pbs 

set EDATE = 1990 
setenv IYE $EDATE
echo "sbatch ./${EXEC} "
set ID = `sbatch ./${EXEC} | awk '{print $4}'`

while ($EDATE < 1990) #2015365 for full year
  sleep 1
  @ EDATE = $EDATE + 1
  setenv IYE $EDATE
  echo "sbatch --dependency=afterok:${ID} ./${EXEC}"
  set ID = `sbatch --dependency=afterok:${ID} ./${EXEC} | awk '{print $4}'`
end

