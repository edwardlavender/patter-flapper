#!/bin/bash

# Update
# * Script batch size is too small
# * Some batches killed (due to memory usage?)
# * Killed batches re-run (see deploy-simulations-2.sh)
# * For future runs, run larger batches (fewer CPUs)

Rscript ./R/010-simulate-algorithms.R 1095.0 "1:10" &
Rscript ./R/010-simulate-algorithms.R 1095.0 "11:20" &
Rscript ./R/010-simulate-algorithms.R 1095.0 "21:30" &
Rscript ./R/010-simulate-algorithms.R 1095.0 "31:40" &
Rscript ./R/010-simulate-algorithms.R 1095.0 "41:50" &
Rscript ./R/010-simulate-algorithms.R 1095.0 "51:60" &
Rscript ./R/010-simulate-algorithms.R 1095.0 "61:70" &
Rscript ./R/010-simulate-algorithms.R 1095.0 "71:80" &
Rscript ./R/010-simulate-algorithms.R 1095.0 "81:90" &
Rscript ./R/010-simulate-algorithms.R 1095.0 "91:100" &
Rscript ./R/010-simulate-algorithms.R 1095.0 "101:110" &
Rscript ./R/010-simulate-algorithms.R 1095.0 "111:120" &
Rscript ./R/010-simulate-algorithms.R 1095.0 "121:130" &
Rscript ./R/010-simulate-algorithms.R 1095.0 "131:140" &
Rscript ./R/010-simulate-algorithms.R 1095.0 "141:150" &
Rscript ./R/010-simulate-algorithms.R 1095.0 "151:160" &
Rscript ./R/010-simulate-algorithms.R 1095.0 "161:170" &
Rscript ./R/010-simulate-algorithms.R 1095.0 "171:180" &
Rscript ./R/010-simulate-algorithms.R 1095.0 "181:190" &
Rscript ./R/010-simulate-algorithms.R 1095.0 "191:200" &
Rscript ./R/010-simulate-algorithms.R 1095.0 "201:210" &
Rscript ./R/010-simulate-algorithms.R 1095.0 "211:220" &
Rscript ./R/010-simulate-algorithms.R 1095.0 "221:230" &
Rscript ./R/010-simulate-algorithms.R 1095.0 "231:240" &
Rscript ./R/010-simulate-algorithms.R 1095.0 "241:250" &
Rscript ./R/010-simulate-algorithms.R 1095.0 "251:260" &
Rscript ./R/010-simulate-algorithms.R 1095.0 "261:270" &
Rscript ./R/010-simulate-algorithms.R 1095.0 "271:280" &
Rscript ./R/010-simulate-algorithms.R 1095.0 "281:290" &
Rscript ./R/010-simulate-algorithms.R 1095.0 "291:300" &
Rscript ./R/010-simulate-algorithms.R 1095.0 "301:310" &
Rscript ./R/010-simulate-algorithms.R 1095.0 "311:320" &
Rscript ./R/010-simulate-algorithms.R 1095.0 "321:330" &
Rscript ./R/010-simulate-algorithms.R 1095.0 "331:340" &
Rscript ./R/010-simulate-algorithms.R 1095.0 "341:350" &
Rscript ./R/010-simulate-algorithms.R 1095.0 "351:360" &
Rscript ./R/010-simulate-algorithms.R 1095.0 "361:370" &
Rscript ./R/010-simulate-algorithms.R 1095.0 "371:380" &
Rscript ./R/010-simulate-algorithms.R 1095.0 "381:390" &
Rscript ./R/010-simulate-algorithms.R 985.5 "1:10" &
Rscript ./R/010-simulate-algorithms.R 985.5 "11:20" &
Rscript ./R/010-simulate-algorithms.R 985.5 "21:27" &
Rscript ./R/010-simulate-algorithms.R 1204.5 "1:10" &
Rscript ./R/010-simulate-algorithms.R 1204.5 "11:20" &
Rscript ./R/010-simulate-algorithms.R 1204.5 "21:27" &
wait
echo "All simulations completed!"