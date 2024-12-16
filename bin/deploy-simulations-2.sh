#!/bin/bash
Rscript ./R/010-simulate-algorithms.R 1095.0 "131:140" &
Rscript ./R/010-simulate-algorithms.R 1095.0 "151:160" &
Rscript ./R/010-simulate-algorithms.R 1095.0 "191:200" &
Rscript ./R/010-simulate-algorithms.R 1095.0 "201:210" &
Rscript ./R/010-simulate-algorithms.R 1095.0 "211:220" &
Rscript ./R/010-simulate-algorithms.R 1095.0 "231:240" &
Rscript ./R/010-simulate-algorithms.R 1095.0 "241:250" &
Rscript ./R/010-simulate-algorithms.R 1095.0 "301:310" &
Rscript ./R/010-simulate-algorithms.R 1095.0 "331:340" &
Rscript ./R/010-simulate-algorithms.R 1095.0 "361:370" &
Rscript ./R/010-simulate-algorithms.R 1095.0 "381:390" &
Rscript ./R/010-simulate-algorithms.R 1204.5 "21:27" &
wait
echo "All simulations completed!"