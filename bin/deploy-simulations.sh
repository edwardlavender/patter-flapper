#!/bin/bash

# Record start time
start_time=$(date)
echo "Start time: $start_time"
start_seconds=$(date +%s)

# For 1095.0, in batches of 25
# * The number of batches is relatively low to handle memory requirements
Rscript ./R/010-simulate-algorithms.R 1095.0 "1:25" &
Rscript ./R/010-simulate-algorithms.R 1095.0 "26:50" &
Rscript ./R/010-simulate-algorithms.R 1095.0 "51:75" &
Rscript ./R/010-simulate-algorithms.R 1095.0 "76:100" &
Rscript ./R/010-simulate-algorithms.R 1095.0 "101:125" &
Rscript ./R/010-simulate-algorithms.R 1095.0 "126:150" &
Rscript ./R/010-simulate-algorithms.R 1095.0 "151:175" &
Rscript ./R/010-simulate-algorithms.R 1095.0 "176:200" &
Rscript ./R/010-simulate-algorithms.R 1095.0 "201:225" &
Rscript ./R/010-simulate-algorithms.R 1095.0 "226:250" &
Rscript ./R/010-simulate-algorithms.R 1095.0 "251:275" &
Rscript ./R/010-simulate-algorithms.R 1095.0 "276:300" &
Rscript ./R/010-simulate-algorithms.R 1095.0 "301:325" &
Rscript ./R/010-simulate-algorithms.R 1095.0 "326:350" &
Rscript ./R/010-simulate-algorithms.R 1095.0 "351:375" &
Rscript ./R/010-simulate-algorithms.R 1095.0 "376:390" &

# For 985.5, one batch stopping at 27
Rscript ./R/010-simulate-algorithms.R 985.5 "1:27" &

# For 1204.5, one batch stopping at 27
Rscript ./R/010-simulate-algorithms.R 1204.5 "1:27" &

# Wait for processes to complete
wait

# Print end time & duration (hours, mins, secs)
end_seconds=$(date +%s)
duration=$((end_seconds - start_seconds))
hours=$((duration / 3600))
minutes=$(((duration % 3600) / 60))
seconds=$((duration % 60))
echo "Duration: $hours hours, $minutes minutes, $seconds seconds"

echo " Simulations completed!"