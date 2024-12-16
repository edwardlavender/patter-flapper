#!/bin/bash

# Record start time
start_time=$(date)
echo "Start time: $start_time"
start_seconds=$(date +%s)

# For 1095.0, in batches of 30, stopping at 528
Rscript ./R/017-run-patter.R 1095.0 "1:30" &
Rscript ./R/017-run-patter.R 1095.0 "31:60" &
Rscript ./R/017-run-patter.R 1095.0 "61:90" &
Rscript ./R/017-run-patter.R 1095.0 "91:120" &
Rscript ./R/017-run-patter.R 1095.0 "121:150" &
Rscript ./R/017-run-patter.R 1095.0 "151:180" &
Rscript ./R/017-run-patter.R 1095.0 "181:210" &
Rscript ./R/017-run-patter.R 1095.0 "211:240" &
Rscript ./R/017-run-patter.R 1095.0 "241:270" &
Rscript ./R/017-run-patter.R 1095.0 "271:300" &
Rscript ./R/017-run-patter.R 1095.0 "301:330" &
Rscript ./R/017-run-patter.R 1095.0 "331:360" &
Rscript ./R/017-run-patter.R 1095.0 "361:390" &
Rscript ./R/017-run-patter.R 1095.0 "391:420" &
Rscript ./R/017-run-patter.R 1095.0 "421:450" &
Rscript ./R/017-run-patter.R 1095.0 "451:480" &
Rscript ./R/017-run-patter.R 1095.0 "481:510" &
Rscript ./R/017-run-patter.R 1095.0 "511:528" &

# For 985.5, in batches of 30, stopping at 144
Rscript ./R/017-run-patter.R 985.5 "1:30" &
Rscript ./R/017-run-patter.R 985.5 "31:60" &
Rscript ./R/017-run-patter.R 985.5 "61:90" &
Rscript ./R/017-run-patter.R 985.5 "91:120" &
Rscript ./R/017-run-patter.R 985.5 "121:144" &

# For 1204.5, in batches of 30, stopping at 144
Rscript ./R/017-run-patter.R 1204.5 "1:30" &
Rscript ./R/017-run-patter.R 1204.5 "31:60" &
Rscript ./R/017-run-patter.R 1204.5 "61:90" &
Rscript ./R/017-run-patter.R 1204.5 "91:120" &
Rscript ./R/017-run-patter.R 1204.5 "121:144" &

# Wait for processes to complete
wait

# Print end time & duration (hours, mins, secs)
end_seconds=$(date +%s)
duration=$((end_seconds - start_seconds))
hours=$((duration / 3600))
minutes=$(((duration % 3600) / 60))
seconds=$((duration % 60))
echo "Duration: $hours hours, $minutes minutes, $seconds seconds"

echo "Simulations completed!"
