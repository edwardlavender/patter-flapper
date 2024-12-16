#!/bin/bash

# Record start time
start_time=$(date)
echo "Start time: $start_time"
start_seconds=$(date +%s)

# For 1095.0, in batches of 50, stopping at 528
Rscript ./R/017-run-patter.R 1095.0 "1:50" &
Rscript ./R/017-run-patter.R 1095.0 "51:100" &
Rscript ./R/017-run-patter.R 1095.0 "101:150" &
Rscript ./R/017-run-patter.R 1095.0 "151:200" &
Rscript ./R/017-run-patter.R 1095.0 "201:250" &
Rscript ./R/017-run-patter.R 1095.0 "251:300" &
Rscript ./R/017-run-patter.R 1095.0 "301:350" &
Rscript ./R/017-run-patter.R 1095.0 "351:400" &
Rscript ./R/017-run-patter.R 1095.0 "401:450" &
Rscript ./R/017-run-patter.R 1095.0 "451:500" &
Rscript ./R/017-run-patter.R 1095.0 "501:528" &

# For 985.5, in batches of 50, stopping at 144
Rscript ./R/017-run-patter.R 985.5 "1:50" &
Rscript ./R/017-run-patter.R 985.5 "51:100" &
Rscript ./R/017-run-patter.R 985.5 "101:144" &

# For 1204.5, in batches of 50, stopping at 144
Rscript ./R/017-run-patter.R 1204.5 "1:50" &
Rscript ./R/017-run-patter.R 1204.5 "51:100" &
Rscript ./R/017-run-patter.R 1204.5 "101:144" &

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