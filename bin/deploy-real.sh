#!/bin/bash

# Record start time
start_time=$(date)
echo "Start time: $start_time"
start_seconds=$(date +%s)

# For 1095.0, in balanced batches, stopping at 528
Rscript ./R/017-run-patter.R 1095.0 "1:88" &
Rscript ./R/017-run-patter.R 1095.0 "89:176" &
Rscript ./R/017-run-patter.R 1095.0 "177:264" &
Rscript ./R/017-run-patter.R 1095.0 "265:352" &
Rscript ./R/017-run-patter.R 1095.0 "353:440" &
Rscript ./R/017-run-patter.R 1095.0 "441:528" &

# For 985.5, in batches of ~72, stopping at 144
Rscript ./R/017-run-patter.R 985.5 "1:72" &
Rscript ./R/017-run-patter.R 985.5 "73:144" &

# For 1204.5, in batches of ~72, stopping at 144
Rscript ./R/017-run-patter.R 1204.5 "1:72" &
Rscript ./R/017-run-patter.R 1204.5 "73:144" &

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
