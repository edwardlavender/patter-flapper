#!/bin/bash

# Instructions to run the forward simulation on siam-linux20

# 1. Copy RStudio Project onto server into documents/projects/patter-flapper

# 2. ssh into server
ssh lavended@siam-linux20 

# 3. Update code
cd documents/projects
git pull https://github.com/edwardlavender/patter-flapper.git

# 3. Open RStudio Project interactively to create project library via {renv}
cd patter-flapper
R
q("no")

# 3. Run the forward simulation
R CMD BATCH R/010-run-forward-sim.R >>> data/output/forward/log.txt