#!/bin/bash

# Instructions to run the forward simulation on siam-linux20

# ---------- SIA-LAVENDED-M test ---------- 

# 1. Navigate to RStudio Project
cd /Users/lavended/Documents/work/packages/flapper/patter-flapper/patter-flapper/

# 2. Run code 
R CMD BATCH ./R/010-run-forward-sim-1.R ./data/output/forward/log.txt

# ---------- siam-linux20 ---------- 

# 1. Copy RStudio Project onto server into documents/projects/patter-flapper

# 2. ssh into server
ssh lavended@siam-linux20 

# 3. Update code
cd documents/projects
git pull https://github.com/edwardlavender/patter-flapper.git

# 3. Open RStudio Project interactively to create project library via {renv}
cd patter-flapper
R
Sys.setenv(GITHUB_PAT = "<insert-pat-here>")
t1 <- Sys.time()
renv::restore()
t2 <- Sys.time()
difftime(t2, t1)
q("no")

# 3. Run the forward simulation

R CMD BATCH ./R/010-run-forward-sim-1.R ./data/output/forward/log.txt