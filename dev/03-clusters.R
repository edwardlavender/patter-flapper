# Set up the workflow on other machines

# Patter.jl installation -------------------------------------------------------

# Windows: set up the number of threads
Threads.nthreads()

# Install Patter.jl from the DEV branch
julia
using Pkg
Pkg.add(url = "https://github.com/edwardlavender/Patter.jl#dev")


# patter installation ----------------------------------------------------------

# Clone patter, navigate to dev branch & install locally 
# Or, istall Patter from the DEV branch
renv::install("edwardlavender/patter@dev")

# Define Julia Home
julia_home <- NULL
julia_home <- "/Applications/Julia-1.10.app/Contents/Resources/julia/bin/"

# Set up patter (test)
library(patter)
julia_connect(JULIA_HOME = julia_home)


# Checks -----------------------------------------------------------------------

# Patter.ModelMove should contain comments about mobility:
julia_help("Patter.ModelMove")

# Run example code
# ?pf_filter


# patter-flapper set up --------------------------------------------------------

# Set up patter-flapper version control branch
# https://github.com/edwardlavender/patter-flapper.git

# Navigate to DEV branch
# * Delete Julia/ folder contents

# Add .Renviron (Windows)
JULIA_PROJ = "C:/Users/lavended/Documents/projects/patter-flapper/Julia"
OMP_NUM_THREADS = 8
JULIA_NUM_THREADS = 8

# Add .Renviron (M1 MacBook)
JULIA_HOME = "/Applications/Julia-1.10.app/Contents/Resources/julia/bin/"
JULIA_PROJ = "/Users/el72/Documents/projects/patter-flapper"
OMP_NUM_THREADS = 12
JULIA_NUM_THREADS = 12

# Project-specific installation of Julia dependencies
library(patter)
julia_connect()

# Update Patter.jl to the DEV branch
# * This is necessary b/c julia_connect() installs from the main branch
# > cd("Julia")
# > ]
# > activate .
# > add https://github.com/edwardlavender/Patter.jl#dev

# Restart RStudio & validate Patter.jl installation 
library(JuliaCall)
library(patter)
julia_connect()
julia_help("ModelMove")

# Directories
# * Copy data/input/ (input.zip)
# * Copy data/input/spatial (spatial.zip, via One Drive)
# * Create data/output/analysis
# * Create data/output/simulation

