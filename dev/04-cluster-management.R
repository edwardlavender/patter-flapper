#### Copy file(s) from siam-linux20 to lavended

#### Define server path

# siam-linux20
# * Real: ACDC (all batches) [DONE], AC (batch 2), AC (batch 3), DC (batch 2), DC (batch 3) [TO DO]
# server <- "/Volumes/homes/documents/projects/patter-flapper"
# server <- "/Volumes/lavended/documents/projects/patter-flapper"
server   <- "/Users/lavended/Desktop/server"

# workstation
# * Real: AC (batch 1) [DONE]
# server <- "/Users/lavended/Desktop/workstation"

# macbook
# * Real: DC (batch 1) [DONE]
# server <- "/Users/lavended/Desktop/macbook"

stopifnot(dir.exists(server))

#### Define iteration
# iteration <- copy(iteration_patter)
nrow(iteration)
head(iteration[, .(folder_coord)])

#### Copy folder_coord
# * ~12 mins (sims, 444 rows)
tic()
from    <- file.path(server, iteration$folder_coord)
to      <- iteration$folder_coord
head(cbind(from, to))
# success <- dirs.copy(from[1], to[1], cl = NULL)
success <- dirs.copy(from, to, cl = 12L)
toc()

#### Review success
table(success$success)

# 2024-11-16 simulations copied from server