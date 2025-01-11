#### Copy file(s) from siam-linux20 to lavended

library(tictoc)

#### Define server path

# siam-linux20
# * Real: ACDC (all batches) [DONE], AC (batch 2), AC (batch 3), DC (batch 2), DC (batch 3) [DONE]
# server <- "/Volumes/homes/documents/projects/patter-flapper"
server <- "/Volumes/lavended/documents/projects/patter-flapper"
# server   <- "/Users/lavended/Desktop/server"

# workstation
# * Real: AC (batch 1) [DONE]
# server <- "/Users/lavended/Desktop/workstation"

# macbook
# * Real: DC (batch 1) [DONE]
# server <- "/Users/lavended/Desktop/macbook"

stopifnot(dir.exists(server))

#### Define iteration
# iteration <- copy(iteration_patter)
# iteration <- iteration[dataset %in% c("ac", "dc") & mobility %in% pars$pmovement$mobility[2:3], ]
nrow(iteration)
head(iteration[, .(folder_coord)])

#### Time trial: copy folder_coord[1]
# Simulations: 22.819 s -> 22.819 * nrow(iteration)/60 = 169 mins
from    <- file.path(server, iteration$folder_coord)
to      <- iteration$folder_coord
head(cbind(from, to))
tic()
success <- dirs.copy(from[1], to[1], cl = NULL)
toc()

#### Copy folder_coord
# Simulations: 2.83 hr (Fritzbox)
tic()
success <- dirs.copy(from, to, cl = 12L)
toc()

#### (optional) Keep connection to server active
# This code should be run in a separate R session
for (i in 1:1e6) {
  message("-------------------------------------")
  print(paste0(i, ": ", Sys.time()))
  tic()
  file.create(file.path(server, "data-raw", "server", paste0(i, ".txt")))
  toc()
  Sys.sleep(60 * 5)
}
unlink(list.files(file.path(server, "data-raw", "server"), full.names = TRUE))

#### Review success
table(success$success)

#### End of code. 