#### Copy file(s) from siam-linux20 to lavended

# Define server path
server <- "/Volumes/lavended/documents/projects/patter-flapper"
stopifnot(dir.exists(server))

# Define iteration
iteration <- copy(iteration_patter)

# Copy folder_coord
# * ~3.5 mins (sims, 153 rows)
tic()
from    <- file.path(server, iteration$folder_coord)
to      <- iteration$folder_coord
# success <- dirs.copy(from[1], to[1], cl = NULL)
success <- dirs.copy(from, to, cl = 10L)
toc()

# Review success
table(success$success)

# 2024-11-12 simulations copied from server