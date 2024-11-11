#### Copy file(s) from siam-linux20 to lavended

# Copy folder_coord
tic()
server <- "/Volumes/siam-linux20/lavended/projects/patter-flapper"
from   <- file.path(server, iteration$folder_coord)
to     <- iteration$folder_coord
# dirs.copy(from, to, cl = 2L)
toc()