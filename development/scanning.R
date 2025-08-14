n_atoms <- 18L
n_lengths <- n_atoms - 1L
n_angles <- n_atoms - 2L
n_dihedrals <- n_atoms - 3L

lengths_bounds <- c(0.9, 5.7)
angles_bounds <- c(0.52, 3.05)
dihedral_bounds <- c(-pi, pi)

create_data <- function(n, zeros_up_to, boundaries) {
  grid_dense <- 10
  lapply(1:n, function(x) {
    print(x)
    if (x <= zeros_up_to) {
      rep(0.0, grid_dense)
    } else {
      seq(boundaries[1], boundaries[2], length.out = grid_dense)
    }
  })
}

lengths <- create_data(n_atoms, 1, lengths_bounds)
angles <- create_data(n_atoms, 2, angles_bounds)
dihedrals <- create_data(n_atoms, 3, dihedral_bounds)

df <- lapply(1:n_atoms, function(x) {
  temp <- do.call(expand.grid, list(lengths[[x]], angles[[x]], dihedrals[[x]]))
  temp[!duplicated(temp), ]
})

# Too many combinations
