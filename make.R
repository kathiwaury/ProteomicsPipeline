source("R/packages.R")  # loads packages
source("R/functions.R") # defines the functions
source("R/plan.R")      # creates the drake plan

make(
  plan,
  verbose = TRUE,
  lock_envir = FALSE
)


