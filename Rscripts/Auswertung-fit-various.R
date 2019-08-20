GH_BA_dir <- normalizePath("~/ethz/BA/")
save_dir  <- normalizePath("~/BScThesis/4MM/fit-various")
##--------------------------------------------------
stopifnot(dir.exists(GH_BA_dir),
          dir.exists(save_dir))


devtools::load_all(file.path(GH_BA_dir, "norMmix"))

data(SMI.12, package="copula")

## str(SMI.12)
#  num [1:141, 1:20] 16.1 15.7 15.7 16.1 16.6 ...
#  - attr(*, "dimnames")=List of 2
#   ..$ : chr [1:141] "2011-09-09" "2011-09-12" "2011-09-13" "2011-09-14" ...
#   ..$ : chr [1:20] "ABBN" "ATLN" "ADEN" "CSGN" ...
# NULL

## parcond(SMI.12, k=10, model="VVV")
# [1] 0.0610654
## very bad
## parcond(SMI.12, k=2, model="EII")
# [1] 3.357143
## no way to get "good" i.e. >5 values


####
##--------------------------------------------------
####
## datasets: SMI.12, iris
