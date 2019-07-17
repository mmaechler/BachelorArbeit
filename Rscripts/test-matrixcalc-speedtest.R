

set.seed(2019)

TT <- matrix(sample(1:247),19,13)


TT[19,] <- 1


D. <- diag(1:13)

TT %*% D.

TT * rep((1:13), each=19)

### testing multiplication

system.time({for (i in 1:100000) {TT%*%D.}})
#    user  system elapsed 
#   0.204   0.000   0.204 


system.time({for (i in 1:100000) {TT*rep((1:13), each=19)}})
#    user  system elapsed 
#   0.816   0.004   0.820 

system.time({for (i in 1:100000) {TT%*%diag(1:13)}})
#    user  system elapsed 
#   0.396   0.000   0.396 


### testing division

system.time({for (i in 1:100000) {TT/rep((1:13), each=19)}})
#    user  system elapsed 
#   0.884   0.000   0.884 

system.time({for (i in 1:100000) {TT%*%diag(1/(1:13))}})
#    user  system elapsed 
#   0.388   0.000   0.387 


### testing with difficult floating point numbers

TT <- TT*pi
system.time({for (i in 1:100000) {TT*rep((1:13), each=19)}})
#    user  system elapsed 
#   0.836   0.000   0.844 

system.time({for (i in 1:100000) {TT%*%diag(1:13)}})
#    user  system elapsed 
#   0.444   0.000   0.444 
