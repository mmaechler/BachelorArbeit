### test script to learn mclust package
## version 5.4.3

# data() lists all available example data sets

library(mclust)

data(diabetes)

class <- diabetes$class
table(class)
# class
# Chemical   Normal    Overt 
#       36       76       33 

X <- diabetes[,-1]

clPairs(X, class)
