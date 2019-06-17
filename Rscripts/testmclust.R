### test script to learn mclust and mixtools package

# data() lists all available example data sets

library(mclust) # version 5.4.3
library(mixtools) #version 1.1.0

data(diabetes)

class <- diabetes$class
table(class)
# class
# Chemical   Normal    Overt 
#       36       76       33 

X <- diabetes[,-1]

clPairs(X, class)

x9 <- rnorMix(MW.nm9, n=5000)

y9 <- Mclust(x9)
summary(y9, par=TRUE)
# ---------------------------------------------------- 
# Gaussian finite mixture model fitted by EM algorithm 
# ---------------------------------------------------- 
# Mclust E (univariate, equal variance) model with 3 components: 
#  log-likelihood    n df      BIC       ICL
#       -7846.796 5000  6 -15744.7 -17150.79
# Clustering table:
#    1    2    3 
# 1978 1019 2003 
# Mixing probabilities:
#         1         2         3 
# 0.3850980 0.2224116 0.3924904 
# Means:
#           1           2           3 
# -1.34250220 -0.05709185  1.31779230 
# Variances:
#         1         2         3 
# 0.2874259 0.2874259 0.2874259 

#somewhat off, how to controll iterations


	z9 <- em(modelName="V", x9,   ) # needs init param

emctl <- emControl(itmax=50000)
y9 <- Mclust(x9, control= emctl)
summary(y9, par=TRUE) # still not accurate at itmax 50000. 
# maybe stuck at loc. opt.



# normalmixEM from mixtools package

library(nor1mix)

x9 <- rnorMix(MW.nm9, n=500)
mixem <- normalmixEM(x9)
mixem$mu
# [1] -1.1039705  0.9953695
mixem$sigma
# [1] 0.6942218 0.7528372
#only found 2 clusters
plot(mixem$all.loglik, type="l") #stagnated after 50 iterations

x9 <- rnorMix(MW.nm9, n=5000) #give more data
mixem$mu
# [1]  1.035134 -1.120661
mixem$sigma
# [1] 0.7204488 0.6556280
#compare to actual:
MW.nm9
# 'Normal Mixture' object 	 ``#9 Trimodal'' 
#        mu sigma    w
# [1,] -1.2  0.60 0.45
# [2,]  1.2  0.60 0.45
# [3,]  0.0  0.25 0.10
plot(mixem$all.loglik, type="l") #similar results



## comparing packages


