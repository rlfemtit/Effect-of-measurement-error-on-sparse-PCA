

set.seed(1)
library(lsa)

# Case1 : no sparsity
# Case2 : with sparsity
load("Psparse2.RData")
load("variables.RData")
Out = list(X = out$X, P =out$P, W = out$W, Z= out$Z, nzeros = out$nzeros,
           k = out$k,Propsparse= out$Propsparse )

X = Out$X # 1000 X 100 matrix / generated data : Z%*%t(P) + Error
P = Out$P # 100 X 5 matrix / 0.2 sparsity : 30% loadings per each variable are zero.
Z = Out$Z # 1000 X 5 matrix / true score

#rownames(P) = paste0("X", 1:dim(X)[2])
#colnames(P) = paste0("PC", 1:2)

X_svd = svd(X)
d = X_svd$d # singular values
u = X_svd$u # 10000 x 10 matrix / each column has the unit length
v = X_svd$v # 10 x 10 matrix / each column is the eigen vector

# data with high sparsity
v[,1] = c(0.422, 0.422, 0.422, 0.422, 0, 0, 0, 0, 0.380, 0.380)
v[,2] = c(0 , 0, 0, 0, 0.489, 0.489, 0.489, 0.489, -0.147, 0.147)

d[1] = d[1] * 1.5 # make distinct eigen values

d = rep(1,100)
d[1] = 5

X = u %*% diag(d) %*% t(v)


!!!!!!!!!!!!!!!!!!! median angle과 angle distribution을 만들자... 하지말고 


## sample loading
loading_sample_data = list()
loading_x1 = c()
loading_x2 = c()
loading_x3 = c()
n = 5
a= 1
b= 5
c= 9

for (i in 1:1000){
  j = 0
  
  # sampling
  sample_loading = svd(X[sample(nrow(X),size=n,replace=FALSE),])$v[,1:2]
  
  # Fixing the sign of sample loadings using cosine similiarity
  j=1
  while (j<3){
    cosine_original = cosine(v[,j], sample_loading[,j])
    cosine_flipped = cosine(v[,j], -sample_loading[,j])
    
    if (cosine_original < cosine_flipped){
      sample_loading[,j] = -sample_loading[,j]
    }
    j=j+1
  }

  loading_sample_data[[i]] = sample_loading
  loading_x1 = append(sample_loading[a,1], loading_x1)
  loading_x2 = append(sample_loading[b,1], loading_x2)
  loading_x3 = append(sample_loading[c,1], loading_x3)
}
expected_loading = Reduce("+",loading_sample_data) / length(loading_sample_data) # the expected loadings
rownames(expected_loading) = paste0("X", 1:10)
colnames(expected_loading) = paste0("PC", 1:2)

# angle between expected and true loading
v1 = expected_loading[,1]
v2 = v[,1]

180 * acos( sum(v1*v2) / ( sqrt(sum(v1 * v1)) * sqrt(sum(v2 * v2)) ) ) / pi 


## loading plot
# true loading
par(mfrow=c(1,2))
plot(v[1:10,1:2],   # x and y data
     pch=21,              # point shape
     bg="black",          # point color
     cex=1,               # point size
     xlim=c(-0.7, 0.7),
     ylim=c(-0.6, 0.6),
     main="True loadings"      # title of plot
)
abline(a=0, b=0, h=0, v=0, lty=2)

text(P[1:10,1:2],             # sets position of labels
     labels=rownames(P[1:10,]),   # print labels
     cex=0.9,
     pos=3
)

# sample loading
plot(expected_loading[1:10,1:2],   # x and y data
     pch=21,              # point shape
     bg="black",          # point color
     cex=1,               # point size
     xlim=c(-0.7, 0.7),
     ylim=c(-0.6, 0.6),
     main=paste0("Sample loadings (n=",n,")")      # title of plot
)
abline(a=0, b=0, h=0, v=0, lty=2)

text(expected_loading[1:10,1:2],             # sets position of labels
     labels=rownames(expected_loading[1:10,]),   # print labels
     cex=0.9,
     pos=3
)

## histogram of sample loadings
# a
par(mfrow=c(1,3))
mean(loading_x1)
hist(loading_x1, breaks=40, main="Histogram of loading (X1)", xlab="X1 loading")
abline(v = v[a,1], # true loading      
       col = "red",
       lwd = 2, lty=2)
abline(v = mean(loading_x1),              
       col = "blue",
       lwd = 2, lty=2)
legend(x = "topleft", # location of legend within plot area
       c("Expected", "Truth"),
       col = c("royalblue", "red"),
       lwd = c(2, 2, 2),
       cex=1)

# b
mean(loading_x2)
hist(loading_x2, breaks=40, main="Histogram of loading (X5)", xlab="X5 loading")
abline(v = v[b,1] , # true loading                
       col = "red",
       lwd = 2, lty=2)
abline(v = mean(loading_x2),                
       col = "blue",
       lwd = 2, lty=2)
legend(x = "topleft", # location of legend within plot area
       c("Expected", "Truth"),
       col = c("royalblue", "red"),
       lwd = c(2, 2, 2),
       cex= 0.8)

# c
mean(loading_x3)
hist(loading_x3, breaks=40, main="Histogram of loading (X9)", xlab="X9 loading")
abline(v = v[c,1] , # true loading                
       col = "red",
       lwd = 2, lty=2)
abline(v = mean(loading_x3),                
       col = "blue",
       lwd = 2, lty=2)
legend(x = "topleft", # location of legend within plot area
       c("Expected", "Truth"),
       col = c("royalblue", "red"),
       lwd = c(2, 2, 2),
       cex= 1)

## variance
loading_se = matrix(100,nrow=10,ncol=2)
loading_se

for (j in 1:10){
  for (k in 1:2){
    loading_vector = c()
    
    for (i in 1:1000){
      loading_vector = append(loading_vector, loading_sample_data[[i]][j,k])
    }
    
    loading_se[j,k] = sd(loading_vector)
  }
}
loading_se




