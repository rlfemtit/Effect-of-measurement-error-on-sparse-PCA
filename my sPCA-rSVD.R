

library(lsa)

# Import data
load("Psparse2.RData")
Out = list(X = out$X, P =out$P, W = out$W, Z= out$Z, nzeros = out$nzeros,
           k = out$k,Propsparse= out$Propsparse )

X = Out$X # 1000 X 100 matrix / generated data : Z%*%t(P) + Error
P = Out$P # 100 X 5 matrix / 0.2 sparsity : 20% loadings per each variable are zero.
Z = Out$Z # 1000 X 5 matrix / true score

# make distinct eigen values
X_svd = svd(X)
d = X_svd$d # singular values
u = X_svd$u # 10000 x 10 matrix / each column has the unit length
v = X_svd$v # 10 x 10 matrix / each column is the eigen vector

# data with high sparsity
v[,1] = c(0.422, 0.422, 0.422, 0.422, 0, 0, 0, 0, 0.380, 0.380)
v[,2] = c(0 , 0, 0, 0, 0.489, 0.489, 0.489, 0.489, -0.147, 0.147)

d[1] = d[1] * 1.5 # distinct eigen value


X = u %*% diag(d) %*% t(v)

colnames(X) = paste0("X", 1:10)

# Procedure
# 1. initialization : loss_old=1, u and v are from rank 1 approximagion of SVD
# 2. while loop stops when : abs(loss_old - loss_new) < eps or iter > max_iter
# 3. update(hard_threshold) : u_new =      v_new = 
# 4. sparsity : low abs values of loading <- 0(in each iteration)
# 5. recursive when k>1 (k is the nr of components)


### sPCA_rSVD
#residual(X:data, u:score, v:loading)
residual <- function(X, u, v){
  X_hat = u%*%t(v)
  res = sum(rowSums((X_hat-X)^2))
  return(res)
}

sPCA_rSVD <- function(X,k, zero_loadingc, max_iter, eps){
  
  loading_vector = matrix(100,nrow=dim(X)[2], ncol=k)
  eigen_value = rep(0,k)
    
  j = 1
  while (j < k+1){
  
  # X : data, k : nr of components
  # zero_loadingc : nr of zero loadings per each PC.
  # (Assume that each PC has the same proportion of sparcity)
  
  iter = 0
  stop = 0
  
  # Initialization of u, v and loss function
  X_svd = svd(X)
  
  d = X_svd$d # singular values
  u = X_svd$u # 10000 x 10 matrix / each column has the unit length
  v = X_svd$v # 10 x 10 matrix / each column is the eigen vector
  
  u_old = u[,1]
  v_old = d[1]*v[,1]
  
  loss_old = residual(X,u_old,v_old)
  loss_vec = c(loss_old) # vector to contain loss from iterations
  
  # Update u and v
  while (stop == 0){
    iter = iter + 1
    
    v_new = t(X)%*%u_old
    ind = sort(abs(v_new), index.return = TRUE) # hard threshold
    v_new[ind$ix[1:zero_loadingc]] = 0
    
    unscaled_u = X%*%v_new
    u_new = unscaled_u / norm(unscaled_u, type="2") # u with unit vector
    
    # Calculate the loss
    loss_new = residual(X,u_new,v_new)
    loss_vec = append(loss_vec, loss_new)
    
    # Stop criteria
    if (abs(loss_old - loss_new) < eps){
      stop = 1
    }
    if (iter > max_iter){
      stop = 1
    }
    else{
      v_old = v_new 
      u_old = u_new
    }
  } # end of inner while loop
  
  X = X - u_new%*%t(v_new) # X - X_hat in order to get the remaining loadings
  
  # Assign an eigen value to a vector
  eigen_value[j] = sum(v_new^2)
  
  # Output the loading vector
  loading = v_new / norm(v_new, type="2")
  
  # Fix the sign of sample loadings using cosine similarity
  cosine_original = cosine(v_true[,j], loading)
  cosine_flipped = cosine(v_true[,j], -loading)
  
  if (cosine_original < cosine_flipped){
    loading_vector[,j] = -loading
  }else{
    loading_vector[,j] = loading
  }

  j = j + 1
  
  } # end of outer while loop
  
  # result = list("loading"=loading_vector, "loss"=loss_vec, "iter"=iter)
  result = list("eigen_value"=eigen_value, "loading"=loading_vector)
  return(result)
}
result = sPCA_rSVD(X,k=2, zero_loadingc=4, max_iter=100, eps=0.00001)
result$loading
v_true[,1:2]



### simulation experiment
loading_sample_data = list()
loading_x1 = c()
loading_x2 = c()
loading_x3 = c()
n = 5
a= 1
b= 5
c= 9

for (i in 1:1000){
  
  # sampling
  sample_data = X[sample(nrow(X),size=n,replace=FALSE),]
  
  # Obtain sample loadings from sPCA_rSVD
  result = sPCA_rSVD(sample_data,k=2, zero_loadingc=4, max_iter=100, eps=0.00001)
  sample_loading = result$loading
  
  loading_sample_data[[i]] = sample_loading
  
  loading_x1 = append(sample_loading[a,1], loading_x1)
  loading_x2 = append(sample_loading[b,1], loading_x2)
  loading_x3 = append(sample_loading[c,1], loading_x3)
}

# the expected loadings
expected_loading = Reduce("+",loading_sample_data) / length(loading_sample_data) 
rownames(expected_loading) = paste0("X", 1:10)
colnames(expected_loading) = paste0("PC",1:2)

# comparison between expected loadings vs true loading
expected_loading
v[,1:2]

# angle between expected and true loading
v1 = expected_loading[,1]
v2 = v[,1]
180 * acos( sum(v1*v2) / ( sqrt(sum(v1 * v1)) * sqrt(sum(v2 * v2)) ) ) / pi 
v1


# Calculate the standard error of loadings from the sample
loading_se = matrix(100,nrow=10,ncol=2)

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

text(v[1:10,1:2],             # sets position of labels
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
       lwd = 1.5, lty=2)
abline(v = mean(loading_x1),              
       col = "blue",
       lwd = 1.5, lty=2)
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
       lwd = 1.2, lty=2)
abline(v = mean(loading_x2),                
       col = "blue",
       lwd = 1.2, lty=2)
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

