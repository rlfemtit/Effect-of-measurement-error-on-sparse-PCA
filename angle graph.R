


set.seed(1)
library(lsa)

##### Data generation (sparsity)


R = diag(10)
mu = rep(0, 10)
data = mvtnorm::rmvnorm(10000, mean = mu, sigma=R)

v = svd(data)$v
u = svd(data)$u
d = svd(data)$d

v_true = matrix(0, nrow=10, ncol=2)
v_true[,1] = c(0.422, 0.422, 0.422, 0.422, 0, 0, 0, 0, 0.380, 0.380)
v_true[,2] = c(0 , 0, 0, 0, 0.489, 0.489, 0.489, 0.489, -0.147, 0.147)

v[,1] = v_true[,1]
v[,2] = v_true[,2]

d = c(200, 100, 50, 50, 6, 5, 4, 3, 2, 1)

v = qr.Q(qr(v)) # QR decomposition to get orthogonal eigen vectors
v
cov_matrix = v%*%diag(d)%*%t(v)
X = mvtnorm::rmvnorm(10000, mean = mu, sigma=cov_matrix)




########### Inconsistency of PCA IN HDLSS case : fixed n, d-> infinite : angle converges to a quantity

# single component spike covariance model
inconsistency_PCA <- function(dimension, lambda){
  n = 10 # sample size
  angle_vec = c() # vector to contain the angles
  
  # data generation
  R = diag(dimension)
  mu = rep(0, dimension)
  data = mvtnorm::rmvnorm(1000, mean = mu, sigma=R)
  
  v_true = svd(data)$v
  u = svd(data)$u
  d = svd(data)$d
  
  d = rep(1,dimension)
  d[1] = lambda
  
  # adjusted data with desired eigen values
  data = u %*% diag(d) %*% t(v_true)
  
  # sampling
  sample_loading = svd(data[sample(nrow(data),size=n,replace=FALSE),])$v[,1:2]
  
  # Fixing the sign of sample loadings using cosine similiarity
  j=1
  while (j<2){ # only for the first component
    cosine_original = cosine(v_true[,j], sample_loading[,j])
    cosine_flipped = cosine(v_true[,j], -sample_loading[,j])
    
    if (cosine_original < cosine_flipped){
      sample_loading[,j] = -sample_loading[,j]
    }
    j=j+1
  }
  
  v1 = sample_loading[,1]
  v2 = v_true[,1]
  angle = 180 * acos( sum(v1*v2) / ( sqrt(sum(v1 * v1)) * sqrt(sum(v2 * v2)) ) ) / pi
  angle_vec = append(angle_vec, angle)
  
  return(angle_vec)
}

# line plot 
nr_variable = c(10,100,200,300,400,500,600,700,800,900,1000)
angle_vec1 = c()
angle_vec2 = c()
angle_vec3 = c()
for (i in nr_variable){
  angle_vec1 = append(angle_vec1, inconsistency_PCA(i,5))
  angle_vec2 = append(angle_vec2, inconsistency_PCA(i,10))
  angle_vec3 = append(angle_vec3, inconsistency_PCA(i,20))
}

plot(x=nr_variable, y=angle_vec1, type='o', col='red', main="Inconsistency in PCA",
     xlab='Number of variables', ylab='Angle', ylim=c(0,90))
lines(x=nr_variable, y=angle_vec2, type='o', col='green') # lambda = 10
lines(x=nr_variable, y=angle_vec3, type='o', col='blue') # lambda = 20
legend(x="topleft", legend=c("λ1=5", "λ1=10", "λ1=20"),
       col=c("red", "green", "blue"), lty=1:1, cex=1)



########## PCA : angles between the truth and the first sample eigen vector

angle_cosine <- function(n){
  loading_sample_data = list()
  angle_vec = c()

 for (i in 1:1000){
   j = 0
  
  # sampling
   sample_loading = svd(X[sample(nrow(X),size=n,replace=FALSE),])$v[,1:2]
  
  # Fixing the sign of sample loadings using cosine similiarity
  j=1
  while (j<3){
    cosine_original = cosine(v_true[,j], sample_loading[,j])
    cosine_flipped = cosine(v_true[,j], -sample_loading[,j])
    
    if (cosine_original < cosine_flipped){
      sample_loading[,j] = -sample_loading[,j]
    }
    j=j+1
  }
  
  loading_sample_data[[i]] = sample_loading
  
  v1 = sample_loading[,1]
  v2 = v_true[,1]
  angle = 180 * acos( sum(v1*v2) / ( sqrt(sum(v1 * v1)) * sqrt(sum(v2 * v2)) ) ) / pi
  angle_vec = append(angle_vec, angle)
 }
  
 expected_loading = Reduce("+",loading_sample_data) / length(loading_sample_data)
 rownames(expected_loading) = paste0("X", 1:10)
 colnames(expected_loading) = paste0("PC", 1:2)
 
 v1 = expected_loading[,1]
 v2 = v_true[,1]
 expected_angle = 180 * acos( sum(v1*v2) / ( sqrt(sum(v1 * v1)) * sqrt(sum(v2 * v2)) ) ) / pi
 result = list("angles"=angle_vec, "expected"=expected_angle)
return(result)
}
n=5
angle_cosine(n)$expected
angles_PCA = angle_cosine(n)$angles
mean(angles_PCA)
var(angles_PCA)

hist(angles_PCA, breaks=30, xlim=c(0,100))

######################### sPCA-rSVD
residual <- function(X, u, v){
  X_hat = u%*%t(v)
  res = sum(rowSums((X_hat-X)^2))
  return(res)
}

sPCA_rSVD <- function(X,k, zero_loadingc, max_iter, eps){
  
  loading_vector = matrix(100,nrow=10,ncol=k)
  
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
        u_old = u_old
      }
    } # end of inner while loop
    
    X = X - u_new%*%t(v_new) # X - X_hat in order to get the remaining loadings
    
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
  result = loading_vector
  return(result)
}
result = sPCA_rSVD(X,k=1, zero_loadingc=4, max_iter=100, eps=0.00001)
result

# Iteration to obtain angles
angle_spca <- function(n, c){
  # n : sample size
  # c : sparsity
  loading_sample_data = list()
  spca_angles = c()
  
  for (i in 1:1000){
    # sampling
    sample_data = X[sample(nrow(X),size=n,replace=FALSE),]
    
    # Obtain sample loadings from sPCA_rSVD
    result = sPCA_rSVD(sample_data,k=1, zero_loadingc=c, max_iter=100, eps=0.00001)
    loading_sample_data[[i]] = result
    
    v1 = result[,1]
    v2 = v_true[,1]
    
    sample_angle = 180 * acos( sum(v1*v2) / ( sqrt(sum(v1 * v1)) * sqrt(sum(v2 * v2)) ) ) / pi 
    spca_angles = append(spca_angles, sample_angle)
  }
  
  # expected values
  expected_loading = Reduce("+",loading_sample_data) / length(loading_sample_data)
  
  # angle
  v1 = expected_loading[,1]
  v2 = v_true[,1]
  expected_angle = 180 * acos( sum(v1*v2) / ( sqrt(sum(v1 * v1)) * sqrt(sum(v2 * v2)) ) ) / pi 
  
  results = list("expected"=expected_angle, "angles"=spca_angles)
  return(results)
}

spca_result = angle_spca(5,4) # (sample size, sparsity) 
spca_result$expected
mean(spca_result$angles)
sd(spca_result$angles)
hist(spca_result$angles, breaks=40, xlim=c(0,100))

angle_spca_vec = c()
for (i in c(5,10,15,20)){
  angle_spca_vec = append(angle_spca_vec, mean(angle_spca(i,c=4)$angles))
}




### histogram of angles
par(mfrow=c(1,1))
hist(angles_PCA, xlab="Angle", main='Histogram of angles (PCA)', breaks=40)
abline(v = mean(angles_PCA),                
       col = "red",
       lwd = 2, lty=2)
legend(x = "topright", # location of legend within plot area
       c("mean angle"),
       col = c("red"),
       lwd = c(2, 2, 2),
       cex= 0.8)

par(mfrow=c(1,2))
hist(angles_PCA, xlab="Angle", main='Histogram of angles (PCA)', breaks=40)
abline(v = mean(angles_PCA),                
       col = "red",
       lwd = 2, lty=2)
legend(x = "topright", # location of legend within plot area
       c("mean angle"),
       col = c("red"),
       lwd = c(2, 2, 2),
       cex= 0.8)



# mean angle comparison per sample size
mean_angle = c()

for (i in c(5,10,15,20)){
  mean_angle = append(mean_angle, mean(angle_cosine(i)$angles))
}

x = c(5,10,15,20)
y = mean_angle
plot(x, y, type='o', col='red', xlab='Sample Size', ylab='Mean angle', ylim=c(0,40))


lines(x=x, y=angle_spca_vec, type='o', col='blue')
legend(x="topright", legend=c("PCA", "sPCA-rSVD"),
       col=c("red", "blue"), lty=1:1, cex=1)


