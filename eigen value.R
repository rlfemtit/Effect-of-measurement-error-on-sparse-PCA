
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

d = c(140, 60, 1, 1, 1, 1, 1, 1, 1, 1) # PEV : around 0.7 and 0.3


v = qr.Q(qr(v)) # QR decomposition to get orthogonal eigen vectors
v
cov_matrix = v%*%diag(d)%*%t(v)
X = mvtnorm::rmvnorm(10000, mean = mu, sigma=cov_matrix)

### population eigen value
par(mfrow=c(1,1))
pca = prcomp(X)
population_eigenvalue = summary(pca)$importance[2,]
plot(population_eigenvalue, type='o',xlab='Principal Components', 
     ylab="Variance explained", col='red',xlim=c(1,10), ylim=c(0, 1),main="Comparison of eigen values")

population_eigenvalue


### sample eigenvalues from PCA
n=5
sample_data = matrix(100,nrow=1000, ncol=n)
for (i in 1:1000){
  rows = sample(1:nrow(X), n, replace = FALSE)
  sample_X = X[rows,]
  
  sample_eigenvalue = summary(prcomp(sample_X))$importance[2,]
  sample_data[i,] = sample_eigenvalue
}
sample_data

# standard error of sample eigen values
se_pca = c() 
for (i in 1:n){
  se_pca[i] = sd(sample_data[,i])
}
se_pca

colMeans(sample_data) # mean of sample eigen values
se_pca # standard error of sample eigen values

# plot for sample eigenvalues from PCA
upper_boundary <- colMeans(sample_data) + se_pca
lower_boundary <- colMeans(sample_data) - se_pca

polygon(x=c(1:n, rev(1:n)), y=c(upper_boundary, rev(lower_boundary)),col='lightblue', density=40)
lines(colMeans(sample_data), type='o', col='blue', )

### sample eigenvalues from spca-rsvd
residual <- function(X, u, v){
  X_hat = u%*%t(v)
  res = sum(rowSums((X_hat-X)^2))
  return(res)
}

sPCA_rSVD <- function(X,k, zero_loadingc, max_iter, eps){
  
  loading_vector = matrix(100,nrow=10,ncol=k)
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
result = sPCA_rSVD(X,k=1, zero_loadingc=4, max_iter=100, eps=0.00001)
result
v_true


### homogeneous measurement error 
sample_eigen_values <- function(n, c, s){
  # n : sample size
  # c : sparsity
  # s : variance of error
  
  # sample eigen values
  eigen_values = matrix(100,nrow=1000, ncol=2)
  
  # homogeneous error
  mu = rep(0, 10)
  R = s * diag(10)
  error = mvtnorm::rmvnorm(10000, mean = mu, sigma=R)
  X = X + error
  
  for (i in 1:1000){
    
    # sampling
    sample_data = X[sample(nrow(X),size=n,replace=FALSE),]
    
    # Obtain sample loadings from sPCA_rSVD
    result = sPCA_rSVD(sample_data,k=2, zero_loadingc=c, max_iter=100, eps=0.00001)
    
    # variance explained
    sum = result$eigen_value[1] + result$eigen_value[2]
    eigen_values[i,1] = result$eigen_value[1]/sum
    eigen_values[i,2] = result$eigen_value[2]/sum
  }

  return(eigen_values)  
}

# bar plot
homo_eigen_value = c()
homo_eigen_value2 = c()
homo_eigen_value3 = c()

for (i in c(10,20,30,40,50)){
  homo_eigen_value = append(homo_eigen_value, colMeans(sample_eigen_values(5, 4, i))[1])
  homo_eigen_value2 = append(homo_eigen_value2, colMeans(sample_eigen_values(10, 4, i))[1])
  homo_eigen_value3 = append(homo_eigen_value3, colMeans(sample_eigen_values(30, 4, i))[1])
}

plot(x=c(10,20,30,40,50), y=homo_eigen_value, type='o',xlab='Error Variance', ylab="PEV of the first PC"
        , main=as.expression(bquote("PEV of sPCA-rSVD per" ~σ^2)), ylim=c(0.5,0.8), col='red')
lines(x=c(10,20,30,40,50), y=homo_eigen_value2, type='o', col='blue')
lines(x=c(10,20,30,40,50), y=homo_eigen_value3, type='o', col='green')

legend(x="topright", legend=c("n = 5", "n = 10", "n = 30"),
       col=c("red", "blue", "green"), lty=1:1, cex=1, title="sample size")

### Heterogeneous measurement error
error_variance_ratio = rep(0.2,10)
error_variance_ratio[1] = 0.4
variance_x = var(X)
R = diag(error_variance_ratio * variance_x[row(variance_x)==col(variance_x)])
mu = rep(0, 10)

error = mvtnorm::rmvnorm(10000, mean = mu, sigma=R)
X = X + error


# sample eigen values from spca-rsvd
eigen_values = matrix(100,nrow=1000, ncol=2)

n=5
for (i in 1:1000){
  # sampling
  sample_data = X[sample(nrow(X),size=n,replace=FALSE),]
  
  # Obtain sample loadings from sPCA_rSVD
  result = sPCA_rSVD(sample_data,k=2, zero_loadingc=4, max_iter=100, eps=0.00001)
  
  # variance explained
  sum = result$eigen_value[1] + result$eigen_value[2]
  eigen_values[i,1] = result$eigen_value[1]/sum
  eigen_values[i,2] = result$eigen_value[2]/sum
}
eigen_values # =variance explained
colMeans(eigen_values)

# standard error of sample eigen values
se_spca = c() 
for (i in 1:2){
  se_spca[i] = sd(eigen_values[,i])
}
se_spca
se_pca

expected_eigen_values = colMeans(eigen_values)
expected_eigen_values = append(expected_eigen_values, rep(0,8)) # mean of sample eigen values
se_spca # standard error of sample eigen values

# plot for sample eigenvalues from PCA
upper_boundary <- colMeans(eigen_values) + se_spca
lower_boundary <- colMeans(eigen_values) - se_spca

polygon(x=c(1:5, rev(1:5)), y=c(upper_boundary, rev(lower_boundary)),col='lightgreen', density=40)
lines(expected_eigen_values, type='o', col='green', )
legend(x="topright", legend=c("Population", "PCA", "sPCA-rSVD"),
       col=c("red", "blue","green"), lty=1:1, cex=1)




