
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

### Correlation change by measurement error
cor(X)

# Homogeneous error
mu = rep(0, 10)
R =50 * diag(10)
error = mvtnorm::rmvnorm(10000, mean = mu, sigma=R)
X_error = X + error

cor(X_error)

# Heterogeneous error

hetero_error = rep(0,10)
hetero_error[6] = 10000

R = diag(hetero_error)
mu = rep(0, 10)

error = mvtnorm::rmvnorm(10000, mean = mu, sigma=R)
X_error = X + error

svd(X_error)$v[,1]



### sPCA_rSVD
#residual(X:data, u:score, v:loading)
residual <- function(X, u, v){
  X_hat = u%*%t(v)
  res = sum(rowSums((X_hat-X)^2))
  return(res)
}

sPCA_rSVD <- function(X,k, zero_loadingc, max_iter, eps){
  
  loading_vector = matrix(100,nrow=dim(X)[2], ncol=k)
  
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
  result = list("loading"=loading_vector)
  return(result)
}
result = sPCA_rSVD(X,k=1, zero_loadingc=4, max_iter=100, eps=0.00001)
result

### homogeneous measurement error 
#error_variance_ratio = 0.4
#variance_x = var(X)
#R = diag(error_variance_ratio * variance_x[row(variance_x)==col(variance_x)])

mu = rep(0, 10)
R =100 * diag(10)
error = mvtnorm::rmvnorm(10000, mean = mu, sigma=R)
X = X + error


# simulation experiment
loading_sample_data = list()
loading_x1 = c()
loading_x2 = c()
loading_x3 = c()
n = 5
a= 1
b= 5
c= 6

for (i in 1:1000){
  
  # sampling
  sample_data = X[sample(nrow(X),size=n,replace=FALSE),]
  
  # Obtain sample loadings from sPCA_rSVD
  result = sPCA_rSVD(sample_data,k=1, zero_loadingc=4, max_iter=100, eps=0.00001)
  sample_loading = result$loading
  
  loading_sample_data[[i]] = sample_loading
  
  loading_x1 = append(sample_loading[a,1], loading_x1)
  loading_x2 = append(sample_loading[b,1], loading_x2)
  loading_x3 = append(sample_loading[c,1], loading_x3)
}

# the expected loadings
expected_loading = Reduce("+",loading_sample_data) / length(loading_sample_data)
expected_loading
v_true[,1:2]

# angle between expected and true loading
v1 = expected_loading[,1]
v2 = v_true[,1]
180 * acos( sum(v1*v2) / ( sqrt(sum(v1 * v1)) * sqrt(sum(v2 * v2)) ) ) / pi 

# proportion to correctly identify the true zero loading
length(which(loading_x2 == 0)) / 1000
length(which(loading_x3 == 0)) / 1000

var(X)
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

# get the SE of the first eigen vector when there is error in X
se_loadings_error <- function(n, c, s){
  # n : sample size
  # c : sparsity
  # s : variance of error
  
  mu = rep(0, 10)
  R = s * diag(10)
  error = mvtnorm::rmvnorm(10000, mean = mu, sigma=R)
  X = X + error
  
  # Get the first sample eigen vectors via 1000 iterations
  loading_sample_data = list()
  
  for (i in 1:1000){
    
    # sampling
    sample_data = X[sample(nrow(X),size=n,replace=FALSE),]
    
    # Obtain sample loadings from sPCA_rSVD
    result = sPCA_rSVD(sample_data,k=1, zero_loadingc=c, max_iter=100, eps=0.00001)
    sample_loading = result$loading
    
    loading_sample_data[[i]] = sample_loading
  }
  
  # Get the SE of the each loading
  loading_se = c()
  
  for (j in 1:10){
      loading_vector = c()
      
      for (i in 1:1000){
        loading_vector = append(loading_vector, loading_sample_data[[i]][j,1])
      }
      
      loading_se = append(loading_se, sd(loading_vector))
  }

return(loading_se)  
}

se_loadings_error(n=10, c=4, s=0)
var(X)

# SE sample size vs measurement error
se_loadings_error(n=10, c=4, s=0)

error_variance = c(10,20,30,40)

loading_se_error1 = c()
loading_se_error2 = c()
loading_se_error3 = c()

for (i in error_variance){
  loading_se_error1 = append(loading_se_error1, se_loadings_error(n=5, c=4, s=i)[1]) # n=5
  loading_se_error2 = append(loading_se_error2, se_loadings_error(n=10, c=4, s=i)[1]) # n=10
  loading_se_error3 = append(loading_se_error3, se_loadings_error(n=20, c=4, s=i)[1]) # n=20
}


plot(x=error_variance, y=loading_se_error1, type='o', col='red', main="Standard Error of X1 loading per Sample Size",
     xlab='Variance of Error', ylab='Standard Error', ylim=c(0.15,0.3))
lines(x=error_variance, y=loading_se_error2, type='o', col='green')
lines(x=error_variance, y=loading_se_error3, type='o', col='blue')

legend(x="topright", legend=c("n = 5", "n = 10", "n = 20"),
       col=c("red", "green", "blue"), lty=1:1, cex=1)

# bar graph
se_matrix = matrix(100,nrow=10, ncol=3)
j=0
for (i in c(0,50,100)){
  j = j+1
  se_matrix[,j] = se_loadings_error(n=5, c=4, s=i)
  
  colnames(se_matrix) = c(0, 50, 100)
  rownames(se_matrix) = paste0("X", 1:10)
}
se_matrix

par(mar=c(5, 4, 4, 8), xpd=TRUE)
barplot(t(se_matrix),
        main = as.expression(bquote("Standard Error of Loadings per" ~σ^2)),
        xlab = "Variable", ylab = "Standard Error",
        col = c("darkgrey", "darkblue", "red"),
        legend.text = colnames(se_matrix),
        args.legend = list(title=as.expression(bquote("Error variance"~σ^2)), x = "topright", inset = c(-0.2, 0)),
        xpd = TRUE,
        beside = TRUE) # Grouped bars



####### Heterogeneous measurement error
hetero_error = rep(0,10)
hetero_error[2] = 100

R = diag(hetero_error)
mu = rep(0, 10)

error = mvtnorm::rmvnorm(10000, mean = mu, sigma=R)
error
X = X + error

## change of expected loading when heterogeneous error
expected_loading_hetero <- function(n, c, s, variable){
  # n : sample size
  # c : sparsity
  # s : variance of error
  # variable : which variable has error
  
  # X = X + Error
  hetero_error = rep(0,10)
  hetero_error[variable] = s
  
  R = diag(hetero_error)
  mu = rep(0, 10)
  
  error = mvtnorm::rmvnorm(10000, mean = mu, sigma=R)
  X_error = X + error
  
  loading_sample_data = list()

  for (i in 1:1000){
    # sampling
    sample_data = X_error[sample(nrow(X_error),size=n,replace=FALSE),]
  
    # Obtain sample loadings from sPCA_rSVD
    result = sPCA_rSVD(sample_data,k=1, zero_loadingc=4, max_iter=100, eps=0.00001)
    sample_loading = result$loading
  
    loading_sample_data[[i]] = sample_loading
 }

 # the expected loadings
 expected_loading = Reduce("+",loading_sample_data) / length(loading_sample_data)
 return(expected_loading)
}
expected_loading_hetereo_matrix = matrix(100, nrow=10, ncol=3)
expected_loading_hetereo_matrix[,1] = expected_loading_hetero(n=5, c=4, s=0, variable=6)
expected_loading_hetereo_matrix[,2] = expected_loading_hetero(n=5, c=4, s=200, variable=6)
expected_loading_hetereo_matrix[,3] = expected_loading_hetero(n=5, c=4, s=400, variable=6)

var(X[,5])
expected_loading_hetereo_matrix
colnames(expected_loading_hetereo_matrix) = c("0", "200", "400")
rownames(expected_loading_hetereo_matrix) = paste0("X", 1:10)

# bar chart of loadings
par(mfrow=c(1,1))
barplot(t(expected_loading_hetereo_matrix),
        main = as.expression(bquote("Expected loadings per" ~σ^2)),
        xlab = "Variable", ylab = "Expected loading",
        col = c("darkgrey", "darkblue", "red"),
        ylim=c(-0.05,0.5),
        legend.text = colnames(expected_loading_hetereo_matrix),
        args.legend = list(title=as.expression(bquote(σ^2)), x = "topright"),
        xpd = FALSE, # legend outside of plot
        beside = TRUE) # Grouped bars

# critical point plot
hetero_error = c(200, 300)
expected_loading_hetereo_matrix = matrix(100, nrow=10, ncol=7)
j=1
for (i in hetero_error){
  expected_loading_hetereo_matrix[,j] = expected_loading_hetero(n=5, c=4, s=i, variable=1)
  j = j+1
}
expected_loading_hetereo_matrix






#### angle and classiication rate after heterogeneous measurement error added
loading_hetero <- function(n, c, s, variable){
  # n : sample size
  # c : sparsity
  # s : variance of error
  # variable : which variable has error
  
  # X = X + Error
  hetero_error = rep(0,10)
  hetero_error[variable] = s
  
  R = diag(hetero_error)
  mu = rep(0, 10)
  
  error = mvtnorm::rmvnorm(10000, mean = mu, sigma=R)
  X = X + error
  
  # Get the first sample eigen vectors via 1000 iterations
  loading_sample_data = list()
  loading_x5 = c()
  loading_x6 = c()
  loading_x7 = c()
  loading_x8 = c()
  
  for (i in 1:1000){
    
    # sampling
    sample_data = X[sample(nrow(X),size=n,replace=FALSE),]
    
    # Obtain sample loadings from sPCA_rSVD
    result = sPCA_rSVD(sample_data,k=1, zero_loadingc=c, max_iter=100, eps=0.00001)
    sample_loading = result$loading
    
    loading_sample_data[[i]] = sample_loading
    loading_x5 = append(sample_loading[5,1], loading_x5)
    loading_x6 = append(sample_loading[6,1], loading_x6)
    loading_x7 = append(sample_loading[7,1], loading_x7)
    loading_x8 = append(sample_loading[8,1], loading_x8)
  }
  
  # Get the SE of the each loading
  loading_se = c()
  
  for (j in 1:10){
    loading_vector = c()
    
    for (i in 1:1000){
      loading_vector = append(loading_vector, loading_sample_data[[i]][j,1])
    }
    
    loading_se = append(loading_se, sd(loading_vector))
  }
  
  # Get the angle between the expected and the truth
  # expected values
  expected_loading = Reduce("+",loading_sample_data) / length(loading_sample_data)
  
  # angle
  v1 = expected_loading[,1]
  v2 = v_true[,1]
  expected_angle = 180 * acos( sum(v1*v2) / ( sqrt(sum(v1 * v1)) * sqrt(sum(v2 * v2)) ) ) / pi 
  
  # classification rate
  cl_rate = matrix(100, nrow=4, ncol=1)
  
  
  cl_rate[1,1] = length(which(loading_x5 == 0)) / 1000
  cl_rate[2,1] = length(which(loading_x6 == 0)) / 1000
  cl_rate[3,1] = length(which(loading_x7 == 0)) / 1000
  cl_rate[4,1] = length(which(loading_x8 == 0)) / 1000
  
  result = list("SE"=loading_se,"expected_angle"=expected_angle, "classification"=cl_rate)
  return(result)
}




# angle graph : error variance vs angle for each of x1, x5, x9
# cl graph : error variance vs cl for each of x1, x5, x9


hetero_errors = c(50,100,200)
hetero_error1 = c()
hetero_error2 = c()
hetero_error3 = c()
hetero_error4 = c()

# x5,6,7,8
for (i in hetero_errors){
  result = loading_hetero(n=5, c=4, s=i, variable=5)$classification
  hetero_error1 = append(hetero_error1, result[1,1])
  hetero_error2 = append(hetero_error2, result[2,1])
  hetero_error3 = append(hetero_error3, result[3,1])
  hetero_error4 = append(hetero_error4, result[4,1])
}
cor(X)

# angle
plot(x=hetero_errors, y=hetero_error1, type='o', col='red', main="Angle between the expected loading and the truth",
     xlab='Variance of Error', ylab='Angle', ylim=c(0,15))
lines(x=hetero_errors, y=hetero_error2, type='o', col='green')
lines(x=hetero_errors, y=hetero_error3, type='o', col='blue')

legend(x="topleft", legend=c("x1 = 0.42", "x5 = 0", "x9 = 0.38"),
       col=c("red", "green", "blue"), lty=1:1, cex=1, title="True loading")

# classification rate
plot(x=hetero_errors, y=hetero_error1, type='o', col='red', main="Classification Rate of X5",
     xlab='Variance of Error', ylab='Classification Rate', ylim=c(0,1))
lines(x=hetero_errors, y=hetero_error2, type='o', col='green')
lines(x=hetero_errors, y=hetero_error3, type='o', col='blue')

legend(x="topright", legend=c("x1 = 0.42", "x5 = 0", "x9 = 0.38"),
       col=c("red", "green", "blue"), lty=1:1, cex=1, title="True loading")

# classification rate : zero loadings
plot(x=hetero_errors, y=hetero_error1, type='o', col='red', main="Classification Rate when error on X5",
     xlab='Variance of Error', ylab='Classification Rate', ylim=c(0,1))
lines(x=hetero_errors, y=hetero_error2, type='o', col='green')
lines(x=hetero_errors, y=hetero_error3, type='o', col='blue')
lines(x=hetero_errors, y=hetero_error4, type='o', col='black')

legend(x="topright", legend=c("x5 = 0", "x6 = 0", "x7 = 0", "x8 = 0"),
       col=c("red", "green", "blue", "black"), lty=1:1, cex=1, title="True loading")


# bar graph to show the change of estimated loadings' variance
# when error on x1 vs x5
se_matrix = matrix(100,nrow=10, ncol=3)
j=0
for (i in c(0,300,400)){
  j = j+1
  se_matrix[,j] = loading_hetero(n=5, c=4, s=i, variable=5)$SE
  
  colnames(se_matrix) = c(0, 50, 100)
  rownames(se_matrix) = paste0("X", 1:10)
}
se_matrix

par(mar=c(5, 4, 4, 8), xpd=TRUE)
barplot(t(se_matrix),
        main = as.expression(bquote("Standard Error of Loadings per" ~σ^2)),
        xlab = "Variable", ylab = "Standard Error",
        col = c("darkgrey", "darkblue", "red"),
        legend.text = colnames(se_matrix),
        args.legend = list(title=as.expression(bquote(σ^2)), x = "topright", inset = c(-0.2, 0)),
        xpd = TRUE,
        beside = TRUE) # Grouped bars





