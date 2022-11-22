

# 1. p=100, 각 p에 n이 1000개인 data를 만듬 -> population line에 쓰인다.
# 2. 각 p를 10개씩 sample해서 나온 데이터에 PCA와 S-PCA를 적용한 뒤, sorted eigenvalues를 구한다
# 3. sample 수인 n을 다양하게 해서 2개의 plot(for PCA / for s-PCA)을 만든다.
library(lsa)

load("Psparse2.RData")
Out = list(X = out$X, P =out$P, W = out$W, Z= out$Z, nzeros = out$nzeros,
            k = out$k,Propsparse= out$Propsparse )

X = Out$X # 1000 X 100 matrix / generated data : Z%*%t(P) + Error
P = Out$P # 100 X 5 matrix / 0.3 sparsity : 30% loadings per each variable are zero.
Z = Out$Z # 1000 X 5 matrix / true score

colnames(X) = paste0("X", 1:10)
rownames(P) = paste0("X", 1:10)
colnames(P) = paste0("PC", 1:2)

P[,1] = P[,1] / norm(P[,1], type="2") 
P[,2] = P[,2] / norm(P[,2], type="2")

j = 0
for (first_element in P[1,]){
  j = j+1
  if (first_element < 0){
    P[,j] = -P[,j]
  }
}
population_loading = P
population_loading

# 1. Sparse loading data (with various conditions)
# 2. Unstructured data

################ when there is no noise : eigen value plot and loadign plot(first PCA and Second PCA)
# measurement error -> bias에 영향을 미치는건 당연, 어떤 좋은 insight를 끌어낼 수 있는가?(variation....?)
# 일단 실험 설계를 다 하고 multiplicative error든 correlated error든 해보자. 어떤게 sPCA에 가장 나쁘게 영향을 미칠까...

### Inconsistency of PCA when p >> n
### Eigen value of PCA when p >> n : mean line + shaed error curve


################### 1. Sparse loading data #####################
pca = prcomp(X)
population_eigenvalue = summary(pca)$importance[2,]
plot(population_eigenvalue, type='o',xlab='Principal Components', ylab="Variance explained", col='red',xlim=c(1,10), ylim=c(0, 0.6))

tuck <- TuckerCoef(loadingsTRUE, loadingsEST)  #use Tucker function
loadingsEST <- loadingsEST[,tuck$perm]           #match order of components to TRUE
for (r in 1:R) {
  corsign <- sign(sum(sign(loadingsTRUE[,r])*sign(loadingsEST[,r])))
  loadingsEST <- loadingsEST*corsign   #match sign to TRUE
}

### Sampling(n=10) from X
sample_data = matrix(,nrow=1000, ncol=10)
for (i in 1:1000){
  rows = sample(1:nrow(X), 10, replace = FALSE)
  sample_X = X[rows,]
  
  sample_eigenvalue = summary(prcomp(sample_X))$importance[2,]
  sample_data[i,] = sample_eigenvalue
}

se_pca = c() # standard error of sample eigen values
for (i in 1:10){
  se_pca[i] = sd(sample_data[,i])
}

colMeans(sample_data) # mean of sample eigen values
se_pca # standard error of sample eigen values

upper_boundary <- colMeans(sample_data) + se_pca
lower_boundary <- colMeans(sample_data) - se_pca

polygon(x=c(1:10, rev(1:10)), y=c(upper_boundary, rev(lower_boundary)),col='lightblue', density=40)
lines(colMeans(sample_data), type='o', col='blue', )
legend(7.5, 0.6, legend=c("Population", "n=10"),
       col=c("red", "blue"), lty=1:2, cex=0.8)


################### 2. Unstructured data(=no correlation between variables) #####################
# data generation
R <- diag(100)
mu <- rep(0, 100)
population_data = mvtnorm::rmvnorm(10000, mean = mu, sigma=R)
colnames(population_data) = c(paste0("x", 1:100)) # assign variable names to each column

# Apply PCA
unstructured_eigenvalue = summary(prcomp(population_data))$importance[2,]

# samples from population_data
n = 10
unstructured_sample_data = matrix(,nrow=100, ncol=n) # row : data from each iteration / column : variance explained in each PCA
for (i in 1:100){
  rows = sample(1:nrow(population_data), n, replace = FALSE)
  sample_data = population_data[rows,]
  
  unstructured_sample_eigenvalue = summary(prcomp(sample_data))$importance[2,]
  unstructured_sample_data[i,] = unstructured_sample_eigenvalue
}

un_se_pca = c() # standard error of unstructured_sample_data's eigen values
for (i in 1:n){
  un_se_pca[i] = sd(unstructured_sample_data[,i])
}

colMeans(unstructured_sample_data) # mean of sample eigen values
un_se_pca # standard error of sample eigen values

upper_boundary_un <- colMeans(unstructured_sample_data) + un_se_pca
lower_boundary_un <- colMeans(unstructured_sample_data) - un_se_pca

# plot
plot(unstructured_eigenvalue, type='o',xlab='Principal Components', ylab="Variance explained", col='red',xlim=c(1,15), ylim=c(0, 0.20))

polygon(x=c(1:n, rev(1:n)), y=c(upper_boundary_un, rev(lower_boundary_un)),col='lightblue', density=40)
lines(colMeans(unstructured_sample_data), type='o', col='blue', )

polygon(x=c(1:n, rev(1:n)), y=c(upper_boundary_un, rev(lower_boundary_un)),col='lightgreen', density=40)
lines(colMeans(unstructured_sample_data), type='o', col='green', )

legend(12, 0.20, legend=c("Population", "n=10", "n=20"), col=c("red", "blue", "green"), lty=1, cex=0.8)

############################################## Loading plot ###################################
##### unstructural data(=spherical data)
spherical_PCA = prcomp(population_data, retx=TRUE)
spherical_PCA_scores = spherical_PCA$x # PC scores
spherical_PCA_loadings = spherical_PCA$rotation # eigen vectors

# spherical_PCA_scores %*% t(spherical_PCA_loadings) # recover the original data from scores and loadings
# plot
par(mfrow=c(1,2))
plot(spherical_PCA_loadings[1:5,1:2],   # x and y data
     pch=21,              # point shape
     bg="black",          # point color
     cex=1,               # point size
     ylim=c(-0.2, 0.2),
     main="Loadings from population"      # title of plot
)
abline(a=0, b=0, h=0, v=0, lty=2)

text(spherical_PCA_loadings[1:5,1:2],             # sets position of labels
     labels=rownames(spherical_PCA_loadings[1:5,]),   # print labels
     cex=1,
     pos=3
)


# loadings from samples
loading_sample_spherical_data = list()
n = 100 # nr of samples

for (i in 1:1000){
  rows = sample(1:nrow(population_data), n, replace = FALSE)
  sample_data = population_data[sort(rows),]
  
  PCA_loadings = prcomp(sample_data)$rotation[1:5,1:2]
  
  j = 0
  for (first_element in PCA_loadings[1,]){
    j = j+1
    if (first_element < 0){
      PCA_loadings[,j] = -PCA_loadings[,j]
    }
  }
  
  loading_sample_spherical_data[[i]] = PCA_loadings
}
loading_sample_spherical_data

# the expected loadings
expected_loadings_spherical = Reduce("+",loading_sample_spherical_data) / length(loading_sample_spherical_data) 
expected_loadings_spherical

# plot from samples
plot(expected_loadings_spherical[1:5,1:2],   # x and y data
     pch=21,              # point shape
     bg="black",          # point color
     cex=1,               # point size
     ylim=c(-0.2, 0.2),
     main=paste("Loadings when n =",n)      # title of plot
)
abline(a=0, b=0, h=0, v=0, lty=2)


text(expected_loadings_spherical[1:5,1:2],       # sets position of labels
     labels=rownames(expected_loadings[1:5,]),   # print labels
     cex=1,
     pos=3
)



########## sparse loading data (when there are only the 2 true components)
# when p << n (n= 1000 / p = 10)
pca = prcomp(X)
PCAloadings = pca$rotation[1:10,1:2]

j = 0
for (first_element in PCAloadings[1,]){
  j = j+1
  if (first_element < 0){
    PCAloadings[,j] = -PCAloadings[,j]
  }
}

PCAloadings
population_loading
cosine(PCAloadings[,1], population_loading[,1])
cosine(PCAloadings[,2], population_loading[,2])

# plot
par(mfrow=c(1,3))
plot(population_loading[1:10,1:2],   # x and y data
     pch=21,              # point shape
     bg="black",          # point color
     cex=1,               # point size
     #ylim=c(-0.2, 0.2),
     main="True loadings"      # title of plot
)
abline(a=0, b=0, h=0, v=0, lty=2)

text(population_loading[1:10,1:2],             # sets position of labels
     labels=rownames(population_loading[1:10,]),   # print labels
     cex=1,
     pos=3
)

# loadings from samples
loading_sample_data = list()
n = 100 # nr of samples

for (i in 1:1000){
  rows = sample(1:nrow(X), n, replace = FALSE)
  sample_X = X[sort(rows),]
  
  PCA_loadings = prcomp(sample_X)$rotation[1:10,1:2]
  j = 0
  for (first_element in PCA_loadings[1,]){
    j = j+1
    if (first_element < 0){
      PCA_loadings[,j] = -PCA_loadings[,j]
    }
  }
  
  loading_sample_data[[i]] = PCA_loadings
}
loading_sample_data

expected_loadings = Reduce("+",loading_sample_data) / length(loading_sample_data) # the expected loadings
expected_loadings

# plot from samples
plot(expected_loadings[1:10,1:2],   # x and y data
     pch=21,              # point shape
     bg="black",          # point color
     cex=1,               # point size
     #ylim=c(-0.2, 0.2),
     main=paste("Loadings when n =",n)      # title of plot
)
abline(a=0, b=0, h=0, v=0, lty=2)


text(expected_loadings[1:10,1:2],             # sets position of labels
     labels=rownames(expected_loadings[1:10,]),   # print labels
     cex=1,
     pos=3
)


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
sample_10_se = loading_se

# sd of loadings in PC1 is greather than PC2 
expected_loadings
population_loading
sample_10_se







