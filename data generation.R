

set.seed(1)
library(lsa)


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

d = c(200, 10, 50, 50, 6, 5, 4, 3, 2, 1)

v = qr.Q(qr(v)) # QR decomposition to get orthogonal eigen vectors

cov_matrix = v%*%diag(d)%*%t(v)
X = mvtnorm::rmvnorm(10000, mean = mu, sigma=cov_matrix)