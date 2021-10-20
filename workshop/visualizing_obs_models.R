# this script visualizes over dispersion

# rho_a*eta = 0.04*0.3 == 0.012
rho_age <- 0.04
eta <- 0.30
C <- 100
m <- rho_age*eta*C 

psi <- 2

v <- m*(1 - rho_age + m*psi^2)

v2 <- m*(1 - rho_age)

set.seed(9867881)
par(mfrow = c(3, 1))

a <- rnorm(1e4, m, sqrt(v))
b <- rnorm(1e4, m, sqrt(v2))
c <- rnorm(1e4, m, sqrt(m))

hist(ifelse(a > 0, a, 0))
hist(ifelse(b > 0, b, 0))
hist(ifelse(c > 0, c, 0))





