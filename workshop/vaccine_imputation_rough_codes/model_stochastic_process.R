# brownian motion

# generate a an empty vector to pop
w_init = 0  
w_current = c(w_init, rep(NA, times = 10000-1))

`T` = 10000

# set.seed(986747881L)
for(i in 2:`T`) {
  
  w_current[i] = w_current[i-1] + rnorm(1, 0, sqrt(i*1e-9))
  
}


# Brownian bridge

kappa = 0.5

bb_current = rep(NA, times = `T`)
for(i in 1:`T`) {
  
  bb_current[i] = (w_current[i] + (i/`T`)^kappa) - (i*w_current[`T`])/`T`
  
}



# plot(1:10000, w_current, type = "l")
plot(1:10000, bb_current, col = "darkslategrey", type = "l")

bb_current[c(1,`T`)]





