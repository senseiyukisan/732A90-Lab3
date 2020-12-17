library(poweRlaw)
library(ggplot2)

#########
### 1 ###
#########

c = 1.5
t_min = 1
alpha = 1.5

plot(1, xlim = c(0,10), ylim = c(0, 0.5), type = "n", xlab = "", ylab = "", main = "f(x) ~ f_p(x)")
curve(eval(c) * (sqrt(2 * pi)^(-1)) * exp(-eval(c)^(2) / (2 * x)) * x^(-3/2), from=0, to=10, add=TRUE, col="black")
curve(ifelse(x > t_min, (eval(alpha) - 1 / eval(t_min)) * (x / eval(t_min))^(-eval(alpha)), 0), from=t_min, to=10, add=TRUE, col="red")

legend("topright", inset=.02, title="functions",
       c("target","majorizing"), horiz=TRUE, cex=0.8, col = 1:2, lty = 1)

#########
### 2 ###
#########

target_fun = function(x, c) {
  res = 0
  res = ((c * (sqrt(2 * pi)^(-1)) * exp((-c^(2)) / (2 * x)) * x^(-3/2)))
  res[x<0] = 0
  return (res)
}

power_law = function(x, alpha, t_min) {
  return(((alpha - 1) / t_min) * (x / t_min)^(-alpha))
}


majorizing_fun = function(x, alpha, t_min) {
  sapply(x, function(y) {
    res = NA
    if (y<0) {
      res = 0
    }
    if ((y>=0) && (y<=t_min)) {
      res = power_law(t_min, alpha, t_min)
    }
    if (y>t_min) {
      res = power_law(y, alpha, t_min)
    }	
    res
  }, simplify = TRUE)
}

alpha = 1.3
c = 1.5
t_min = 1

# Find out maximum value max(f(x)) for given parameters
x_max_target = c * (sqrt(2 * pi)^(-1)) * exp(-c^(2) / (2 * (c^2/3))) * (c^2/3)^(-3/2)

# Plot both the target function and the combined majorizing function
vx = c(seq(0, t_min, t_min/10000), seq(t_min, 30, 30/10000))
plot(vx, ylim=c(0, 0.5), (target_fun(vx, c)), pch=19, cex=0.4, xlab="x", ylab="density", main="Truncated normal and majorizing densities")
points(vx, (majorizing_fun(vx, alpha, t_min)), pch=19, cex=0.4, col="pink")


#########
### 3 ###
#########

Nsample=10000

rmajorizing=function(n){
  sapply(1:n,function(i){
    res=rplcon(1, t_min, alpha)
    return(res)
  })
}


fgentruncnormal=function(c_target){
  x=NA
  major_c_uniform = c_target^(-2) * 1 / sqrt(2 * pi) * 3^(3/2) * exp(-3/2)
  major_c_pow_law = 2 * c_target / sqrt(2 * pi)
  prob = integrate(target_fun, 0, t_min, c_target)$value
  num_reject=0
  num_tries=0
  while (is.na(x)){
    num_tries = num_tries + 1
    u=runif(1)
    if (u > prob) {
      y=rmajorizing(1)
      if (u<=target_fun(y, c_target)/(major_c_pow_law*majorizing_fun(y, alpha, t_min))){
        x=y
      } else{
        num_reject=num_reject+1
      }
    } else {
      y=runif(1)
      if (u<=target_fun(y, c_target)/(major_c_uniform*majorizing_fun(y, alpha, t_min))){
        x=y
      } else{
        num_reject=num_reject+1
      }
    }
  }
  c(x,num_reject,num_tries)
}

c_vals_target = c(1,5,10)
for (i in 1:length(c_vals_target)) {
  vtruncnormal_acceptreject = sapply(rep(c_vals_target[i],Nsample), fgentruncnormal)
  vtruncnormal_acceptreject_density_vals = vtruncnormal_acceptreject[1,]
  vtruncnormal_acceptreject_density_vals = vtruncnormal_acceptreject_density_vals[vtruncnormal_acceptreject_density_vals <= 500]
  
  hist(vtruncnormal_acceptreject_density_vals, col="green", breaks=100, xlab="", ylab="sample density", freq=FALSE, main="")
  
  vtruncnormal_acceptreject_mean = mean(vtruncnormal_acceptreject[1,])
  vtruncnormal_acceptreject_var = var(vtruncnormal_acceptreject[1,])
  
  cat("\nMean: ", vtruncnormal_acceptreject_mean)
  cat("\nVar: " , vtruncnormal_acceptreject_var)
  
  num_rejections = sum(vtruncnormal_acceptreject[2,])
  num_tries = sum(vtruncnormal_acceptreject[3,])
  rejection_rate = num_rejections/num_tries
  cat("\nRejection rate: ", rejection_rate)
}
