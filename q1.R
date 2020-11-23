library(poweRlaw)
library(ggplot2)

#########
### 1 ###
#########

c = 2
t_min = 1.5
alpha = 1.1

plot(1, xlim = c(0,75), ylim = c(0, 1), type = "n", xlab = "", ylab = "", main = "f(x) ~ f_p(x)")
curve(eval(c) * (sqrt(2 * pi)^(-1)) * exp(-eval(c)^(2) / (2 * x)) * x^(-3/2), from=0, to=75, add=TRUE, col="black")
curve((eval(alpha) - 1 / eval(t_min)) * (x / eval(t_min))^(-eval(alpha)), from=0, to=75, add=TRUE, col="red")

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
      res = x_max_target
    }
    if (y>t_min) {
      res = power_law(y, alpha, t_min)
    }	
    res
  }, simplify = TRUE)
}

alpha = 2
c = 2
t_min = 8.6

# Find out maximum value max(f(x)) for given parameters
x_max_target = c * (sqrt(2 * pi)^(-1)) * exp(-c^(2) / (2 * (c^2/3))) * (c^2/3)^(-3/2)

# Plot both the target function and the combined majorizing function
vx = c(seq(0, t_min, t_min/10000), seq(t_min, 30, 30/10000))
plot(vx, (target_fun(vx, c)), pch=19, cex=0.4, xlab="x", ylab="density", main="Truncated normal and majorizing densities")
points(vx, (majorizing_fun(vx, alpha, t_min)), pch=19, cex=0.2, col="pink")


#########
### 3 ###
#########

Nsample=10000

rmajorizing=function(n){
  sapply(1:n,function(i){
    res=NA
    component=sample(1:2,1,prob=c(2/3,1/3))
    if(component==1) {
      res=runif(1, 0, t_min)
    }
    if(component==2){
      res=rplcon(1, t_min, alpha)
    }
    res
  })
}


fgentruncnormal=function(majorizing_constant){
  x=NA
  num_reject=0
  while (is.na(x)){
    y=rmajorizing(1)
    u=runif(1)
    if (u<=target_fun(y, 2)/(majorizing_constant*majorizing_fun(y, 2, 8.6))){
      x=y
    }
    else{
      num_reject=num_reject+1
    }
  }
  c(x,num_reject)
}



vtruncnormal_acceptreject = sapply(rep(x_max_target,Nsample),fgentruncnormal)[1,]
vtruncnormal_acceptreject = vtruncnormal_acceptreject[vtruncnormal_acceptreject <= 30]
vtruncnormal_direct = rnorm(2*Nsample)
vtruncnormal_direct = vtruncnormal_direct[vtruncnormal_direct>=0]

hist(vtruncnormal_acceptreject, col="green", breaks=100, xlab="", ylab="sample density", freq=FALSE, main="")

vtruncnormal_acceptreject_mean = mean(vtruncnormal_acceptreject)
vtruncnormal_acceptreject_var = var(vtruncnormal_acceptreject)

cat("Mean: ", vtruncnormal_acceptreject_mean)
cat("\nVar: " , vtruncnormal_acceptreject_var)

num_rejections = sum(sapply(rep(x_max_target,Nsample),fgentruncnormal)[2,])
rejection_rate = num_rejections/Nsample
cat("\nRejection rate: ", rejection_rate)

