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

target_fun <- function(x, c) {
  res = 0
  res = ((c * (sqrt(2 * pi)^(-1)) * exp(-c^(2) / (2 * x)) * x^(-3/2)))
  res[x<0] = 0
  return (res)
}

power_law <- function(x, alpha, t_min) {
  return((alpha - 1 / t_min) * (x / t_min)^(-alpha))
}


majorizing_fun <- function(x, alpha, t_min, support_border) {
  sapply(x, function(y) {
    res = NA
    if (y<0) {
      res = 0
    }
    if ((y>=0) && (y<=support_border)) {
      res = x_max
    }
    if (y>support_border) {
      res = power_law(y, alpha, t_min)
    }	
    res
  }, simplify = TRUE)
}

alpha = 1.1
c = 2
t_min = 1.5

# Find out maximum value max(f(x)) for given parameters
x_max_target = c * (sqrt(2 * pi)^(-1)) * exp(-c^(2) / (2 * (c^2/3))) * (c^2/3)^(-3/2)

# Find x value of majorizing function with same y value as x_max_target
x_values = c(seq(0, 10, by=0.001))
for (x in x_values) {
  if (x>1) {
    res = power_law(x, alpha, t_min)
    if (round(res, 5) == round(x_max, 5)) {
      x_max_power_law = x
    }
  }
}

# If we search for x value of majorzing function with same y value as x_max_target we get an x value of ~5. This is our
# border until we use the Unif distribution and afterwards the power-law but if this is our border this is also
# our T_min or not? So our initial T_min of 1.5 is not the same as the support border? Something is wrong in our logic.
support_border = x_max_power_law

vx = c(seq(0, t_min, t_min/10000), seq(t_min, 50, 50/10000))
plot(vx, (target_fun(vx, c)), pch=19, cex=0.4, xlab="x", ylab="density", main="Truncated normal and majorizing densities")
points(vx, (majorizing_fun(vx, alpha, t_min, support_border)), pch=19, cex=0.2, col="gray")

#########
### 3 ###
#########

rmajorizing<-function(n) {
  sapply(1:n,function(i) {
    res<-NA
    # What should the probabilities be? We need the integral of f(x) between 0 and T_min to get the probability
    # The probability from T_min to Inf is gonna be 1-p then. (Microsoft Teams answer to that question)
    component<-sample(1:2,1,prob=c(5/10,5/10))
    if(component==1){res<-runif(1)}
    # Here I want to sample from our power-law function as described in number 2.
    if(component==2){res<-rplcon(1, support_border, alpha)}
    res
  })
}

Nsample<-100000
num_histbreaks<-1000
hist(rmajorizing(Nsample),breaks=1000,col="black",xlab="",ylab="",main="majorizing density",freq=FALSE)

fgentruncnormal<-function(c){
  x<-NA
  num_reject<-0
  while (is.na(x)){
    y<-rmajorizing(1)
    u<-runif(1)
    if (u<=target_fun(y)/(c*majorizing_fun(y))){x<-y}
    else{num_reject<-num_reject+1}
  }
  c(x,num_reject)
}

c = 1
vtruncnormal_acceptreject = sapply(rep(c,Nsample),fgentruncnormal)[1,]
vtruncnormal_direct = rnorm(2*Nsample)
vtruncnormal_direct = vtruncnormal_direct[vtruncnormal_direct>=0]

hist(vtruncnormal_acceptreject,col="black",breaks=100,xlab="",ylab="",freq=FALSE,main="")
hist(vtruncnormal_direct,col=gray(0.8),breaks=100,xlab="",ylab="",freq=FALSE,main="",add=TRUE)
legend("topright",pch=19,cex=1.5,legend=c("acceptance/rejection algorithm","direct sampling"),col=c("black",gray(0.8)),bty="n")
