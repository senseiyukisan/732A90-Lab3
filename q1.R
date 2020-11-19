target_fun <- function(x) {
  return ()
}

majorizing_fun <- function(alpha) {

}

#########
### 1 ###
#########

# we set c = ?
#
c = 1
curve(eval(c) * (sqrt(2 * pi)^(-1)) * exp(-eval(c)^(2) / (2 * x)) * x^(-3/2), from=0, to=50)

# we set t_min = ?
# we split f_p(x) into:
#
# f_p1(x) with support [0,1]
# f_p2(x) with support [1,Inf]
#
# f_p1(x) is going to be uniform distribution [0,1]
# we choose our majorizing constant at biggest value of target function
#
# f_p2(x) is going to be power-law distribution [1, Inf]
# we choose our
t_min = 1
alpha = 2
curve((eval(alpha) - 1 / eval(t_min)) * (x / eval(t_min))^(-eval(alpha)), from=0, to=50, add=TRUE)

