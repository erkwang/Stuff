#STA 250 HOMEWORK 3
#PROBLEM 1

#function for bisection
bsf = function(func, interval, tol = 1e-06, maxit = 1000, printout = FALSE, up = TRUE){
  #printout = TRUE can be used for debugging
  #check if indeed a function was given
  if (!is.function(func)) return("Error: no function to be evaluated!")
  l = min(interval)
  u = max(interval)
  iter = 0
  root = NA
  step = (u-l)/(maxit/2)
  while (iter <= maxit) {
    iter = iter + 1
    med = u/2+l/2
    if (abs(func(med))<=tol) {
      root = med;break()
    }
    else {
      if (func(med)*func(l)<0) {
        u = med
      }
      else if (func(med)*func(u)<0) {
        l = med
      }
      else {
        #when the function is not linear
        #we may need to force one side of the interval to update
        #up controls how med moves
        if (up) l = l+step else u = u-step
      }
    }
    if (printout & iter%%100 == 0) {
      print(paste("l =", l, "u =", u, "value =", func(med), "iteration:", iter))
    }
  }
  if (is.na(root)) return("No root found at maximum iteration")
  return(list(root = root, iteration = iter))
}

#===========================================================#

#function for Newton-Raphson
nrf = function(func, deriv, start, tol = 1e-06, maxit = 1000, printout = FALSE){
  #check if indeed a function was given
  if (!is.function(func)) return("Error: no function to be evaluated!")
  x = start
  iter = 0
  if (norm(func(x),"2")<=tol) return(list(root = x, iteration = iter)) 
  while (iter <= maxit) {
    iter = iter + 1
    x = x - func(x)*(deriv(x))^-1
    if (printout & iter%%100 == 0) {
      print(paste("x =", x, "value =", func(x), "iteration =", iter))
    }
    if (norm(func(x),"2")<=tol) {
      return(list(root = x, iteration = iter))
    }
  }
  return("No root found at maximum iteration")
}

#===========================================================#

#solve part c of problem 1
llh = function(lambda) {
  125*log(2+lambda) + 38*log(1-lambda) + 34*log(lambda)
}

llh.deriv = function(lambda) {
  express = D(expression(125*log(2+lambda) + 38*log(1-lambda) + 34*log(lambda)),
              "lambda")
  val = eval(express)
  val
}

llh.sec.deriv = function(lambda) {
  express = D(D(expression(125*log(2+lambda) + 38*log(1-lambda) + 34*log(lambda)),
              "lambda"), "lambda")
  val = eval(express)
  val
}

#the interval/starting value has to be specially adjusted to avoid trap
#around 0/1
solution.bs = bsf(llh.deriv, c(0.01, 0.99), maxit = 10000, 
                  tol = 1e-10)
solution.nr = nrf(llh.deriv, llh.sec.deriv, 0.1, maxit = 10000,
                  tol = 1e-10)