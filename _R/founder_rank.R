##gaussian kernel
sigma = 0.1
gKernel <- function(x) 1/(sigma*sqrt(2*pi)) * exp(-(x^2)/2*sigma^2)
kernelG <- function(x, mean=0, sd=1) dnorm(x, mean = mean, sd = sd)

##kernel density estimate
dens <- function(x, bw = bw.nrd0, kernel = kernelG, n = 4096,
                from = min(x) - 3*sd, to = max(x) + 3*sd, adjust = 1,
                ...) {
  if(has.na <- any(is.na(x))) {
    x <- na.omit(x)
    if(length(x) == 0)
        stop("no finite or non-missing data!")
  }
  sd <- (if(is.numeric(bw)) bw[1] else bw(x)) * adjust
  X <- seq(from, to, len = n)
  M <- outer(X, x, kernel, sd = sd, ...)
  structure(list(x = X, y = rowMeans(M), bw = sd,
                 call = match.call(), n = length(x),
                 data.name = deparse(substitute(x)),
                 has.na = has.na), class =  "density")
}

integr <- function(x, f)
{
       
       # var is numeric
       if (!is.numeric(x))
       {
              stop('The variable of integration "x" is not numeric.')
       }

       # integrand is numeric
       if (!is.numeric(f))
       {
              stop('The integrand "f" is not numeric.')
       }

       # length(var)=length(int)
       if (length(x) != length(f))
       {
              stop('The lengths of the variable of integration and the integrand do not match.')
       }

      # get lengths of var and integrand
       n = length(x)

       # trapezoidal integration
       integral = 0.5*sum((x[2:n] - x[1:(n-1)]) * (f[2:n] + f[1:(n-1)]))

       # print definite integral
       return(integral)
}
