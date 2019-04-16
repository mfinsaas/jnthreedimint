library(MASS)
library(stats)

# panel a
#sample size
n <- 100
# regression coefficients
beta0 <- 0.8
beta1 <- 0.2
# generate covariate values
x <- runif(n=n, min=0, max=15)
# compute mu's
mu <- exp(beta0 + beta1 * x)
#generate Y-values
y <- rpois(n=n, lambda=mu)

z <- rnorm(50, mean=0, n=n)
hist(z)

summary(lm(y ~ x*z))

# data set
data <- data.frame(y=y, x=x)

# histogram
hist(data$y,
	xlab = "Scores on Depression Measure (Clinical Cut-off = 30)",
	main = "Histogram of Generated Count Data \nIllustrating the Impact of Skew \non Standard Deviations")

# add -1 SD line
abline(v = mean(data$y) - sd(data$y),
 lwd = 2,
 lty = 3)

# add mean line
abline(v = mean(data$y),
 lty = 1,
 lwd = 2)

# add +1 SD line
abline(v = mean(data$y) + sd(data$y),
 lty = 2,
 lwd = 2)

a <- mean(data$y) - sd(data$y)
a <- format(round(a, digits = 2), nsmall = 2)
b <- mean(data$y)
c <- mean(data$y) + sd(data$y)
c <- format(round(c, digits = 2), nsmall = 2)

legend(x = "topright", # location of legend within plot area
 c(paste0("-1 SD (", a, ")"), 
 	paste0("Mean (", b, ")"), paste0("+1 SD (", c, ")")),
 lty = c(3, 1, 2),
 lwd = c(2, 2, 2))

# panel b 
x <- rnbinom(100, mu = 10, size = 1)
x <- x/10
x <- round(x, digits = 0)
hist(x, main = "Histogram of Generated Data \nIllustrating -1 SD Falling Outside \nRange of Observed Data",
	xlab = "Household Income Category", xlim = c(-.3, 6))

abline(v = mean(x) - sd(x),
 lty = 3,
 lwd = 2)

abline(v = mean(x),
 lty = 1,
 lwd = 2)

abline(v = mean(x) + sd(x),
 lty = 2,
 lwd = 2)

a <- mean(x) - sd(x)
a <- format(round(a, digits = 2), nsmall = 2)
b <- mean(x)
b <- format(round(b, digits = 2), nsmall = 2)
c <- mean(x) + sd(x)
c <- format(round(c, digits = 2), nsmall = 2)

legend(x = "topright", # location of legend within plot area
 c(paste0("-1 SD (", a, ")"), 
 	paste0("Mean (", b, ")"), paste0("+1 SD (", c, ")")),
 lty = c(3, 1, 2),
 lwd = c(2, 2, 2))

