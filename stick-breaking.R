c <-1
G_0 <- function(n) rbeta(n, 2, 5)
n <- 7
b <- rbeta(n, 1, c)
p <- numeric(n)
p[1] <- b[1]
p[2:n] <- sapply(2:n, function(i) b[i] * prod(1 - b[1:(i-1)]))
p
y <- G_0(n)
theta <- sample(y, prob = p, replace = TRUE)
hist(theta)
quantile(theta, na.rm = FALSE,
         names = TRUE)
