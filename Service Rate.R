df <- read.table(pipe("pbpaste"), sep = "\t", header = TRUE)

# Assume df has two columns: df$de and df$ar,
# and that df$deviation = df$de - df$ar already exists
mu  <- mean(df$deviation)
sdv <- sd(df$deviation)

# Parameters (tunable)
p         <- 0.6   # central plateau proportion (e.g. 0.6 or 0.5)
M         <- 2     # plateau base value (theta_base at mean)
k         <- 0.05  # controls steepness of exponential tails in deviation (larger = steeper)
theta_min <- 0     # lower bound for theta (0 is fine, theta can be fractional)

# 1) Compute baseline theta_base from deviation (continuous)
z  <- (df$deviation - mu) / sdv
zL <- qnorm((1 - p) / 2)   # lower z-bound for plateau
zU <- qnorm((1 + p) / 2)   # upper z-bound for plateau

theta_base <- numeric(length(z))
for (i in seq_along(z)) {
  zi <- z[i]
  if (zi >= zL && zi <= zU) {
    # inside central plateau
    theta_base[i] <- M
  } else if (zi > zU) {
    # right tail: exponential decay as deviation increases
    theta_base[i] <- M * exp(-k * (zi - zU))
  } else {
    # left tail: exponential decay as deviation decreases
    theta_base[i] <- M * exp(-k * (zL - zi))
  }
}

# 2) Compute s = de + ar and standardise it
s      <- df$de + df$ar
s_mean <- mean(s)
s_sd   <- sd(s)
z_s    <- (s - s_mean) / s_sd

# 3) Define multiplicative factor g(s)
#    g(s) is scaled so that at the 90th percentile of s, g(s) = f_target
f_target <- 1.5
s90_z    <- qnorm(0.90)
beta     <- log(f_target) / s90_z

g_s <- exp(beta * z_s)

# 4) Combine: no rounding, no fixed global maximum
theta_cont_final <- theta_base * g_s

# 5) Enforce theta <= max(de, ar) row-wise
row_max   <- pmax(df$de, df$ar)
df$theta  <- pmin(row_max, pmax(theta_min, theta_cont_final))

# 6) Plot
plot(df$deviation, df$theta, pch = 20,
     xlab = "deviation", ylab = "theta",
     main = "theta with max(de, ar) constraint (continuous)")
lines(sort(df$deviation), df$theta[order(df$deviation)], lwd = 2)

# Copy theta back to clipboard
write.table(df$theta,
            file      = pipe("pbcopy"),
            row.names = FALSE,
            col.names = FALSE)
