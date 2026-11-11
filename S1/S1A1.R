library(tidyr)

# a) 
table <- data.frame(read.csv("mda.csv", sep=",", header = TRUE))
first_hist <- unlist(table[,2:21])

mu <- mean(first_hist)
mu
sigma <- sqrt(sum((first_hist - mu)^2) / length(first_hist))
sigma
#Auch mgl. aber n-1:
#sigma <- sd(first_hist)

ran <- 20
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))
hist(first_hist, 
     main = "Relatives Histogramm über alle Werte", 
     freq = FALSE, ylab = "Relative Häufigkeit", 
     xlab = "Logarithmierte Expressionswerte", 
     col = "skyblue",
     xlim = range(0,20), 
     ylim = range(0,0.25),
     breaks=ran*2)
x <- seq(0, ran, length.out = length(first_hist))
y <- dnorm(x, mean = mu, sd = sigma)
lines(x, y, col = "red", lwd = 2)
legend("topright", legend = c("Empirisch", "Normalverteilung"),
       col = c("skyblue", "red"), lwd = 2, bty = "n")

# b)
par(mfrow = c(5, 4),mar = c(2, 2.5, 1, 1))
for (i in c(2:21)) {
  mu_i <- mean(table[,i])
  sigma_i <- sqrt(sum((table[,i] - mu)^2) / length(table[,i]))
  hist(table[,i], 
       freq = FALSE, 
       main=paste("Probe:", names(table)[i]), 
       #ylab = "Rel. Häufigkeit", 
       #xlab = "Log. Expressionswerte", 
       col = "skyblue",
       xlim = range(0,20), 
       ylim = range(0,0.25),
       breaks = ran*2)
  x_i <- seq(0, ran, length.out = length(table[,i]))
  y_i <- dnorm(x_i, mean = mu_i, sd = sigma_i)
  lines(x_i, y_i, col = "red", lwd = 2)
}

# c)

par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))
plot(ecdf(first_hist), 
     main = "Verteilungsfunktion über alle Werte", 
     ylab = "F(x)", 
     xlim = range(0,20),
     xlab = "Logarithmierte Expressionswerte", 
     col = "skyblue",
     lwd=2)
x <- seq(0, ran, length.out = length(first_hist))
y <- pnorm(x, mean = mu, sd = sigma)
lines(x, y, col = "red", lwd = 1)
legend("bottomright", legend = c("Empirisch", "Normalverteilung"),
       col = c("skyblue", "red"), lwd = 2, bty = "n")

# d)
par(mfrow = c(5, 4),mar = c(2, 2.5, 1, 1))
for (i in c(2:21)) {
  mu_i <- mean(table[,i])
  sigma_i <- sqrt(sum((table[,i] - mu)^2) / length(table[,i]))
  plot(ecdf(table[,i]), 
       main=paste("Probe:", names(table)[i]), 
       #ylab = "F(x)", 
       #xlab = "Log. Expressionswerte", 
       col = "skyblue",
       xlim = range(0,20),
       lwd =2)
  x_i <- seq(0, ran, length.out = length(table[,i]))
  y_i <- pnorm(x_i, mean = mu_i, sd = sigma_i)
  lines(x_i, y_i, col = "red", lwd = 1)
}

# e)
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))
mu_list <- c()
sigma_list <- c()
probe_names <- c()
for (i in 2:21) {
  mu_i <- mean(table[,i])
  sigma_i <- sqrt(sum((table[,i] - mu_i)^2) / length(table[,i]))  # MLE für σ
  mu_list <- c(mu_list, mu_i)
  sigma_list <- c(sigma_list, sigma_i)
  probe_names <- c(probe_names, names(table)[i])
}
params_df <- data.frame(probe = probe_names, mu = mu_list, sigma = sigma_list)
groups <- rep(c("DMSO", "BA", "MA", "CAI3", "MSBAS"), each = 4)
params_df$treatment <- groups

boxplot(mu ~ treatment,
        data = params_df,
        main = "Verteilung der MLE-Schätzer für μ je Behandlung",
        xlab = "Behandlung",
        ylab = "μ",
        col = "lightblue")

boxplot(sigma ~ treatment,
        data = params_df,
        main = "Verteilung der MLE-Schätzer für σ je Behandlung",
        xlab = "Behandlung",
        ylab = "σ",
        col = "lightblue")

# f)
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))
expression_long <- data.frame(
  value = unlist(table[, 2:21]),
  probe = rep(names(table)[2:21], each = nrow(table))
)
expression_long$Behandlung <- rep(c("DMSO", "BA", "MA", "CAI3", "MSBAS"), each = 4 * nrow(table))
farben <- rep(c("skyblue", "lightgreen", "grey", "lightyellow", "salmon"), each = 4)
boxplot(value ~ probe,
        data = expression_long,
        col = farben,
        las = 2,
        main = "Boxplots gruppiert nach Behandlung",
        ylab = "Logarithmierte Expressionswerte",
        xlab ="",
        ylim = c(0, 20))

