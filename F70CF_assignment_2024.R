# F70CF 2024/25
# Assignment - Questions 3 and 4

##################################################
####  REPLACE 1234 IN THE FOLLOWING LINE WITH ####
####  THE LAST FOUR DIGITS OF YOUR STUDENT ID ####
set.seed(3379)
##################################################
rm(list = ls())    # delete everything from the workspace
graphics.off()     # close all existing graphics/plots
source("F70CF.R")  # load functions for implicit FD
####  some constants  ############################
K = 2             # strike price
T = 4
r = 0.05         # interest rate p.a.
sigma = 0.4 + rexp(1, 10)      # volatility p.a.
S_max = 8         # max share price for FD
deltaS = 0.05
M = S_max / deltaS  # number of share price intervals for FD
deltaT = 1 / 120     # length of a time step in FD method
N = T / deltaT      # number of time intervls for FD
##################################################

#### YOUR CODE SHOULD BEGIN AFTER THIS LINE   ####
#################### # Q3 a) and explaination in the report #####################
# Black-Scholes Formula for European Put Option

# Fonction pour calculer le prix d'une option de vente européenne
black_scholes_put <- function(S, K, T, r, sigma) {
  # S: current stock price
  # K: strike price
  # T: time to maturity in years
  # r: risk-free interest rate
  # sigma: volatility of the underlying asset
  
  # Calcul de d1 et d2
  d1 <- (log(S / K) + (r + sigma^2 / 2) * T) / (sigma * sqrt(T))
  d2 <- d1 - sigma * sqrt(T)
  
  # Calcul du prix de l'option de vente
  put_price <- K * exp(-r * T) * pnorm(-d2) - S * pnorm(-d1)
  
  return(put_price)
}

#################### # Q3 b) #####################
# in the report

################### Q3 c) and explaination in the report #################
# Fonction de Black-Scholes pour une option d'achat européenne (call)
black_scholes_call <- function(S0, K, T, sigma, r) {
  # S0: Prix actuel de l'action
  # K: Prix d'exercice
  # T: Maturité (4 ans)
  # sigma: Volatilité annuelle
  # r: Taux sans risque
  
  d1 <- (log(S0 / K) + (r + 0.5 * sigma^2) * T) / (sigma * sqrt(T))
  d2 <- d1 - sigma * sqrt(T)
  call_price <- S0 * pnorm(d1) - K * exp(-r * T) * pnorm(d2)
  return(call_price)
}

# Fonction de Black-Scholes pour une option de vente européenne (put)
black_scholes_put <- function(S0, K, T, sigma, r) {
  d1 <- (log(S0 / K) + (r + 0.5 * sigma^2) * T) / (sigma * sqrt(T))
  d2 <- d1 - sigma * sqrt(T)
  put_price <- K * exp(-r * T) * pnorm(-d2) - S0 * pnorm(-d1)
  return(put_price)
}

# Calcul du prix d'une Strap
strap_price <- function(S0, K, T, sigma, r) {
  # Prix de 2 options de vente et 3 options d'achat
  put_price <- black_scholes_put(S0, K, T, sigma, r)
  call_price <- black_scholes_call(S0, K, T, sigma, r)
  
  strap_price <- 2 * put_price + 3 * call_price
  return(strap_price)
}

#################### Q3 d) ####################
# Vecteur de valeurs de S0
S0_values <- seq(0, S_max, by = deltaS)

# Calcul des prix de la Strap pour chaque valeur de S0
BS_x0 <- sapply(S0_values, strap_price, K = K, T = T, sigma = sigma, r = r)

# Tracé du graphique
plot(S0_values, BS_x0, type = "l", col = "black", lwd = 2,
     xlab = expression(S[0]), ylab = expression(X[0](S[0])),
     main = "Prix de la Strap à t = 0 pour différentes valeurs de S0",
     ylim = c(0, max(BS_x0)))

#################### Q4 a) ####################
# Definition of the stock price range and initialization of the option value matrix
stock_prices <- seq(0, S_max, by = deltaS)  # Generating stock prices between 0 and S_max
option_matrix <- matrix(0, nrow = length(stock_prices), ncol = N + 1)  # Initializing the matrix for option values

# Initialization of the final payoff for the Strap at maturity time
option_matrix[, N + 1] <- 3 * pmax(stock_prices - K, 0) + 2 * pmax(K - stock_prices, 0)

# Definition of boundary conditions: stock price at 0 and at S_max
limit_low <- 2 * K * exp(-r * (T - (0:(N)) * deltaT))  # Pre-computation of the limit for S = 0 at each time step
option_matrix[1, ] <- limit_low  # Applying the lower boundary condition on the first row

limit_high <- rep(3 * S_max - K, N + 1)  # Pre-computation of the limit for S = S_max
option_matrix[length(stock_prices), ] <- limit_high  # Applying the upper boundary condition on the last row

# Applying the finite difference method to calculate the price of the European Strap
# Choosing the implicit method with the argument explicit = FALSE
option_matrix <- FD_Euro(option_matrix, r, sigma, deltaT, explicit = FALSE)  # Calculating values with the implicit method
option_matrix

#################### Q4 b) ####################
# Vecteur de valeurs de S0
S0_values <- seq(0, S_max, by = deltaS)

# Calcul des prix de la Strap pour chaque valeur de S0
BS_x0 <- sapply(S0_values, strap_price, K = K, T = T, sigma = sigma, r = r)

# Tracé du graphique initial
plot(S0_values, BS_x0, type = "l", col = "black", lwd = 2,
     xlab = expression(S[0]), ylab = expression(X[0](S[0])),
     main = "Prix de la Strap à t = 0 pour différentes valeurs de S0",
     ylim = c(0, max(BS_x0)))

# Ajout des prix pour t = 0 obtenus dans la Question 4a
lines(S0_values, option_matrix[, 1], col = "blue", lwd = 2)  # Remplacer `option_matrix[, 1]` par les prix obtenus en Question 4a


