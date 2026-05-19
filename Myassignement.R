set.seed(3379)
#################### # Q1 a) #####################

# Configuration des paramètres
T <- 4                       # Nombre d'unités de temps
n <- 500                     # Nombre de pas de temps par unité de temps
M <- 1000                    # Nombre de chemins
deltaT <- T / (T * n)        # Taille de chaque pas de temps

# Matrice pour stocker les chemins de Brownian Motion
W <- matrix(0, nrow = M, ncol = T * n + 1) # M chemins, chacun avec T*n + 1 points

# Simulation des chemins
for (i in 1:M) {
  W[i, 2:(T * n + 1)] <- cumsum(sqrt(deltaT) * rnorm(T * n))
}

#################### # Q1 b) #####################

# Extraire les valeurs de W2, W3 et W4 pour chaque chemin
W2 <- W[, 2 * n]        # W_t pour t = 2
W3 <- W[, 3 * n]        # W_t pour t = 3
W4 <- W[, 4 * n]        # W_t pour t = 4

# (i) Calcul de la corrélation Corr(W2, W4)
corr_W2_W4 <- cor(W2, W4)

# (ii) Calcul de la corrélation Corr(W2, (W4 - W3)^3)
diff_W4_W3_cubed <- (W4 - W3)^3
corr_W2_diffW4W3_cubed <- cor(W2, diff_W4_W3_cubed)

# Affichage des résultats
cat("Corr(W2, W4) =", corr_W2_W4, "\n")
cat("Corr(W2, (W4 - W3)^3) =", corr_W2_diffW4W3_cubed, "\n")

#################### # Q1 c) #####################
# in the report

#################### # Q1 d) #####################
# Calcul de la variable Z = 2 * W2 - W4 pour chaque chemin
Z <- 2 * W2 - W4

# Tracé de l'histogramme de Z
hist(Z, breaks = 30, col = "lightblue", main = "Histogramme de Z = 2 * W2 - W4",
     xlab = "Valeur de Z", ylab = "Fréquence", border = "darkblue")

#################### # Q1 e) #####################
num_simulations <- 100       # Nombre de simulations Monte Carlo

# Initialiser un vecteur pour stocker la variance de Z pour chaque simulation
variances_Z <- numeric(num_simulations)

# Boucle de simulation Monte Carlo
for (sim in 1:num_simulations) {
  
  # Simulation de M chemins de Brownian Motion
  W <- matrix(0, nrow = M, ncol = T * n + 1)
  for (i in 1:M) {
    W[i, 2:(T * n + 1)] <- cumsum(sqrt(deltaT) * rnorm(T * n))
  }
  
  # Extraction de W2 et W4 pour calculer Z
  W2 <- W[, 2 * n]
  W4 <- W[, 4 * n]
  Z <- 2 * W2 - W4
  
  # Calcul de la variance de Z pour cette simulation
  variances_Z[sim] <- var(Z)
}

# Estimation finale de la variance de Z par la moyenne des variances de chaque simulation
variance_Z_MC <- mean(variances_Z)

# Affichage du résultat
cat("Estimation Monte Carlo de la variance de Z :", variance_Z_MC, "\n")

#################### # Q1 f) #####################
# in the report ( juste demander à abubakar si il faut exprimer les mouvements browniens
# en fonction d'une certaine filtration F?)

#################### # Q2 a) #####################
Delta_t <- 1/n  # Taille du pas de temps
Z_0_a <- 0.25     # Valeur initiale de Z
W_0_a <- 0        # Valeur initiale de W

# Séquence de temps
t <- Delta_t * (0:(T * n))  # Sequence de temps pour le plot

# Initialisation des vecteurs pour Z et W
Z_a <- numeric(T * n + 1)  # Vecteur de stockage pour Z
W_a <- numeric(T * n + 1)  # Vecteur de stockage pour W
Z_a[1] <- Z_0_a              # Condition initiale pour Z
W_a[1] <- W_0_a              # Condition initiale pour W

# Simulation de W et de Z
for (k in 1:(T * n)) {
  # Incrément de W avec un Brownien standard
  dW_a <- rnorm(1, mean = 0, sd = sqrt(Delta_t))
  W_a[k + 1] <- W_a[k] + dW_a  # Mise à jour de W
  
  # Mise à jour de Z en utilisant la méthode d'Euler-Maruyama
  Z_a[k + 1] <- Z_a[k] + (-1/2 * sin(W_a[k])) * Delta_t + cos(W_a[k]) * dW_a
}

#################### # Q2 b) and explication in the report #####################
# Paramètres de simulation
Z_0_b <- 0.25           # Valeur initiale de Z

# Initialisation des vecteurs pour Z et W
Z_b <- numeric(T * n + 1)     # Vecteur de stockage pour Z
Z_b[1] <- Z_0_b                 # Condition initiale pour Z
W_b <- numeric(T * n + 1)     # Vecteur de stockage pour W
W_b[1] <- 0                   # Condition initiale pour W

# Simulation plus précise de W et Z
Z_increments <- rnorm(T * n)  # Incréments iid de N(0,1) pour simuler W
for (k in 1:(T * n)) {
  # Incrément de W avec un Brownien standard
  dW_b <- sqrt(Delta_t) * Z_increments[k]
  W_b[k + 1] <- W_b[k] + dW_b
  
  # Mise à jour de Z avec une discrétisation exponentielle
  Z_b[k + 1] <- Z_b[k] * exp((-1/2 * sin(W_b[k])) * Delta_t + cos(W_b[k]) * dW_b)
}

#################### # Q2 c) #####################
# Tracé de la trajectoire de Z_a
plot(t, Z_a, type = "l", col = "blue", lwd = 2,
     xlab = "Temps", ylab = expression(Z[t]), main = "Simulation des trajectoires de Z_a_t et Z_b_t")

# Ajout de la trajectoire de Z_b sur le même graphique
lines(t, Z_b, col = "red", lwd = 2)

# Ajouter une légende pour différencier les trajectoires
legend("topright", legend = c(expression(Z[a](t)), expression(Z[b](t))),
       col = c("blue", "red"), lwd = 2, cex = 0.8)

#################### # Q2 d) #####################
# in the report

#################### Questions 3 and 4 are in F70CF_assignment_2024.R #####################








