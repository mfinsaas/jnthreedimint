library(stats)
library(MASS)

# setting correlation matrix for 2 predictors and outcome
Sigma <- matrix(c(1, 0.06, 0.11, 0.06, 1, 0.22, 0.11, 0.22, 1), 
 	ncol = 3, nrow = 3)

# generating data
mu <- rep(0,3)
rawvars <- mvrnorm(n=500, mu=mu, Sigma=Sigma)
pvars <- pnorm(rawvars)
binomvars <- qpois(1-pvars, 3, .25) # count/skewed variable
sim.data <- data.frame(rawvars[,1],binomvars[,2], rawvars[,3] )
colnames(sim.data) <- c("Neuroticism", "LifeStress", "Depression")

# center variables
sim.data$Neuroticism.c <- scale(sim.data$Neuroticism, center = TRUE, scale = FALSE)
sim.data$LifeStress.c <- scale(sim.data$LifeStress, center = TRUE, scale = FALSE)
sim.data$Depression.c <- scale(sim.data$Depression, center = TRUE, scale = FALSE)
sim.data$Neuroticism.c <- as.numeric(sim.data$Neuroticism.c)
sim.data$LifeStress.c <- as.numeric(sim.data$LifeStress.c)
sim.data$Depression.c <- as.numeric(sim.data$Depression.c)

# create dichotomized "sex" variable from life stress
mean(sim.data$LifeStress)
sim.data$Sex <- ifelse(sim.data$LifeStress <3, 0, 1)

# check correlation matrix
cor(sim.data)

# write to csv
write.csv(sim.data, "/Users/mfinsaas/klein_lab_projects/jn_threedimint/data/simulateddata.final.csv")
