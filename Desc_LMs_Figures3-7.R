library(jtools)

sim.data <- read.csv("/Users/mfinsaas/klein_lab_projects/jn_threedimint/data/simulateddata.final.csv")

# descriptives
mean(sim.data$Neuroticism.c.c)
mean(sim.data$LifeStress.c.c)
mean(sim.data$Depression.c)

sd(sim.data$Neuroticism.c.c)
sd(sim.data$LifeStress.c.c)
sd(sim.data$Depression.c)

range(data$Neuroticism.c.c)
range(data$LifeStress.c.c)
range(data$Depression.c)

# estimate continuous x dichotomous model
mod.dichotomous <- summary(lm(Depression.c ~ Sex * Neuroticism.c, data = sim.data))
# simple slopes for appendix
sim_slopes(mod.d, pred = Neuroticism.c, modx = Sex)
sim_slopes(mod.d, pred = Sex, modx = Neuroticism.c)

# estimate continuous x continuous model
mod <- lm(Depression.c ~ Neuroticism.c * LifeStress.c, data = sim.data)
summary(mod)

# Figure 3. simple slopes plots with centered data no 95% CIs
interact_plot(mod, pred = Neuroticism.c, modx = LifeStress.c, 
	modx.labels = c("- 1 SD (-1.77)", "Mean (0.00)", "+ 1 SD (1.77)"), 
	color.class = c("darkgray","darkgray","red"),
	line.thickness = 1.3)

interact_plot(mod, pred = LifeStress.c, modx = Neuroticism.c, 
	modx.labels = c("- 1 SD (-1.02)", "Mean (0.00)", "+ 1 SD (1.02)"), 
	color.class = c("darkgray","blue","blue"),
	line.thickness = 1.3)

# Figure 5. marginal effects plots
johnson_neyman(mod, pred = Neuroticism.c, modx = LifeStress.c, sig.color = "red", 
	insig.color = "darkgray", mod.range = c(min(sim.data$LifeStress.c), max(sim.data$LifeStress.c)),
	title = NULL)
johnson_neyman(mod, pred = LifeStress.c, modx = Neuroticism.c, sig.color = "blue", 
	insig.color = "darkgray", mod.range = c(min(sim.data$Neuroticism.c), max(sim.data$Neuroticism.c)),
	title = NULL)

# Figure 6. simple slopes with J-N lower and upper bounds
interact_plot(mod, pred = Neuroticism.c, modx = LifeStress.c, modx.values = c(-1.77, 0, 0.11, 1.77, 6.05),
	modx.labels = c("SS: - 1 SD (-1.77)", "SS: Mean (0.00)", "J-N: LB (0.11)", 
		"SS: + 1 SD (1.77)", "J-N: UB (6.05)"), 
	color.class = c("darkgray","darkgray","red","red","red"),
	line.thickness = 1.3)

# confidence band table 
johnson_neyman(mod, pred = Neuroticism.c, modx = LifeStress.c, 
	mod.range = c(min(sim.data$LifeStress.c), max(sim.data$LifeStress.c)),
	title = NULL)$cband
# used output from PROCESS for figure 6

# Figure 7. regression lines in 2D space  
mod.noint <- lm(Depression.c ~ Neuroticism.c + LifeStress.c, data = sim.data)
summary(mod.noint)
plot(sim.data[c("Neuroticism.c", "Depression.c")])
abline(a = -0.00000000000000003163, b = 0.07960638525025363899, col = "red", lwd = 3)

plot(sim.data[c("LifeStress.c", "Depression.c")])
abline(a = -0.00000000000000003163, b = 0.14206996236015134727, col = "blue", lwd = 3)

