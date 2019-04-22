# Figure 8. 
# Used with 1-FunctionCode_LivePlot.R code.
plot3dInt(file.loc.name = "simulateddata.final.csv",
	plot.name = "Life Stress x Neuroticism Predicting Depression", 
	y.var = "Depression", 
	y.printed.name = "Depression",
	x1.var = "Neuroticism",
	x1.printed.name = "Neuroticism",
	x2.var = "LifeStress", 
	x2.printed.name = "Life Stress",
	cov = NULL,
	std.y = "center",
	std.x1 = "center",
	std.x2 = "center",
	regplane = "on", 
	scatter = "off", 
	pred.val.ci = "off", 
	jnROSx1mod = "solid", 
	jnROSx2mod = "solid", 
	crossover.x1mod = "on",
	crossover.x2mod = "on") 


