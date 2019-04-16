# Figure 8. 
plot3dInt(file.loc.name = "simulateddata.final.csv", # Data file name. Include extension. Must be .csv, .xls, or .sav. If the dataset is not in the same location at the R project you will need to also supply a file path
	plot.name = "Life Stress x Neuroticism Predicting Depression", # Plot name and name of file that will be saved to your computer. 
	y.var = "Depression", # Outcome variable. This should match the variable name in your dataset exactly, including capitalization.
	y.name = "Depression", # Name of outcome variable to appear in plot. This name can be anything you want it to be. It will appear on the Y axis in the plot.
	x1.var = "LifeStress", # X1 predictor variable. This should match the variable name in your dataset exactly, including capitalization.
	x1.name = "Life Stress", # Name of outcome variable to appear in plot. This name can be anything you want it to be. It will appear on the X/Z axis in the plot.
	x2.var = "Neuroticism", # X2 predictor variable. This should match the variable name in your dataset exactly, including capitalization.
	x2.name = "Neuroticism", # Name of outcome variable to appear in plot. This name can be anything you want it to be. It will appear on the X/Z axis in the plot.
	cov = NULL, # Covariates. Multiple covariates accepted. If you don't have any, write NULL (with no quotations). These should match the variable name in your dataset exactly, including capitalization.
	std.y = "center", # Transform Y variable. Other options are: "standardize" and "raw".
	std.x1 = "center", # Transform X1 variable. Other options are: "standardize" and "raw".
	std.x2 = "center", # Transform X2 variable. Other options are: "standardize" and "raw".
	regplane = "on", # If "on", this will add the regression plane depicting the form of the interaction to the plot. 
	scatter = "off", # If "on", this will add a scatter plot of the observed data to the plot. 
	pred.val.ci = "off", # If "on", this will add semi-opaque regression planes reflecting a 95% confidence interval around the predicted values above and below the main regression plane.
	jnROSx1mod = "on", # If "on", this will shade the ROS for when X1 is the moderator on the regression plane. Only applicable if ROS falls within range of data.
	jnROSx2mod = "on", # If "on", this will shade the ROS for when X2 is the moderator on the regression plane.  Only applicable if ROS falls within range of data.
	jn.gradient.mx1 = "off", # If "on", this will add a gradient to the ROS when X1 is the moderator reflecting the width of the 95% confidence band around the slope term. It will also add a legend. Only relevant if ROS when X1 is moderator is within range of data. 
	jn.gradient.mx2 = "off", # If "on", this will add a gradient to the ROS when X2 is the moderator reflecting the width of the 95% confidence band around the slope term. It will also add a legend. Only relevant if ROS when X2 is moderator is within range of data. 
	crossover.x1mod = "on", # If "on", this will add a simple slope marking where the slope of X2 with the outcome transitions from positive to negative or vice versa. Only relevant if crossover point within range of X1.
	crossover.x2mod = "on")  # If "on", this will add a simple slope marking where the slope of X1 with the outcome transitions from positive to negative or vice versa. Only relevant if crossover point within range of X2.


