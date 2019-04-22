# Input Code for Live 3D Interaction Plot
# Be sure that you have run the function code in this R session before running this input code.
# The plot will appear in the RStudio Viewer window and a copy will be saved to your computer. 
# To publish the plot online, click the blue "Publish" button in the Viewer window. You will need to make a free RPubs account. 

# Here is an input code example. Update the arguments in quotations and copy and paste it to run in R to produce a live plot. 
plot3dInt(file.loc.name = "simulateddata.final.csv", # Data file name. Include extension. Must be .csv or .sav. If the dataset is not in the same location at the R project must also supply a file path.
	plot.name = "Life Stress x Neuroticism Predicting Depression", # Plot name and name of file that will be saved to your computer. 
	y.var = "Depression.centered", # Outcome variable. Match variable name in dataset exactly, including capitalization.
	y.printed.name = "Depression", # Name of outcome variable to appear in plot. 
	x1.var = "LifeStress.centered", # X1 predictor variable. Match variable name in dataset exactly, including capitalization.
	x1.printed.name = "Life Stress", # Name of outcome variable to appear in plot. 
	x2.var = "Neuroticism.centered", # X2 predictor variable. Match variable name in dataset exactly, including capitalization.
	x2.printed.name = "Neuroticism", # Name of outcome variable to appear in plot. 
	cov = NULL, # Covariates. Multiple covariates accepted (e.g., "c("Sex", "Age")"). If none, write NULL (with no quotations). Match variable name(s) in your dataset exactly, including capitalization.
	std.y = "center", # Centers Y variable. Other options are: "standardize" and "raw".
	std.x1 = "center", # Centers X1 variable. Other options are: "standardize" and "raw".
	std.x2 = "center", # Centers X2 variable. Other options are: "standardize" and "raw".
	regplane = "on", # If "on", regression plane of interaction effect added to plot. 
	scatter = "off", # If "on", scatter plot of observed data added to plot. 
	pred.val.ci = "off", # If "on", regression planes of lower and upper bounds of 95% confidence interval around predicted values added to plot.
	jnROSx1mod = "gradient", # If "solid", ROS for X1 as moderator shaded in one color on the regression plane. If "gradient", ROS shaded according to width of 95% confidence band. If "off", no ROS shaded. Ignored if ROS falls outside range of data.
	jnROSx2mod = "solid", # If "solid", ROS for X1 as moderator shaded in one color on the regression plane. If "gradient", ROS shaded according to width of 95% confidence band. If "off", no ROS shaded. Ignored if ROS falls outside range of data.
	crossover.x1mod = "on", # If "on", this will add a simple slope marking where slope of X2 with the outcome transitions from positive to negative or vice versa. Ignored if crossover point outside range of X1.
	crossover.x2mod = "off")  # If "on", this will add a simple slope marking where the slope of X1 with the outcome transitions from positive to negative or vice versa. Ignored if crossover point outside range of X2.

# Here are all the arguments for the function.
# plot3dInt <- function(file.loc.name, 	plot.name,
# 	y.var, # text (matching variable name in dataset)
# 	y.printed.name, # text
# 	x1.var,  # text (matching variable name in dataset)
# 	x1.printed.name,  # text
# 	x2.var,  # text (matching variable name in dataset)
# 	x2.printed. name,  # text
# 	cov = NULL,  # text. (matching variable name in dataset). can take multiple covariates. 
# 	std.y = c("raw", "center", "standardize"), 
# 	std.x1 = c("raw", "center", "standardize"),
# 	std.x2 = c("raw", "center", "standardize"), 
# 	regplane = c("on", "off"),
# 	scatter = c("on", "off"),
# 	pred.val.ci = c("on", "off"), 
# 	jnROSx1mod = c("solid", "gradient","off"),
# 	jnROSx2mod = c("solid", "gradient","off"),
# 	crossover.x1mod = c("on", "off"), 
# 	crossover.x2mod = c("on", "off")