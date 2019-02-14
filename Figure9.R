# file = file.name
# y variable = y.var
# y plot label = y.name
# x variable = x.var
# x plot label = x1.name
# w variable = w.var
# w plot label = x2.name
# plot name = plot.name


# 3D plot function --------------------------------------------------------
plot3dInt <- function(file.loc.name, file.type = c("SPSS", "csv", "xls"), y.var, y.name, x1.var, x1.name, x2.var, x2.name, 
	std.var = c("standardized", "unstandardized", "centered"), plot.name) {
	
	# packages required
	require(plotly)
	require(htmlwidgets)
	require(jtools)
	require(plyr)
	require(foreign)
	require(xlsx)
	
	# read in file ------------------------------------------------------------
	
	if (file.type == "SPSS")
		df <- read.spss(file.loc.name, to.data.frame = TRUE)
	
	if (file.type == "csv")
		df <- read.csv(file.loc.name) 
	
	if (file.type == "xlsx")
		df <- read.xlsx(file.loc.name)
	
	# scale predictors and outcome --------------------------------------------
	scale.c <- function(s) {
		s <- scale(s)
		c(s)
	}
	
	center.c <- function(s) {
		s <- scale(s, scale = FALSE)
		c(s)
	}
	
	df["y.var.new"]	<- if (std.var == "standardized") 
		scale.c(df[c(y.var)]) 
	else df[c(y.var)]
	
	df["x1.var.new"]	<- if (std.var == "standardized") 
		scale.c(df[c(x1.var)]) 
	else df[c(x1.var)]
	
	df["x2.var.new"]	<- if (std.var == "standardized") 
		scale.c(df[c(x2.var)]) 
	else df[c(x2.var)]
	
	
	df["y.var.new"]	<- if (std.var == "centered") 
		center.c(df[c(y.var)]) 
	else df[c(y.var)]
	
	df["x1.var.new"]	<- if (std.var == "centered") 
		center.c(df[c(x1.var)]) 
	else df[c(x1.var)]
	
	df["x2.var.new"]	<- if (std.var == "centered") 
		center.c(df[c(x2.var)]) 
	else df[c(x2.var)]
	
	# regression model containing significant interaction term ----------------
	regmod.int <- lm(y.var.new ~ x1.var.new * x2.var.new,
		data = df)
	
	regmod.int.summary <- summary(regmod.int)
	
	# create matrix of predicted values for regression plane -----------------------------------
	# create vectors of equally distanced values across 
	# range of actual data for x1 and x2 variables
	vecx1 <- seq(min(df["x1.var.new"], na.rm = TRUE),
		max(df["x1.var.new"], na.rm = TRUE), length.out = 100)
	
	vecx2 <- seq(min(df["x2.var.new"], na.rm = TRUE),
		max(df["x2.var.new"], na.rm = TRUE),length.out = 100)
	
	# create matrix containing predicted values based on regression model 
	# with equally spaced range of values for x and z as inputs
	pred.mat <- outer(vecx1, vecx2, function(x, z) {
		predict.lm(regmod.int, data.frame(x1.var.new = x, 
			x2.var.new = z))
	}
	)
	
	# add vector values to row and column names of matrix
	rownames(pred.mat) <- vecx1
	colnames(pred.mat) <- vecx2
	
	print(vecx1)
	print(vecx2)
	
	# identify min and max of x1 and x2 and round for later printing
	minx1 <- min(df["x1.var.new"], na.rm = TRUE)
	maxx1 <- max(df["x1.var.new"], na.rm = TRUE)
	
	minx2 <- min(df["x2.var.new"], na.rm = TRUE)
	maxx2 <- max(df["x2.var.new"], na.rm = TRUE)
	
	r.minx1 <- format(round(min(df["x1.var.new"], na.rm = TRUE), digits = 2), nsmall = 2)
	r.maxx1 <- format(round(max(df["x1.var.new"], na.rm = TRUE), digits = 2), nsmall = 2)
	
	r.minx2 <- format(round(min(df["x2.var.new"], na.rm = TRUE), digits = 2), nsmall = 2)
	r.maxx2 <- format(round(max(df["x2.var.new"], na.rm = TRUE), digits = 2), nsmall = 2)
	
	
	
	# j-n when x1 is mod ------------------------------------------------------
	# test j-n
	jn.modx1 <- johnson_neyman(regmod.int, pred = x2.var.new, modx = x1.var.new, 
		alpha = 0.05)
	print(jn.modx1)
	
	# revalue factor levels
	jn.modx1$cbands$Significance <- revalue(jn.modx1$cbands$Significance,
		c("Insignificant" = "0", "Significant" = "1"))
	
	# change to numeric
	jn.modx1$cbands$Significance <- as.numeric(levels(jn.modx1$cbands$Significance))[jn.modx1$cbands$Significance]
	
	# subset jn df to only include rows with x1 values within range of actual data
	jn.testdf.modx1 <- jn.modx1$cbands[jn.modx1$cband$x1.var.new >= minx1 &
			jn.modx1$cband$x1.var.new <= maxx1, ]
	
	# identify jn inflection points
	jn.inf.modx1 <- mapply(identical, jn.testdf.modx1[1:(nrow(jn.testdf.modx1)-1), "Significance"], 
		jn.testdf.modx1[2:nrow(jn.testdf.modx1), "Significance"])
	
	# count jn inflection points (0, 1, 2)
	jn.infcount.modx1 <- as.data.frame(table(factor(jn.inf.modx1, levels = c("FALSE", "TRUE"))))[1,"Freq"]
	print(c("count of JN inflection points when x1 is mod", jn.infcount.modx1))
	
	# row number of inflection points
	jn.inf.modx1.rnum <- which(jn.inf.modx1 == FALSE)
	
	# make j-n planes
	# no ROS
	if (jn.infcount.modx1 == 0 & jn.testdf.modx1[1, "Significance"] == 0) { 
		plot.jn.plane.x1mod <- pred.mat
		plot.jn.plane.x1mod[!is.na(plot.jn.plane.x1mod)] <- NA
		print("no ROS")
		
		# print summary
		jn.summary.x1mod <- paste0("When ", x1.name, " is acting as the moderator, the slope of ",
			"the relationship between ", x2.name, " and ", y.name, " is not ",
			"significantly different from zero at any value of ", x1.name, ".", 
			"(Actual values of ", x1.name, ": ", r.minx1, " to ", r.maxx1, ").")
		print(jn.summary.x1mod)
	} 
	# all ROS
	if (jn.infcount.modx1 == 0 & jn.testdf.modx1[1, "Significance"] == 1) { 
		plot.jn.plane.x1mod <- pred.mat
		print("all ROS")
		
		# print summary
		jn.inf.modx1.sl <- ifelse(sign(jn.testdf.modx1[1, 'Slope of x2.var.new']) > 0,
			'positive', 'negative')
		jn.summary.x1mod <- paste0("When ", x1.name, " is acting as the moderator, the slope ",
			"of the relationship between ", x2.name, " and ", y.name, " is ",
			jn.inf.modx1.sl,
			" and significantly different from zero across all values of ", x1.name,
			" (Actual values of ", x1.name, ": ", r.minx1, " to ", r.maxx1, ").")
		print(jn.summary.x1mod)
	} 
	
	# 1 band
	if (jn.infcount.modx1 == 1) {
		jn.pattern.modx1 <- paste0(c(jn.testdf.modx1[jn.inf.modx1.rnum, "Significance"], 
			jn.testdf.modx1[(jn.inf.modx1.rnum + 1), "Significance"]), collapse = "") 
		jn.compare.modx1 <- paste0(c(1, 0), collapse = "")
		print(jn.pattern.modx1)
		print("1 band")
		if (jn.pattern.modx1 == jn.compare.modx1) { # low band 
			print("low band") 
			plot.jn.plane.x1mod <- pred.mat[which(abs(vecx1 - minx1) ==
					min(abs(vecx1 - minx1))):(which(abs(vecx1 - jn.modx1$bounds[1]) ==
							min(abs(vecx1 - jn.modx1$bounds[1])))), ]
			
			# print summary
			jn.inf.modx1.sl <- ifelse(sign(jn.testdf.modx1[jn.inf.modx1.rnum, 'Slope of x2.var.new']) > 0,
				'positive', 'negative')
			jn.modx1.perc <- prop.table(table(factor(df["x1.var.new"] < jn.modx1$bounds[1])))
			jn.modx1.perc <- format(round(100*jn.modx1.perc["TRUE"], digits = 2), nsmall = 2)
			jn.summary.x1mod <- paste0("When ", x1.name, " is acting as the moderator, the slope ",
				"of the relationship between ", x2.name, " and ", y.name, " is ",
				jn.inf.modx1.sl,
				" and significantly different from zero at values of ", 
				x1.name, " between ",
				r.minx1, " and ", 
				format(round(jn.modx1$bounds[1], digits = 2), nsmall = 2), ", which includes ", jn.modx1.perc,
				"% of the sample.",
				" (Actual values of ", x1.name, ": ", r.minx1, " to ", r.maxx1, ").")
			print(jn.summary.x1mod)
		} else { # high band
			print("high band")
			plot.jn.plane.x1mod <- pred.mat[which(abs(vecx1 - jn.modx1$bounds[2]) ==
					min(abs(vecx1 - jn.modx1$bounds[2]))):which(abs(vecx1 - maxx1) ==
							min(abs(vecx1 - maxx1))), ]
			print(pred.mat[which(abs(vecx1 - jn.modx1$bounds[2]) ==
					min(abs(vecx1 - jn.modx1$bounds[2]))):which(abs(vecx1 - maxx1) ==
							min(abs(vecx1 - maxx1))), ])
			jn.plane.x1mod.diff <- nrow(pred.mat) - nrow(plot.jn.plane.x1mod) # pad j-n, match pred.mat
			jn.plane.x1mod.pad <- matrix(c(rep(0, jn.plane.x1mod.diff*100)), ncol = 100)
			jn.plane.x1mod.pad[jn.plane.x1mod.pad==0] <- NA
			plot.jn.plane.x1mod <- rbind(jn.plane.x1mod.pad, plot.jn.plane.x1mod)
			
			# print summary
			jn.inf.modx1.sl <- ifelse(sign(jn.testdf.modx1[(jn.inf.modx1.rnum + 1), 'Slope of x2.var.new']) > 0,
				'positive', 'negative')
			jn.modx1.perc <- prop.table(table(factor(df["x1.var.new"] > jn.modx1$bounds[2])))
			jn.modx1.perc <- round(100*jn.modx1.perc["TRUE"], digits = 2)
			jn.summary.x1mod <- paste0("When ", x1.name, " is acting as the moderator, the slope ",
				"of the relationship between ", x2.name, " and ", y.name, " is ",
				jn.inf.modx1.sl,
				" and significantly different from zero at values of ", 
				x1.name, " between ",
				format(round(jn.modx1$bounds[2], digits = 2), nsmall = 2), " and ", 
				r.maxx1, ", which includes ", jn.modx1.perc,
				"% of the sample.",
				" (Actual values of ", x1.name, ": ", r.minx1, " to ", r.maxx1, ").")
			print(jn.summary.x1mod)
		}
	}
	
	# 2 bands
	if (jn.infcount.modx1 == 2) {
		jn.pattern.modx1 <- paste0(c(jn.testdf.modx1[c(jn.inf.modx1.rnum[1], jn.inf.modx1.rnum[1] + 1, 
			jn.inf.modx1.rnum[2], jn.inf.modx1.rnum[2] + 1), "Significance"]), collapse = "")
		print(jn.pattern.modx1)
		jn.compare.modx1 <- paste0(c(1, 0, 0, 1), collapse = "")
		print("2 bands")
		if (jn.pattern.modx1 == jn.compare.modx1) { # low & high band
			# low band
			plot.jn.low.plane.x1mod <- pred.mat[which(abs(vecx1 - minx1) ==
					min(abs(vecx1 - minx1))):(which(abs(vecx1 - jn.modx1$bounds[1]) ==
							min(abs(vecx1 - jn.modx1$bounds[1])))), ]
			
			# high band
			plot.jn.high.plane.x1mod <- pred.mat[which(abs(vecx1 - jn.modx1$bounds[2]) ==
					min(abs(vecx1 - jn.modx1$bounds[2]))):which(abs(vecx1 - maxx1) ==
							min(abs(vecx1 - maxx1))), ]
			
			plot.jn.plane.x1mod <- pred.mat
			plot.jn.plane.x1mod[nrow(plot.jn.low.plane.x1mod):(100-nrow(plot.jn.high.plane.x1mod)),] <- NA
			
			# print summary
			jn.inf.modx1.sl.low <- ifelse(sign(jn.testdf.modx1[jn.inf.modx1.rnum, 'Slope of x2.var.new']) > 0,
				'positive', 'negative')
			jn.modx1.perc.low <- prop.table(table(factor(df["x1.var.new"] < jn.modx1$bounds[1])))
			jn.modx1.perc.low <- format(round(100*jn.modx1.perc.low["TRUE"], digits = 2), nsmall = 2)
			
			jn.inf.modx1.sl.high <- ifelse(sign(jn.testdf.modx1[(jn.inf.modx1.rnum + 1), 'Slope of x2.var.new']) > 0,
				'positive', 'negative')
			jn.modx1.perc.high <- prop.table(table(factor(df["x1.var.new"] > jn.modx1$bounds[2])))
			jn.modx1.perc.high <- format(round(100*jn.modx1.perc.high["TRUE"], digits = 2), nsmall = 2)
			
			jn.summary.x1mod <- paste0("When ", x1.name, " is acting as the moderator, the slope ",
				"of the relationship between ", x2.name, " and ", y.name, " is ",
				jn.inf.modx1.sl.low,
				" and significantly different from zero at values of ", 
				x1.name, " between ",
				r.minx1, " and ", 
				format(round(jn.modx1$bounds[1], digits = 2), nsmall = 2), ", which includes ", jn.modx1.perc.low,
				"% of the sample,", " and the slope is ", jn.inf.modx1.sl.high, " and",
				" significantly different from zero at values of ",
				x1.name, " between ", format(round(jn.modx1$bounds[2], digits = 2), nsmall = 2), " and ",
				r.maxx1, ", which includes ", jn.modx1.perc.high, "% of the sample.",
				" (Actual values of ", x1.name, ": ", r.minx1, " to ", r.maxx1, ").")
			print(jn.summary.x1mod[1])			
		} else { # center band
			print("center band")
			
			# low empty
			plot.jn.lowempty.plane.x1mod <- pred.mat[which(abs(vecx1 - minx1) ==
					min(abs(vecx1 - minx1))):(which(abs(vecx1 - jn.modx1$bounds[1]) ==
							min(abs(vecx1 - jn.modx1$bounds[1])))), ]
			
			# high empty
			plot.jn.high.plane.x1mod <- pred.mat[which(abs(vecx1 - jn.modx1$bounds[2]) ==
					min(abs(vecx1 - jn.modx1$bounds[2]))):which(abs(vecx1 - maxx1) ==
							min(abs(vecx1 - maxx1))), ]
			plot.jn.plane.x1mod <- pred.mat[nrow(plot.jn.lowempty.plane.x1mod):(100-nrow(plot.jn.high.plane.x1mod)),]
			
			jn.plane.x1mod.diff <- nrow(pred.mat) - nrow(plot.jn.lowempty.plane.x1mod) - nrow(plot.jn.high.plane.x1mod) # pad j-n, match pred.mat
			jn.plane.x1mod.pad <- matrix(c(rep(0, jn.plane.x1mod.diff*100)), ncol = 100)
			jn.plane.x1mod.pad[jn.plane.x1mod.pad==0] <- NA
			plot.jn.plane.x1mod <- rbind(jn.plane.x1mod.pad, plot.jn.plane.x1mod)
			
			# print summary
			jn.inf.modx1.sl <- ifelse(sign(jn.testdf.modx1[(jn.inf.modx1.rnum + 1), 'Slope of x2.var.new']) > 0,
				'positive', 'negative')
			jn.modx1.perc <- prop.table(table(factor(df["x1.var.new"] > jn.modx1$bounds[1] & df["x1.var.new"] < jn.modx1$bounds[2])))
			jn.modx1.perc <- format(round(100*jn.modx1.perc["TRUE"], digits = 2), nsmall = 2)
			
			jn.summary.x1mod <- paste0("When ", x1.name, " is acting as the moderator, the slope ",
				"of the relationship between ", x2.name, " and ", y.name, " is ",
				jn.inf.modx1.sl,
				" and significantly different from zero at values of ", 
				x1.name, " between ",
				format(round(jn.modx1$bounds[1], digits = 2), nsmall = 2), " and ", 
				format(round(jn.modx1$bounds[2], digits = 2), nsmall = 2), ", which includes ", jn.modx1.perc,
				"% of the sample.",
				" (Actual values of ", x1.name, ": ", r.minx1, " to ", r.maxx1, ").")
			print(jn.summary.x1mod[1])	
		}
	}
	
	# j-n when x2 is mod ------------------------------------------------------
	# test j-n
	jn.modx2 <- johnson_neyman(regmod.int, pred = x1.var.new, modx = x2.var.new, 
		alpha = 0.05)
	print(jn.modx2)
	
	# revalue factor levels
	jn.modx2$cbands$Significance <- revalue(jn.modx2$cbands$Significance,
		c("Insignificant" = "0", "Significant" = "1"))
	
	# change to numeric
	jn.modx2$cbands$Significance <- as.numeric(levels(jn.modx2$cbands$Significance))[jn.modx2$cbands$Significance]
	
	# subset jn df to only include rows with x2 values within range of actual data
	jn.testdf.modx2 <- jn.modx2$cbands[jn.modx2$cband$x2.var.new >= minx2 &
			jn.modx2$cband$x2.var.new <= maxx2, ]
	
	# identify jn inflection points
	jn.inf.modx2 <- mapply(identical, jn.testdf.modx2[1:(nrow(jn.testdf.modx2)-1), "Significance"], 
		jn.testdf.modx2[2:nrow(jn.testdf.modx2), "Significance"])
	
	# count jn inflection points (0, 1, 2)
	jn.infcount.modx2 <- as.data.frame(table(factor(jn.inf.modx2, levels = c("FALSE", "TRUE"))))[1,"Freq"]
	print(c("count of JN inflection points when x2 is mod", jn.infcount.modx2))
	
	# row number of inflection points
	jn.inf.modx2.rnum <- which(jn.inf.modx2 == FALSE)
	
	# make j-n planes
	# no ROS
	if (jn.infcount.modx2 == 0 & jn.testdf.modx2[1, "Significance"] == 0) { 
		plot.jn.plane.x2mod <- pred.mat
		plot.jn.plane.x2mod[!is.na(plot.jn.plane.x2mod)] <- NA
		print("no ROS")
		
		# print summary
		jn.summary.x2mod <- paste0("When ", x2.name, " is acting as the moderator, the slope of ",
			"the relationship between ", x1.name, " and ", y.name, " is not ",
			"significantly different from zero at any value of ", x2.name, ". ", 
			"(Actual values of ", x2.name, ": ", r.minx2, " to ", r.maxx2, ").")
		print(jn.summary.x2mod)
	} 
	# all ROS
	if (jn.infcount.modx2 == 0 & jn.testdf.modx2[1, "Significance"] == 1) { 
		plot.jn.plane.x2mod <- pred.mat
		print("all ROS")
		
		# print summary
		jn.summary.x2mod <- paste0("When ", x2.name, " is acting as the moderator, the slope ",
			"of the relationship between ", x1.name, " and ", y.name, " is ",
			jn.inf.modx2.sl,
			" and significantly different from zero across all values of ", x2.name,
			" (Actual values of ", x2.name, ": ", minx2, " to ", maxx2, ").")
		print(jn.summary.x2mod)
	} 
	
	# 1 band (significant is 2, insignificant is 1)
	if (jn.infcount.modx2 == 1) {
		jn.pattern.modx2 <- paste0(c(jn.testdf.modx2[jn.inf.modx2.rnum, "Significance"], 
			jn.testdf.modx2[(jn.inf.modx2.rnum + 1), "Significance"]), collapse = "") 
		jn.compare.modx2 <- paste0(c(1, 0), collapse = "")
		print(jn.pattern.modx2)
		print("1 band")
		if (jn.pattern.modx2 == jn.compare.modx2) { # low band 
			print("low band") 
			plot.jn.plane.x2mod <- pred.mat[, which(abs(vecx2 - minx2) ==
					min(abs(vecx2 - minx2))):(which(abs(vecx2 - jn.modx2$bounds[1]) ==
							min(abs(vecx2 - jn.modx2$bounds[1]))))]
			
			# print summary
			jn.inf.modx2.sl <- ifelse(sign(jn.testdf.modx2[jn.inf.modx2.rnum, 'Slope of x1.var.new']) > 0,
				'positive', 'negative')
			jn.modx2.perc <- prop.table(table(factor(df["x2.var.new"] < jn.modx2$bounds[1])))
			jn.modx2.perc <- format(round(100*jn.modx2.perc["TRUE"], digits = 2), nsmall = 2)
			jn.summary.x2mod <- paste0("When ", x2.name, " is acting as the moderator, the slope ",
				"of the relationship between ", x1.name, " and ", y.name, " is ",
				jn.inf.modx2.sl,
				" and significantly different from zero at values of ", 
				x2.name, " between ",
				r.minx2, " and ", 
				format(round(jn.modx2$bounds[1], digits = 2), nsmall = 2), ", which includes ", jn.modx2.perc,
				"% of the sample.",
				" (Actual values of ", x2.name, ": ", r.minx2, " to ", r.maxx2, ").")
			print(jn.summary.x2mod)
			
			
		} else { # high band
			print("high band")
			plot.jn.plane.x2mod.original <- pred.mat[, which(abs(vecx2 - jn.modx2$bounds[2]) ==
					min(abs(vecx2 - jn.modx2$bounds[2]))):which(abs(vecx2 - maxx2) ==
							min(abs(vecx2 - maxx2)))]
			jn.plane.x2mod.diff <- ncol(pred.mat) - ncol(plot.jn.plane.x2mod.original) # pad j-n, match pred.mat
			jn.plane.x2mod.pad <- matrix(c(rep(0, jn.plane.x2mod.diff*100)), nrow = 100)
			jn.plane.x2mod.pad[jn.plane.x2mod.pad==0] <- NA
			plot.jn.plane.x2mod <- cbind(jn.plane.x2mod.pad, plot.jn.plane.x2mod.original)
			
			# print summary
			jn.inf.modx2.sl <- ifelse(sign(jn.testdf.modx2[(jn.inf.modx2.rnum + 1), 'Slope of x1.var.new']) > 0,
				'positive', 'negative')
			jn.modx2.perc <- prop.table(table(factor(df["x2.var.new"] > jn.modx2$bounds[2])))
			jn.modx2.perc <- round(100*jn.modx2.perc["TRUE"], digits = 2)
			jn.summary.x2mod <- paste0("When ", x2.name, " is acting as the moderator, the slope ",
				"of the relationship between ", x1.name, " and ", y.name, " is ",
				jn.inf.modx2.sl,
				" and significantly different from zero at values of ", 
				x2.name, " between ",
				format(round(jn.modx2$bounds[2], digits = 2), nsmall = 2), " and ", 
				r.maxx2, ", which includes ", jn.modx2.perc,
				"% of the sample.",
				" (Actual values of ", x2.name, ": ", r.minx2, " to ", r.maxx2, ").")
			print(jn.summary.x2mod)
			
		}
	}
	
	# 2 bands
	if (jn.infcount.modx2 == 2) {
		jn.pattern.modx2 <- paste0(c(jn.testdf.modx2[c(jn.inf.modx2.rnum[1], jn.inf.modx2.rnum[1] + 1, 
			jn.inf.modx2.rnum[2], jn.inf.modx2.rnum[2] + 1), "Significance"]), collapse = "")
		print(jn.pattern.modx2)
		jn.compare.modx2 <- paste0(c(1, 0, 0, 1), collapse = "")
		print("2 bands")
		if (jn.pattern.modx2 == jn.compare.modx2) { # low & high band
			# low band
			plot.jn.low.plane.x2mod <- pred.mat[ , which(abs(vecx2 - minx2) ==
					min(abs(vecx2 - minx2))):(which(abs(vecx2 - jn.modx2$bounds[1]) ==
							min(abs(vecx2 - jn.modx2$bounds[1]))))]
			
			# high band
			plot.jn.high.plane.x2mod <- pred.mat[ , which(abs(vecx2 - jn.modx2$bounds[2]) ==
					min(abs(vecx2 - jn.modx2$bounds[2]))):which(abs(vecx2 - maxx2) ==
							min(abs(vecx2 - maxx2)))]
			
			plot.jn.plane.x2mod <- pred.mat
			plot.jn.plane.x2mod[ , ncol(plot.jn.low.plane.x2mod):(100-ncol(plot.jn.high.plane.x2mod))] <- NA
			
			# print summary
			jn.inf.modx2.sl.low <- ifelse(sign(jn.testdf.modx2[jn.inf.modx2.rnum, 'Slope of x1.var.new']) > 0,
				'positive', 'negative')
			jn.modx2.perc.low <- prop.table(table(factor(df["x2.var.new"] < jn.modx2$bounds[1])))
			jn.modx2.perc.low <- format(round(100*jn.modx2.perc.low["TRUE"], digits = 2), nsmall = 2)
			
			jn.inf.modx2.sl.high <- ifelse(sign(jn.testdf.modx2[(jn.inf.modx2.rnum + 1), 'Slope of x1.var.new']) > 0,
				'positive', 'negative')
			jn.modx2.perc.high <- prop.table(table(factor(df["x2.var.new"] > jn.modx2$bounds[2])))
			jn.modx2.perc.high <- format(round(100*jn.modx2.perc.high["TRUE"], digits = 2), nsmall = 2)
			
			jn.summary.x2mod <- paste0("When ", x2.name, " is acting as the moderator, the slope ",
				"of the relationship between ", x1.name, " and ", y.name, " is ",
				jn.inf.modx2.sl.low,
				" and significantly different from zero at values of ", 
				x2.name, " between ",
				r.minx2, " and ", 
				format(round(jn.modx2$bounds[1], digits = 2), nsmall = 2), ", which includes ", jn.modx2.perc.low,
				"% of the sample,", " and the slope is ", jn.inf.modx2.sl.high, " and",
				" significantly different from zero at values of ",
				x2.name, " between ", format(round(jn.modx2$bounds[2], digits = 2), nsmall = 2), " and ",
				r.maxx2, ", which includes ", jn.modx2.perc.high, "% of the sample.",
				" (Actual values of ", x2.name, ": ", r.minx2, " to ", r.maxx2, ").")
			print(jn.summary.x2mod[1])
			
		} else { # center band
			print("center band")
			
			# low empty
			plot.jn.lowempty.plane.x2mod <- pred.mat[, which(abs(vecx2 - minx2) ==
					min(abs(vecx2 - minx2))):(which(abs(vecx2 - jn.modx2$bounds[1]) ==
							min(abs(vecx2 - jn.modx2$bounds[1]))))]
			
			# high empty
			plot.jn.high.plane.x2mod <- pred.mat[, which(abs(vecx2 - jn.modx2$bounds[2]) ==
					min(abs(vecx2 - jn.modx2$bounds[2]))):which(abs(vecx2 - maxx2) ==
							min(abs(vecx2 - maxx2)))]
			plot.jn.plane.x2mod <- pred.mat[ , ncol(plot.jn.lowempty.plane.x2mod):(100-ncol(plot.jn.high.plane.x2mod))]
			
			jn.plane.x2mod.diff <- ncol(pred.mat) - ncol(plot.jn.lowempty.plane.x2mod) - ncol(plot.jn.high.plane.x2mod) # pad j-n, match pred.mat
			jn.plane.x2mod.pad <- matrix(c(rep(0, jn.plane.x2mod.diff*100)), nrow = 100)
			jn.plane.x2mod.pad[jn.plane.x2mod.pad==0] <- NA
			plot.jn.plane.x2mod <- cbind(jn.plane.x2mod.pad, plot.jn.plane.x2mod)
			
			# print summary
			jn.inf.modx2.sl <- ifelse(sign(jn.testdf.modx2[(jn.inf.modx2.rnum + 1), 'Slope of x1.var.new']) > 0,
				'positive', 'negative')
			jn.modx2.perc <- prop.table(table(factor(df["x2.var.new"] > jn.modx2$bounds[1] & df["x2.var.new"] < jn.modx2$bounds[2])))
			jn.modx2.perc <- format(round(100*jn.modx2.perc["TRUE"], digits = 2), nsmall = 2)
			
			jn.summary.x2mod <- paste0("When ", x2.name, " is acting as the moderator, the slope ",
				"of the relationship between ", x1.name, " and ", y.name, " is ",
				jn.inf.modx2.sl,
				" and significantly different from zero at values of ", 
				x2.name, " between ",
				format(round(jn.modx2$bounds[1], digits = 2), nsmall = 2), " and ", 
				format(round(jn.modx2$bounds[2], digits = 2), nsmall = 2), ", which includes ", jn.modx2.perc,
				"% of the sample.",
				" (Actual values of ", x2.name, ": ", r.minx2, " to ", r.maxx2, ").")
			print(jn.summary.x2mod[1])	
			
		}
	}
	
	color <- rep(0.5, length(df$y.var.new))
	dim(color) <- dim(df$y.var.new)
	
	color2 <- rep(1, length(t(plot.jn.plane.x1mod)))
	dim(color2) <- dim(t(plot.jn.plane.x1mod))
	
	color3 <- rep(0, length(t(plot.jn.plane.x2mod)))
	dim(color3) <- dim(t(plot.jn.plane.x2mod))
	
	# plot --------------------------------------------------------------------
	jn.plot <- plot_ly(colors = c('darkred', 'darkblue')) %>%
		# create scatterplot
		add_trace(data = df, x = ~x1.var.new,
			y = ~x2.var.new,
			z = ~y.var.new,
			type = "scatter3d",
			mode = "markers",
			marker = list(size = 0.1, opacity = 0)) %>%
		# add regression plane
		add_surface(x = ~vecx1, y = ~vecx2, z = ~t(pred.mat),
			type = "surface",
			surfacecolor=color, 
			cauto=F,
			cmax=1,
			cmin=0,
			opacity = 0.3,
			showscale = FALSE)
	# add j-n plane for x1 as mod
	jn.plot <-	add_surface(jn.plot, x = ~vecx1, y = ~vecx2, z = ~t(plot.jn.plane.x1mod),
		type = "surface",
		surfacecolor=color2,
		cauto=F,
		cmax=1,
		cmin=0,
		opacity = 0.5,
		showscale = FALSE) %>%
		# add simple slopes when x1 is mod (Figure 9)
		add_trace(x = ~vecx1[100], # jn upper bound
			y = ~vecx2, z = ~pred.mat[100,],
			type = "scatter3d",
			mode = "lines",
			line = list(width = 3, color = "red")) %>%
		add_trace(x = ~vecx1[35],
			y = ~vecx2, z = ~pred.mat[35,],
			type = "scatter3d",
			mode = "lines",
			line = list(width = 3, color = "red")) %>%
		add_trace(x = ~vecx1[53],
			y = ~vecx2, z = ~pred.mat[53,],
			type = "scatter3d",
			mode = "lines",
			line = list(width = 3, color = "red", dash = "longdash")) %>%
		add_trace(x = ~vecx1[33],
			y = ~vecx2, z = ~pred.mat[33,],
			type = "scatter3d",
			mode = "lines",
			line = list(width = 3, color = "gray", dash = "longdash")) %>%
		add_trace(x = ~vecx1[14],
			y = ~vecx2, z = ~pred.mat[14,],
			type = "scatter3d",
			mode = "lines",
			line = list(width = 3, color = "gray", dash = "longdash")) %>%
		# format
		layout(title = plot.name,
			scene = list(
				camera = list(eye = list(x = 1.6, y = -1.6, z = 1.25)),
				xaxis = list(title = x1.name, color = "blue"),
				yaxis = list(title = x2.name, color = "red"),
				zaxis = list(title = y.name, color = "black"))) %>%
		hide_legend()
	htmlwidgets::saveWidget(widget = jn.plot, "3Dplot.html", selfcontained = FALSE)
		print(vecx1)
}

plot3dInt(file.loc.name = "/Users/mfinsaas/klein_lab_projects/interactions/simulateddata.final.csv",
	file.type = "csv",
	y.var = "Depression", y.name = "Depression",
	x1.var = "LifeStress", x1.name = "Life Stress",
	x2.var = "Neuroticism", x2.name = "Neuroticism",
	std.var = "centered",
	plot.name = "Neuroticism x Life Stress Predicting Depression")
