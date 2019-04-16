# 3D plot function --------------------------------------------------------
figure9 <- function(file.loc.name, y.var, y.name, x1.var, x1.name, x2.var, x2.name, 
	cov = NULL, cov.name = NULL, 
	std.y = c("raw", "center", "standardize"), 
	std.x1 = c("raw", "center", "standardize"),
	std.x2 = c("raw", "center", "standardize"), 
	jn.gradient.mx1 = c("on", "off"),
	jn.gradient.mx2 = c("on", "off"),	
	pred.val.ci = c("on", "off"), 
	crossover.x1mod = c("on", "off"), 
	crossover.x2mod = c("on", "off"),
	scatter = c("on", "off"),
	regplane = c("on", "off"),
	jnROSx1mod = c("on", "off"),
	jnROSx2mod = c("on", "off"),
	plot.name) {
	
	# packages required
	require(plotly)
	require(htmlwidgets)
	require(jtools)
	require(plyr)
	require(foreign)
	require(xlsx)
	require(psych)
	
	# helper functions --------------------------------------------------------
	# returns whether ROS falues between J-N values
	johnson_neyman_inside <- source("/Users/mfinsaas/klein_lab_projects/jn_threedimint/src/johnson_neyman_inside.R")	
	
	# adapts johnson-neyman from jtools by using vec as the moderator values and 
	# testing slopes/cbs across these values
	johnson_neyman_veccb2 <- source("/Users/mfinsaas/klein_lab_projects/jn_threedimint/src/johnson_neyman_veccb2.R")	
	
	# read in file ------------------------------------------------------------
	df <- if(tools::file_ext(file.loc.name) == "sav") {
		foreign::read.spss(file.loc.name, to.data.frame = TRUE)
	} else {
		if(tools::file_ext(file.loc.name) == "xlsx") {
			read.xlsx(file.loc.name)
		} else {
			if(tools::file_ext(file.loc.name) == "xls") {
				read.xls(file.loc.nameh)
			} else {
				if(tools::file_ext(file.loc.name) == "csv") {
					read.csv(file.loc.name, header = TRUE)
				}
			}
		}
	}
	
	# if(file.type == "SPSS") {df <- read.spss(file.loc.name, to.data.frame = TRUE)}
	# if(file.type == "csv") {df <- read.csv(file.loc.name)} 
	# if(file.type == "xlsx") {df <- read.xlsx(file.loc.name)}
	
	df <- df[,c(y.var, x1.var, x2.var, cov)] 
	
	if (!is.null(cov)) {
		covn <- paste0("c",seq_along(cov))
		colnames(df) <- c("y", "x1", "x2", covn)
		df[complete.cases(df), ] # keep only cases with data on all variables
	} else {
		colnames(df) <- c("y", "x1", "x2")
		df[complete.cases(df), ] # keep only cases with data on all variables
	}
	
	# descriptives ------------------------------------------------------------
	ds <- describe(df)[c("n","mean","sd","median","min","max","skew","kurtosis")]
	print(ds)
	
	# scale predictors and outcome --------------------------------------------
	if (std.x1 == "standardize") {
		df$x1 <- c(scale(df$x1))
	}  else {
		if (std.x1 == "raw") {
			df$x1 <- df$x1
		} else {
			if (std.x1 == "center") {
				df$x1 <- c(scale(df$x1, scale = FALSE))
			}
		}
	}
	if (std.x2 == "standardize") {
		df$x2 <- c(scale(df$x2))
	}  else {
		if (std.x2 == "raw") {
			df$x2 <- df$x2
		} else {
			if (std.x2 == "center") {
				df$x2 <- c(scale(df$x2, scale = FALSE))
			}
		}
	}
	if (std.y == "standardize") {
		df$y <- c(scale(df$y))
	}  else {
		if (std.y == "raw") {
			df$y <- df$y
		} else {
			if (std.y == "center") {
				df$y <- c(scale(df$y, scale = FALSE))
			}
		}
	}
	
	if (!is.null(cov)) {
		df[,4:ncol(df)] <- c(scale(df[,4:ncol(df)] , scale = FALSE)) # covariates always centered
	}
	
	# fit regression model ----------------
	mod <- lm(y ~ . + x1 * x2, data = df)
	mod.summary <- summary(mod)
	print(mod.summary)
	
	if(mod.summary$coefficients[nrow(mod.summary$coefficients), 4] < 0.05) {
		intsig <- 1 
	} else {
		intsig <- 0
	}
	
	# identify min and max of x1 and x2 and round for later printing
	minx1 <- min(df$x1, na.rm = TRUE)
	maxx1 <- max(df$x1, na.rm = TRUE)
	minx2 <- min(df$x2, na.rm = TRUE)
	maxx2 <- max(df$x2, na.rm = TRUE)
	
	r.minx1 <- format(round(minx1, digits = 2), nsmall = 2)
	r.maxx1 <- format(round(maxx1, digits = 2), nsmall = 2)
	r.minx2 <- format(round(minx2, digits = 2), nsmall = 2)
	r.maxx2 <- format(round(maxx2, digits = 2), nsmall = 2)
	
	
	
	# johnson neyman analysis -------------------------------------------------
	# when x1 is mod ------------------------------------------------------
	jn.mx1 <- johnson_neyman(mod, pred = x2, modx = x1, 
		mod.range = c(minx1, maxx1), alpha = 0.05, plot = FALSE)
	print(jn.mx1)
	
	ss.mx1 <- sim_slopes(mod, pred = x2, modx = x1)
	print(ss.mx1)
	
	# save jn values if within mod range
	if(jn.mx1$bounds[1] >= minx1 & jn.mx1$bounds[1] <= maxx1) {
		jn.mx1.1 <- as.numeric(jn.mx1$bounds[1])
	} else {jn.mx1.1 <- NA}
	
	if(jn.mx1$bounds[2] >= minx1 & jn.mx1$bounds[2] <= maxx1) {
		jn.mx1.2 <- as.numeric(jn.mx1$bounds[2])
	} else {jn.mx1.2 <- NA}
	
	# identify whether ROS is inside j-n values
	jn.mx1.inside <- johnson_neyman_inside(mod, pred = x2, modx = x1, 
		mod.range = c(minx1, maxx1), alpha = 0.05)
	
	# when x2 is mod ------------------------------------------------------
	jn.mx2 <- johnson_neyman(mod, pred = x1, modx = x2, 
		mod.range = c(minx2, maxx2), alpha = 0.05, plot = FALSE)
	print(jn.mx2)
	
	# save jn values if within mod range
	if(jn.mx2$bounds[1] >= minx2 & jn.mx2$bounds[1] <= maxx2) {
		jn.mx2.1 <- as.numeric(jn.mx2$bounds[1])
	} else {jn.mx2.1 <- NA}
	
	if(jn.mx2$bounds[2] >= minx2 & jn.mx2$bounds[2] <= maxx2) { 
		jn.mx2.2 <- as.numeric(jn.mx2$bounds[2])
	} else {jn.mx2.2 <- NA}
	
	# identify whether ROS is inside j-n values
	jn.mx2.inside <- johnson_neyman_inside(mod, pred = x1, modx = x2, 
		mod.range = c(minx2, maxx2), alpha = 0.05)
	
	# vectors of x1 and x2 -----------------------------------
	vecx1 <- seq(minx1, maxx1, length.out = 100)
	vecx2 <- seq(minx2, maxx2, length.out = 100)
	
	# add in j-n values 
	vecx1 <- c(vecx1, jn.mx1.1, jn.mx1.2)
	vecx1 <- sort(vecx1)
	vecx2 <- c(vecx2, jn.mx2.1, jn.mx2.2)
	vecx2 <- sort(vecx2)
	
	# predicted values matrix   -----------
	pm <- outer(vecx1, vecx2, function(pmx1, pmx2) {
		if (is.null(cov)) {
			newdatapm <- data.frame(pmx1, pmx2)
			names(newdatapm) <- c("x1", "x2")
		} else {
			if (length(cov) == 1) {
				x <- mean(df[,4])
				pmcov.l <- list(x)
			} else {
				pmcov.l <- sapply(c(df[,4:ncol(df)]), mean, simplify = FALSE)
			}
			newdatapm <- data.frame(pmx1, pmx2, pmcov.l)
			names(newdatapm) <- c("x1", "x2", covn)
		}
		predict.lm(mod, newdatapm)
	})
	rownames(pm) <- vecx1
	colnames(pm) <- vecx2
	
	# confidence band matrices  -------------------------------------
	# lower ci of predicted values matrix -------------------------------------------------
	pmcil <- outer(vecx1, vecx2, function(pmx1, pmx2) {
		if (is.null(cov)) {
			newdatapm <- data.frame(pmx1, pmx2)
			names(newdatapm) <- c("x1", "x2")
		} else {
			if (length(cov) == 1) {
				x <- mean(df[,4])
				pmcov.l <- list(x) 
			} else {
				pmcov.l <- sapply(c(df[,4:ncol(df)]), mean, simplify = FALSE)
			}
			newdatapm <- data.frame(pmx1, pmx2, pmcov.l)
			names(newdatapm) <- c("x1", "x2", covn)
		}
		predict.lm(mod, newdatapm, se.fit = TRUE, interval = "confidence",
			level = 0.95)$fit[,"lwr"]
	})
	rownames(pmcil) <- vecx1
	colnames(pmcil) <- vecx2
	
	# upper ci of predicted values matrix -------------------------------------------------
	pmciu <- outer(vecx1, vecx2, function(pmx1, pmx2) {
		if (is.null(cov)) {
			newdatapm <- data.frame(pmx1, pmx2)
			names(newdatapm) <- c("x1", "x2")
		} else {
			if (length(cov) == 1) {
				x <- mean(df[,4])
				pmcov.l <- list(x) 
			} else {
				pmcov.l <- sapply(c(df[,4:ncol(df)]), mean, simplify = FALSE)
			}
			newdatapm <- data.frame(pmx1, pmx2, pmcov.l)
			names(newdatapm) <- c("x1", "x2", covn)
		}
		predict.lm(mod, newdatapm, se.fit = TRUE, interval = "confidence",
			level = 0.95)$fit[,"upr"]
	})
	rownames(pmciu) <- vecx1
	colnames(pmciu) <- vecx2
	
	# ROS matrices --------------------------------------------
	# for x1 as mod
	pm.rn <- rownames(pm)
	pm.rn <- as.numeric(pm.rn)
	pm.nrow <- nrow(pm)
	pm.ncol <- ncol(pm)
	jn.mx1.1rn <- match(jn.mx1.1, rownames(pm))
	jn.mx1.2rn <- match(jn.mx1.2, rownames(pm))
	
	
	# if inside j-n values
	pm.mx1.ROS1i <- if(jn.mx1.inside == TRUE) {
		pm[jn.mx1.1rn:jn.mx1.2rn,]
	} 
	
	# if outside j-n values
	pm.mx1.ROS1o <- if(jn.mx1.inside == FALSE & !is.na(jn.mx1.1)) {
		pm[pm.rn <= jn.mx1.1,] 
	}
	
	pm.mx1.ROS2o <- if(jn.mx1.inside == FALSE & !is.na(jn.mx1.2)) {
		pm[pm.rn >= jn.mx1.2,] 
	}
	if(!is.na(jn.mx1.1) & !is.na(jn.mx1.2)) {
		mx1.pad <- matrix(nrow = pm.nrow - nrow(pm.mx1.ROS1o) - nrow(pm.mx1.ROS2o), ncol = pm.ncol)
		pm.mx1.wROS <- rbind(pm.mx1.ROS1o, mx1.pad, pm.mx1.ROS2o)
	} else {
		if (!is.na(jn.mx1.1) & is.na(jn.mx1.2)) {
			mx1.pad <- matrix(nrow = pm.nrow - nrow(pm.mx1.ROS1o), ncol = pm.ncol)
			pm.mx1.wROS <- rbind(pm.mx1.ROS1o, mx1.pad)
		} else {
			if (is.na(jn.mx1.1) & !is.na(jn.mx1.2)) {
				mx1.pad <- matrix(nrow = pm.nrow - nrow(pm.mx1.ROS2o), ncol = pm.ncol)
				pm.mx1.wROS <- rbind(mx1.pad, pm.mx1.ROS2o)
			}
		}
	}
	
	rownames(pm.mx1.wROS) <- vecx1
	colnames(pm.mx1.wROS) <- vecx2
	
	# for x2 as mod
	pm.cn <- colnames(pm)
	pm.cn <- as.numeric(pm.cn)
	jn.mx2.1cn <- match(jn.mx2.1, colnames(pm))
	jn.mx2.2cn <- match(jn.mx2.2, colnames(pm))
	
	# if inside j-n values
	pm.mx2.ROS1i <- if(jn.mx2.inside == TRUE) {
		pm[ ,jn.mx2.1cn:jn.mx2.2cn]
	}
	
	# if outside j-n values
	pm.mx2.ROS1o <- if(jn.mx2.inside == FALSE) {
		pm[,pm.cn <= jn.mx2.1] 
	}
	
	pm.mx2.ROS2o <- if(jn.mx2.inside == FALSE) {
		pm[,pm.cn >= jn.mx2.2] 
	}
	if(!is.na(jn.mx2.1) & !is.na(jn.mx2.2)) {
		mx2.pad <- matrix(nrow = pm.nrow, ncol = pm.ncol - ncol(pm.mx2.ROS1o) - ncol(pm.mx2.ROS2o))
		pm.mx2.wROS <- cbind(pm.mx2.ROS1o, mx2.pad, pm.mx2.ROS2o)
	} else {
		if (!is.na(jn.mx2.1) & is.na(jn.mx2.2)) {
			mx2.pad <- matrix(nrow = pm.nrow, ncol = pm.ncol - ncol(pm.mx2.ROS1o))
			pm.mx2.wROS <- cbind(pm.mx2.ROS1o, mx2.pad)
		} else {
			if (is.na(jn.mx2.1) & !is.na(jn.mx2.2)) {
				mx2.pad <- matrix(nrow = pm.nrow, ncol = pm.ncol - ncol(pm.mx2.ROS2o))
				pm.mx2.wROS <- cbind(mx2.pad, pm.mx2.ROS2o)
			}
		}
	}
	
	rownames(pm.mx2.wROS) <- vecx1
	colnames(pm.mx2.wROS) <- vecx2
	
	# plot formatting ---------------------------------------------------------
	
	# gradient matrix reflecting j-n slope confidence bands ------------
	# x1 as mod
	if(jn.gradient.mx1 == "on") {
		jn.mx1cb <- johnson_neyman_veccb2(mod, pred = x2, modx = x1, 
			mod.range = c(minx1, maxx1), alpha = 0.05, plot = FALSE, vec = vecx1)$cbands
		jn.mx1cb <- jn.mx1cb[pm.rn <= jn.mx1.1 | pm.rn >= jn.mx1.2,]
		jn.mx1cb <- jn.mx1cb[!is.na(jn.mx1cb$x1),]
		jn.mx1cbd <- (jn.mx1cb[,4] - jn.mx1cb[,2])*2 # upper CI minus slope times 2
		pm.mx1.wROSg <- pm.mx1.wROS
		pm.mx1.wROSg[!is.na(pm.mx1.wROSg)] <- jn.mx1cbd
		shscmx1 <- TRUE
	} else {
		pm.mx1.wROSg <- pm.mx1.wROS
		pm.mx1.wROSg[!is.na(pm.mx1.wROSg)] <- 1	
		shscmx1 <- FALSE
	}
	
	# x2 as mod
	if(jn.gradient.mx2 == "on") {
		jn.mx2cb <- johnson_neyman_veccb2(mod, pred = x1, modx = x2, 
			mod.range = c(minx2, maxx2), alpha = 0.05, plot = FALSE, vec = vecx2)$cbands
		jn.mx2cb <- jn.mx2cb[pm.cn <= jn.mx2.1 | pm.cn >= jn.mx2.2,]
		jn.mx2cb <- jn.mx2cb[!is.na(jn.mx2cb$x2),]
		jn.mx2cbd <- (jn.mx2cb[,4] - jn.mx2cb[,2])*2 # upper CI minus slope times 2
		pm.mx2.wROSg <- t(pm.mx2.wROS)
		pm.mx2.wROSg[!is.na(pm.mx2.wROSg)] <- jn.mx2cbd
		pm.mx2.wROSg <- t(pm.mx2.wROSg)
		shscmx2 <- TRUE
	} else {
		pm.mx2.wROSg <- pm.mx2.wROS
		pm.mx2.wROSg[!is.na(pm.mx2.wROSg)] <- 1
		shscmx2 <- FALSE
	}
	
	# identify crossover points -----------------------------------------------
	ncoef <- length(coef(mod)) # number of coefficients
	
	# x1 as mod: -Bx1/Bx1*x3, plot on x2
	coonx2 <- -coef(mod)[2]/coef(mod)[ncoef]
	
	if(coonx2 >= minx2 & coonx2 <= maxx2) {
		coonx2 <- as.numeric(coonx2)
	} else {
		coonx2 <- NA
	}
	
	# predicted y at crossover
	coonx2y <- if (!is.na(coonx2)) {
		if (is.null(cov)) {
			newdatacoonx2y <- data.frame(vecx1, coonx2)
			names(newdatacoonx2y) <- c("x1", "x2")
		} else {
			if (length(cov) == 1) {
				x <- mean(df[,4])
				pmcov.l <- list(x) 
			} else {
				pmcov.l <- sapply(c(df[,4:ncol(df)]), mean, simplify = FALSE)
			}
			newdatacoonx2y <- data.frame(vecx1, coonx2, pmcov.l)
			names(newdatacoonx2y) <- c("x1", "x2", covn)
		}
		predict.lm(mod, newdatacoonx2y)
	}
	
	# x2 as mod: # x2 as mod: -Bx2/Bx1*x3, plot on x1
	coonx1 <- -coef(mod)[3]/coef(mod)[ncoef]
	
	if(coonx1 >= minx1 & coonx1 <= maxx1) {
		coonx1 <- as.numeric(coonx1)
	} else {
		coonx1 <- NA
	}
	
	# predicted y at crossover
	coonx1y <- 	if (!is.na(coonx1)) {
		if (is.null(cov)) {
			newdatacoonx1y <- data.frame(coonx1, vecx2)
			names(newdatacoonx1y) <- c("x1", "x2")
		} else {
			if (length(cov) == 1) {
				x <- mean(df[,4])
				pmcov.l <- list(x) 
			} else {
				pmcov.l <- sapply(c(df[,4:ncol(df)]), mean, simplify = FALSE)
			}
			newdatacoonx1y <- data.frame(coonx1, vecx2, pmcov.l)
			names(newdatacoonx1y) <- c("x1", "x2", covn)
		}
		predict.lm(mod, newdatacoonx1y)
	}
	
	# matrix for regression plane color --------------------------------
	pmcolor <- array(1L, dim(pm))
	
	# set opacity of ci and j-n planes
	# interaction significant? j-n values within range of observed data?
	if(intsig == 1) {
		if(!is.na(jn.mx1.2) | !is.na(jn.mx1.1)) {
			pm.mx1.wROSop <- 0.7
		} else {
			pm.mx1.wROSop <- 0
		}
		if(!is.na(jn.mx2.2) | !is.na(jn.mx2.1)) {
			pm.mx2.wROSop <- 0.7
		} else {
			pm.mx2.wROSop <- 0
		}
	} else {
		pm.mx1.wROSop <- 0
		pm.mx2.wROSop <- 0
	}
	
	# ci
	if(pred.val.ci == "on") {
		pmciop <- 0.4
	} else {
		pmciop <- 0
	}
	
	# set opacity of scatter plot
	if(scatter == "on") {
		scatterop <- 1
	} else {
		scatterop <- 0
	}
	
	
	# plot --------------------------------------------------------------------
	jn.plot <- plot_ly() %>%
		# scatterplot
		add_trace(data = df, x = ~x1,
			y = ~x2,
			z = ~y,
			type = "scatter3d",
			mode = "markers",
			marker = list(size = 1.5, color = "grey"),
			opacity = scatterop) %>%
		# format
		layout(title = plot.name,
			titlefont = list(family = "Helvetica",
				color = "#484747"),
			scene = list(
				camera = list(eye = list(x = 1.6, y = -1.6, z = 1.25)),
				xaxis = list(title = x1.name, 
					titlefont = list(family = "Arial, sans-serif", color = "#6666ff", size = 13),
					tickfont = list(family = "Arial, sans-serif", size = 11)),
				yaxis = list(title = x2.name, 
					titlefont = list(family = "Arial, sans-serif", color = "#86b300", size = 13),
					tickfont = list(family = "Arial, sans-serif", size = 11)),
				zaxis = list(title = y.name, 
					titlefont = list(family = "Arial, sans-serif", color = "#484747", size = 13),
					tickfont = list(family = "Arial, sans-serif",	color = "#484747", size = 11)))) %>%
		hide_legend()
	
	# add regression plane
	if (regplane == "on") {
		jn.plot <- add_surface(jn.plot, x = ~vecx1, y = ~vecx2, z = ~t(pm),
			type = "surface",
			surfacecolor= t(pmcolor), 
			colorscale = list(c(0, 1), c('gray', 'gray')),
			opacity = 0.55,
			showscale = FALSE) 
	}
	
	# add lower confidence band regression plane
	if (pred.val.ci == "on") {
		jn.plot <- add_surface(jn.plot, x = ~vecx1, y = ~vecx2, z = ~t(pmcil),
			type = "surface",
			surfacecolor=t(pmcolor),
			colorscale = list(c(0, 1), c('gray', 'gray')),
			opacity = pmciop,
			showscale = FALSE) %>%
			# add upper confidence band regression plane
			add_surface(x = ~vecx1, y = ~vecx2, z = ~t(pmciu),
				type = "surface",
				surfacecolor= t(pmcolor),
				colorscale = list(c(0, 1), c('gray', 'gray')),
				opacity = pmciop,
				showscale = FALSE) 
	}
	
	# add j-n plane for x2 as mod (x2 first in code because of how coloring is affected)
	if (jnROSx2mod == "on") {
		jn.plot <- add_surface(jn.plot, x = ~vecx1, y = ~vecx2, z = ~t(pm.mx2.wROS),
			type = "surface",
			surfacecolor=t(pm.mx2.wROSg),
			colorscale = list(c(0, 1), c("#0000ff", 'white')), 
			opacity = pm.mx2.wROSop,
			showscale = shscmx2,
			colorbar = list(xanchor = "right", yanchor = "middle",
				title = paste0("Width of<br>95% CB Around<br>", x1.name, " Slope"),
				titlefont = list(family = "Arial, sans-serif", size = 11),
				tickfont = list(family = "Arial, sans-serif", size = 11)))
	}
	
	# add j-n plane for x1 as mod
	if (jnROSx1mod == "on") {
		jn.plot <- add_surface(jn.plot, x = ~vecx1, y = ~vecx2, z = ~t(pm.mx1.wROS),
			type = "surface",
			surfacecolor=t(pm.mx1.wROSg),
			colorscale = list(c(0, 1), c("#99cc00", 'white')),
			opacity = pm.mx1.wROSop,
			showscale = shscmx1,
			colorbar = list(xanchor = "left", yanchor = "middle",
				title = paste0("Width of<br>95% CB Around<br>", x2.name, " Slope"),
				titlefont = list(family = "Arial, sans-serif", size = 11),
				tickfont = list(family = "Arial, sans-serif", size = 11)))
	}
	
	
	# add simple slopes for figure 9
	# add simple slopes when x1 is mod (Figure 9)
	jn.plot <- add_trace(jn.plot, x = ~vecx1[101], # jn upper bound
		y = ~vecx2, z = ~pm[101,],
		type = "scatter3d",
		mode = "lines",
		line = list(width = 3, color = "#336600")) %>%
		add_trace(x = ~vecx1[54],
			y = ~vecx2, z = ~pm[54,],
			type = "scatter3d",
			mode = "lines",
			line = list(width = 3, color = "#336600", dash = "longdash")) %>%
		add_trace(x = ~vecx1[36],
			y = ~vecx2, z = ~pm[36,],
			type = "scatter3d",
			mode = "lines",
			line = list(width = 3, color = "#336600")) %>%
		add_trace(x = ~vecx1[33],
			y = ~vecx2, z = ~pm[33,],
			type = "scatter3d",
			mode = "lines",
			line = list(width = 3, color = "gray", dash = "longdash")) %>%
		add_trace(x = ~vecx1[14],
			y = ~vecx2, z = ~pm[14,],
			type = "scatter3d",
			mode = "lines",
			line = list(width = 3, color = "gray", dash = "longdash"))
	
	htmlwidgets::saveWidget(widget = jn.plot, file = paste0(plot.name, ".html"), selfcontained = TRUE)
	jn.plot
}

figure9(file.loc.name = "simulateddata.final.testing.csv", y.var = "Depression", y.name = "Depression",
	x1.var = "LifeStress", x1.name = "Life Stress",
	x2.var = "Neuroticism", x2.name = "Neuroticism", 
	cov = NULL, cov.name = NULL,
	std.y = "center", 
	std.x1 = "center", 
	std.x2 = "center", 
	jn.gradient.mx1 = "off",
	jn.gradient.mx2 = "off",
	pred.val.ci = "off",
	crossover.x1mod = "off",
	crossover.x2mod = "off",
	scatter = "off",
	regplane = "on",
	jnROSx1mod = "on",
	jnROSx2mod = "off",
	plot.name = "Life Stress x Neuroticism Predicting Depression")
