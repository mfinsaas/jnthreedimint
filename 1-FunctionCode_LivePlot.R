# Function Code for Live 3D Interaction Plot
# Run this code in your R console. Then run the Input code after updating the arguments to create the 3D plot.
# 3D plot function --------------------------------------------------------
plot3dInt <- function(file.loc.name, 	plot.name,
	y.var, 
	y.printed.name, 
	x1.var, 
	x1.printed.name, 
	x2.var, 
	x2.printed.name, 
	cov = NULL, 
	std.y = c("raw", "center", "standardize"), 
	std.x1 = c("raw", "center", "standardize"),
	std.x2 = c("raw", "center", "standardize"), 
	regplane = c("on", "off"),
	scatter = c("on", "off"),
	pred.val.ci = c("on", "off"), 
	jnROSx1mod = c("solid", "gradient", "off"),
	jnROSx2mod = c("solid", "gradient", "off"),
	crossover.x1mod = c("on", "off"), 
	crossover.x2mod = c("on", "off")) {
	
	# packages required
	require(plotly)
	require(htmlwidgets)
	require(jtools)
	require(plyr)
	require(foreign)
	require(xlsx)
	require(psych)
	
	# helper functions --------------------------------------------------------
	# returns whether ROS falls between J-N bounds
	johnson_neyman_inside <- source("/Users/mfinsaas/klein_lab_projects/jn_threedimint/src/johnson_neyman_inside.R")	
	
	# adapts johnson-neyman from jtools by using vec as the moderator values and 
	# testing slopes/cbs across these values
	johnson_neyman_veccb2 <- source("/Users/mfinsaas/klein_lab_projects/jn_threedimint/src/johnson_neyman_veccb2.R")	
	
	# read in file ------------------------------------------------------------
	df <- if(tools::file_ext(file.loc.name) == "sav") {
		foreign::read.spss(file.loc.name, to.data.frame = TRUE)
	} else {
		if(tools::file_ext(file.loc.name) == "csv") {
			read.csv(file.loc.name, header = TRUE)
		}
	}
	
	# subset to variables in model
	df <- df[,c(y.var, x1.var, x2.var, cov)] 
	
	# rename variables and keep only complete cases
	if (!is.null(cov)) {
		covn <- paste0("c",seq_along(cov))
		colnames(df) <- c("y", "x1", "x2", covn)
		df[complete.cases(df), ] 
	} else {
		colnames(df) <- c("y", "x1", "x2")
		df[complete.cases(df), ] 
	}
	
	# descriptives ------------------------------------------------------------
	ds <- describe(df)[c("n","mean","sd","median","min","max","skew","kurtosis")]
	print("Descriptives")
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
	
	# covariates always centered
	if (!is.null(cov)) {
		df[,4:ncol(df)] <- c(scale(df[,4:ncol(df)] , scale = FALSE)) 
	}
	
	# identify min and max of x1 and x2 and round for later printing ----------
	
	minx1 <- min(df$x1, na.rm = TRUE)
	maxx1 <- max(df$x1, na.rm = TRUE)
	minx2 <- min(df$x2, na.rm = TRUE)
	maxx2 <- max(df$x2, na.rm = TRUE)
	
	r.minx1 <- format(round(minx1, digits = 2), nsmall = 2)
	r.maxx1 <- format(round(maxx1, digits = 2), nsmall = 2)
	r.minx2 <- format(round(minx2, digits = 2), nsmall = 2)
	r.maxx2 <- format(round(maxx2, digits = 2), nsmall = 2)
	
	# scatterplot ----------------------------------------
	
	# set opacity of scatter plot
	if(scatter == "on") {
		scatterop <- 1
	} else {
		scatterop <- 0
	}
	
	# add scatter to plot (if "off" will be transluscent but still provide range for axes)
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
				xaxis = list(title = x1.printed.name, 
					titlefont = list(family = "Arial, sans-serif", color = "#567001", size = 13),
					tickfont = list(family = "Arial, sans-serif", size = 11)),
				yaxis = list(title = x2.printed.name, 
					titlefont = list(family = "Arial, sans-serif", color = "#3366ff", size = 13),
					tickfont = list(family = "Arial, sans-serif", size = 11)),
				zaxis = list(title = y.printed.name, 
					titlefont = list(family = "Arial, sans-serif", color = "#484747", size = 13),
					tickfont = list(family = "Arial, sans-serif",	color = "#484747", size = 11)))) %>%
		hide_legend()	
	
	# fit regression model ----------------
	
	mod <- lm(y ~ . + x1 * x2, data = df)
	mod.summary <- summary(mod)
	print("Regression Model Summary")
	print(mod.summary)
	
	if(mod.summary$coefficients[nrow(mod.summary$coefficients), 4] < 0.05) {
		intsig <- 1 
	} else {
		intsig <- 0
	}
	
	# johnson neyman analysis -------------------------------------------------
	
	# when x1 is mod ------------------------------------------------------
	jn.mx1 <- johnson_neyman(mod, pred = x2, modx = x1, 
		mod.range = c(minx1, maxx1), alpha = 0.05, plot = FALSE)
	print("Johnson Neyman Results When X1 is Moderator")
	print(jn.mx1)
	
	# save jn values if within mod range
	if(jn.mx1$bounds[1] >= minx1 & jn.mx1$bounds[1] <= maxx1) {
		jn.mx1.1 <- as.numeric(jn.mx1$bounds[1])
	} else {
		jn.mx1.1 <- NA
	}
	
	if(jn.mx1$bounds[2] >= minx1 & jn.mx1$bounds[2] <= maxx1) {
		jn.mx1.2 <- as.numeric(jn.mx1$bounds[2])
	} else {
		jn.mx1.2 <- NA
	}
	
	# when x2 is mod ------------------------------------------------------
	jn.mx2 <- johnson_neyman(mod, pred = x1, modx = x2, 
		mod.range = c(minx2, maxx2), alpha = 0.05, plot = FALSE)
	print("Johnson-Neyman Results When X2 is Moderator")
	print(jn.mx2)
	
	# save jn values if within mod range
	if(jn.mx2$bounds[1] >= minx2 & jn.mx2$bounds[1] <= maxx2) {
		jn.mx2.1 <- as.numeric(jn.mx2$bounds[1])
	} else {
		jn.mx2.1 <- NA
	}
	
	if(jn.mx2$bounds[2] >= minx2 & jn.mx2$bounds[2] <= maxx2) { 
		jn.mx2.2 <- as.numeric(jn.mx2$bounds[2])
	} else {
		jn.mx2.2 <- NA
	}
	
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
	
	# add regression plane to plot ----------------------------------------------------
	
	# matrix for regression plane color 
	pmcolor <- array(1L, dim(pm))
	
	if (regplane == "on") {
		jn.plot <- add_surface(jn.plot, x = ~vecx1, y = ~vecx2, z = ~t(pm),
			type = "surface",
			surfacecolor= t(pmcolor), 
			colorscale = list(c(0, 1), c('gray', 'gray')),
			opacity = 0.55,
			showscale = FALSE) 
	}
	
	# confidence band matrices  -------------------------------------
	
	if (pred.val.ci == "on") {
		# lower ci 
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
		
		# upper ci 
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
		
		# add pred val CI to plot -------------------------------------------------
		
		jn.plot <- add_surface(jn.plot, x = ~vecx1, y = ~vecx2, z = ~t(pmcil), # lower
			type = "surface",
			surfacecolor=t(pmcolor),
			colorscale = list(c(0, 1), c('gray', 'gray')),
			opacity = 0.4,
			showscale = FALSE) %>%
			add_surface(x = ~vecx1, y = ~vecx2, z = ~t(pmciu), # upper
				type = "surface",
				surfacecolor= t(pmcolor),
				colorscale = list(c(0, 1), c('gray', 'gray')),
				opacity = 0.4,
				showscale = FALSE) 
	}
	
	# ROS matrices --------------------------------------------
	
	# for x1 as mod
	# identify whether ROS is inside j-n values
	jn.mx1.inside <- johnson_neyman_inside(mod, pred = x2, modx = x1, 
		mod.range = c(minx1, maxx1), alpha = 0.05)
	
	# identify whether ROS covers whole range of x1
	if((jn.mx1.inside == FALSE & (jn.mx1$bounds[2] < minx1 | jn.mx1$bounds[1] > maxx1)) |
			(jn.mx1.inside == TRUE & jn.mx1$bounds[1] < minx1 & jn.mx1$bounds[2] > maxx1)) {
		pm.mx1.wROS <- pm 
		mx1.ROS.all <- TRUE
	} else {
		mx1.ROS.all <- FALSE
	}
	
	# identify whether ROS within range of x1
	if((is.na(jn.mx1.1) & is.na(jn.mx1.2)) & mx1.ROS.all == FALSE) {
		mx1.ROS.none <- TRUE
	} else {
		mx1.ROS.none <- FALSE
	}
	
	# make x1 ROS planes
	if((jnROSx1mod == "solid" | jnROSx1mod == "gradient") & intsig == 1 & mx1.ROS.none == FALSE) {
		# ROS inside JN bounds
		if(jn.mx1.inside == TRUE) {
			# inside both j-n values (center band) (not run, program does not accept squared terms)
			if(!is.na(jn.mx1.1) & !is.na(jn.mx1.2)) {
				pm.mx1.ROS1i <- pm[vecx1 >= jn.mx1.1 & vecx1 <= jn.mx1.2,]
				jn.mx1.1cn <- match(jn.mx1.1, vecx1)
				jn.mx1.2cn <- match(jn.mx1.2, vecx1)
				mx1.pad1 <- matrix(nrow = jn.mx1.1cn, ncol = ncol(pm)) # lower pad made of NAs
				mx1.pad2 <- matrix(nrow = nrow(pm) - jn.mx1.2cn, ncol = ncol(pm)) # upper pad  made of NAs
				pm.mx1.wROS <- rbind(mx1.pad1, pm.mx1.ROS1i, mx1.pad2)
			}
			# lower JN bound within range (high band)
			if(!is.na(jn.mx1.1) & is.na(jn.mx1.2)) {
				pm.mx1.ROS1i <- pm[vecx1 >= jn.mx1.1,]  # ROS covers first JN bound through maximum 
				mx1.pad <- matrix(nrow = nrow(pm) - nrow(pm.mx1.ROS1i), ncol = ncol(pm))
				pm.mx1.wROS <- rbind(mx1.pad, pm.mx1.ROS1i)
			}
			# upper JN bound within range (low band)
			if(!is.na(jn.mx1.2) & is.na(jn.mx1.1)) {
				pm.mx1.ROS2i <- pm[vecx1 <= jn.mx1.2,] # ROS covers minimum through second JN bound
				mx1.pad <- matrix(nrow = nrow(pm) - nrow(pm.mx1.ROS2i), ncol = ncol(pm))
				pm.mx1.wROS <- rbind(pm.mx1.ROS2i, mx1.pad)
			}
		}
		
		# ROS outside JN bounds
		if(jn.mx1.inside == FALSE) {
			# first JN bound within range (low band)
			if(!is.na(jn.mx1.1)) { 
				pm.mx1.ROS1o <- pm[vecx1 <= jn.mx1.1,] # ROS covers minimum through first JN bound
			}
			# second JN bound within range (high band)
			if(!is.na(jn.mx1.2)) { 
				pm.mx1.ROS2o <- pm[vecx1 >= jn.mx1.2,] # ROS covers JN bound through maximum
			}
			
			# combining with pads
			# both within range (low and high)
			if(!is.na(jn.mx1.1) & !is.na(jn.mx1.2)) { 
				mx1.pad <- matrix(nrow = nrow(pm) - nrow(pm.mx1.ROS1o) - nrow(pm.mx1.ROS2o), ncol = ncol(pm)) # pad made of NAs
				pm.mx1.wROS <- rbind(pm.mx1.ROS1o, mx1.pad, pm.mx1.ROS2o) # pad between ROS with extra rows
			} 
			# first within range only (low band only)
			if (!is.na(jn.mx1.1) & is.na(jn.mx1.2)) { 
				mx1.pad <- matrix(nrow = nrow(pm) - nrow(pm.mx1.ROS1o), ncol = ncol(pm))
				pm.mx1.wROS <- rbind(pm.mx1.ROS1o, mx1.pad)
			}
			# second within range only (high band only)
			if (is.na(jn.mx1.1) & !is.na(jn.mx1.2)) { 
				mx1.pad <- matrix(nrow = nrow(pm) - nrow(pm.mx1.ROS2o), ncol = ncol(pm))
				pm.mx1.wROS <- rbind(mx1.pad, pm.mx1.ROS2o)
			}
		}
		
		# add vector values to row and column names of matrix
		rownames(pm.mx1.wROS) <- vecx1
		colnames(pm.mx1.wROS) <- vecx2
		
		# shade ROS plane
		if(jnROSx1mod == "solid") {
			pm.mx1.wROSg <- pm.mx1.wROS
			pm.mx1.wROSg[!is.na(pm.mx1.wROSg)] <- 1	
			shscmx1 <- FALSE
		}
		
		if(jnROSx1mod == "gradient") {
			jn.mx1cb <- johnson_neyman_veccb2(mod, pred = x2, modx = x1, 
				mod.range = c(minx1, maxx1), alpha = 0.05, plot = FALSE, vec = vecx1)$cbands
			if(mx1.ROS.all == TRUE) {
				jn.mx1cb <- jn.mx1cb[!is.na(jn.mx1cb$x1),]
			} else {
				if(jn.mx1.inside == FALSE) {
					jn.mx1cb <- jn.mx1cb[vecx1 <= jn.mx1.1 | vecx1 >= jn.mx1.2,]
				}
				if(jn.mx1.inside == TRUE) {
					jn.mx1cb <- jn.mx1cb[vecx1 >= jn.mx1.1 | vecx1 <= jn.mx1.2,]
				}
				jn.mx1cb <- jn.mx1cb[!is.na(jn.mx1cb$x1),]
			}
			jn.mx1cbd <- (jn.mx1cb[,4] - jn.mx1cb[,2])*2 # upper CI minus slope times 2
			pm.mx1.wROSg <- pm.mx1.wROS
			pm.mx1.wROSg[!is.na(pm.mx1.wROSg)] <- jn.mx1cbd
			shscmx1 <- TRUE
		}
		
		# add j-n plane for x1 as mod to plot
		jn.plot <- add_surface(jn.plot, x = ~vecx1, y = ~vecx2, z = ~t(pm.mx1.wROS),
			type = "surface",
			surfacecolor=t(pm.mx1.wROSg),
			colorscale = list(c(0, 1), c("#3366ff", 'white')),
			opacity = 0.7,
			showscale = shscmx1,
			colorbar = list(xanchor = "left", yanchor = "middle",
				title = paste0("Width of<br>95% CB Around<br>", x2.printed.name, " Slope"),
				titlefont = list(family = "Arial, sans-serif", size = 11),
				tickfont = list(family = "Arial, sans-serif", size = 11)))
	}
	
		# for x2 as mod 
	# identify whether ROS is inside j-n values
	jn.mx2.inside <- johnson_neyman_inside(mod, pred = x1, modx = x2, 
		mod.range = c(minx2, maxx2), alpha = 0.05)
	
	# identify whether ROS covers whole range of x2
	if((jn.mx2.inside == FALSE & (jn.mx2$bounds[2] < minx2 | jn.mx2$bounds[1] > maxx2)) |
			(jn.mx2.inside == TRUE & jn.mx2$bounds[1] < minx2 & jn.mx2$bounds[2] > maxx2)) {
		pm.mx2.wROS <- pm 
		mx2.ROS.all <- TRUE
	} else {
		mx2.ROS.all <- FALSE
	}
	
	# identify whether ROS within range of x2
	if((is.na(jn.mx2.1) & is.na(jn.mx2.2)) & mx2.ROS.all == FALSE) {
		mx2.ROS.none <- TRUE
	} else {
		mx2.ROS.none <- FALSE
	}
	
	# make x2 ROS planes 
	if((jnROSx2mod == "solid" | jnROSx2mod == "gradient") & intsig == 1 & mx2.ROS.none == FALSE) {
		
		# ROS inside JN bounds
		if(jn.mx2.inside == TRUE) {
			# inside both j-n values (center band) (not run, program does not accept squared terms)
			if(!is.na(jn.mx2.1) & !is.na(jn.mx2.2)) {
				pm.mx2.ROS1i <- pm[,vecx2 >= jn.mx2.1 & vecx2 <= jn.mx2.2]
				jn.mx2.1cn <- match(jn.mx2.1, vecx2)
				jn.mx2.2cn <- match(jn.mx2.2, vecx2)
				mx2.pad1 <- matrix(nrow = nrow(pm), ncol = jn.mx2.1cn) # lower pad made of NAs
				mx2.pad2 <- matrix(nrow = nrow(pm), ncol = ncol(pm)- jn.mx2.2cn) # upper pad  made of NAs
				pm.mx2.wROS <- cbind(mx2.pad1, pm.mx2.ROS1i, mx2.pad2)
			}
			# lower JN bound within range (high band)
			if(!is.na(jn.mx2.1) & is.na(jn.mx2.2)) {
				pm.mx2.ROS1i <- pm[,vecx2 >= jn.mx2.1]  # ROS covers first JN bound through maximum 
				mx2.pad <- matrix(nrow = nrow(pm), ncol = ncol(pm) - ncol(pm.mx2.ROS1i))
				pm.mx2.wROS <- cbind(mx2.pad, pm.mx2.ROS1i)
			}
			# upper JN bound within range (low band)
			if(!is.na(jn.mx2.2) & is.na(jn.mx2.1)) {
				pm.mx2.ROS2i <-	pm[,vecx2 <= jn.mx2.2] # ROS covers minimum through second JN bound
				mx2.pad <- matrix(nrow = nrow(pm), ncol = ncol(pm) - ncol(pm.mx2.ROS2i))
				pm.mx2.wROS <- cbind(pm.mx2.ROS2i, mx2.pad)
			}
		}
		
		# ROS outside JN bounds
		if(jn.mx2.inside == FALSE) {
			
			# first JN bound within range (low band)
			if(!is.na(jn.mx2.1)) { 
				pm.mx2.ROS1o <- pm[,vecx2 <= jn.mx2.1] # ROS covers minimum through first JN bound
			}
			# second JN bound within range (high band)
			if(!is.na(jn.mx2.2)) { 
				pm.mx2.ROS2o <- pm[,vecx2 >= jn.mx2.2] # ROS covers JN bound through maximum
			}
			
			# combining with pads
			# both within range (low and high)
			if(!is.na(jn.mx2.1) & !is.na(jn.mx2.2)) { 
				mx2.pad <- matrix(nrow = nrow(pm), ncol = ncol(pm) - ncol(pm.mx2.ROS1o) - ncol(pm.mx2.ROS2o))
				pm.mx2.wROS <- cbind(pm.mx2.ROS1o, mx2.pad, pm.mx2.ROS2o)
			} 
			# first within range only (low band only)
			if (!is.na(jn.mx2.1) & is.na(jn.mx2.2)) { 
				mx2.pad <- matrix(nrow = nrow(pm), ncol = ncol(pm) - ncol(pm.mx2.ROS1o))
				pm.mx2.wROS <- cbind(pm.mx2.ROS1o, mx2.pad)
			}
			# second within range only (high band only)
			if (is.na(jn.mx2.1) & !is.na(jn.mx2.2)) { 
				mx2.pad <- matrix(nrow = nrow(pm), ncol = ncol(pm) - ncol(pm.mx2.ROS2o))
				pm.mx2.wROS <- cbind(mx2.pad, pm.mx2.ROS2o)
			}
		}
		
		# add vector values to row and column names of matrix
		rownames(pm.mx2.wROS) <- vecx1
		colnames(pm.mx2.wROS) <- vecx2
		
		# shade ROS plane
		if(jnROSx2mod == "solid") {
			pm.mx2.wROSg <- pm.mx2.wROS
			pm.mx2.wROSg[!is.na(pm.mx2.wROSg)] <- 1	
			shscmx2 <- FALSE
		}
		
		if(jnROSx2mod == "gradient") {
			jn.mx2cb <- johnson_neyman_veccb2(mod, pred = x1, modx = x2, 
				mod.range = c(minx2, maxx2), alpha = 0.05, plot = FALSE, vec = vecx2)$cbands
			if(mx2.ROS.all == TRUE) {
				jn.mx2cb <- jn.mx2cb[!is.na(jn.mx2cb$x2),]
			} else {
				if(jn.mx2.inside == FALSE) {
					jn.mx2cb <- jn.mx2cb[vecx2 <= jn.mx2.1 | vecx2 >= jn.mx2.2,]
				} 
				if (jn.mx2.inside == TRUE) {
					jn.mx2cb <- jn.mx2cb[vecx2 >= jn.mx2.1 | vecx2 <= jn.mx2.2,]
				}
				jn.mx2cb <- jn.mx2cb[!is.na(jn.mx2cb$x2),]
			}
			jn.mx2cbd <- (jn.mx2cb[,4] - jn.mx2cb[,2])*2 # upper CI minus slope times 2
			pm.mx2.wROSg <- t(pm.mx2.wROS)
			pm.mx2.wROSg[!is.na(pm.mx2.wROSg)] <- jn.mx2cbd
			pm.mx2.wROSg <- t(pm.mx2.wROSg)
			shscmx2 <- TRUE
		}
		
		jn.plot <- add_surface(jn.plot, x = ~vecx1, y = ~vecx2, z = ~t(pm.mx2.wROS),
			type = "surface",
			surfacecolor=t(pm.mx2.wROSg),
			colorscale = list(c(0, 1), c("#567001", 'white')), 
			opacity = 0.7,
			showscale = shscmx2,
			colorbar = list(xanchor = "right", yanchor = "middle",
				title = paste0("Width of<br>95% CB Around<br>", x1.printed.name, " Slope"),
				titlefont = list(family = "Arial, sans-serif", size = 11),
				tickfont = list(family = "Arial, sans-serif", size = 11)))
	}
	
	# identify crossover points -----------------------------------------------
	
	ncoef <- length(coef(mod)) # number of coefficients
	
	# for x1 as mod
	if(crossover.x1mod == "on") {
		coonx2 <- -coef(mod)[2]/coef(mod)[ncoef] # -Bx1/Bx1*x2, plot on x2
		
		# identify whether it falls within range
		if(coonx2 >= minx2 & coonx2 <= maxx2) {
			coonx2 <- as.numeric(coonx2)
			print(c("Crossover when X1 is moderator", round(coonx2, digits = 2)))
		} else {
			coonx2 <- NA
			print(c("Crossover when X1 is moderator", "Not within range of data"))
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
		
		# add crossover line to plot for x1 as mod
		if(!is.na(coonx2)) {
			jn.plot <- add_trace(jn.plot, x = vecx1, y = ~coonx2, z = ~coonx2y,
				type = "scatter3d",
				mode = "lines",
				line = list(width = 3, dash = "solid", color = "#567001"),
				opacity = .8)
		}
	}
	
	# for x2 as mod
	if(crossover.x2mod == "on") {
		coonx1 <- -coef(mod)[3]/coef(mod)[ncoef] #-Bx2/Bx1*x2, plot on x1
		
		# identify whether it falls within range
		if(coonx1 >= minx1 & coonx1 <= maxx1) {
			coonx1 <- as.numeric(coonx1)
			print(c("Crossover when X2 is moderator", round(coonx1, digits = 2)))
		} else {
			coonx1 <- NA
			print(c("Crossover when X2 is moderator", "Not within range of data"))
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
		# add crossover line to plot for x2 as mod
		if(!is.na(coonx1)) {
			jn.plot <- add_trace(jn.plot, x = coonx1, y = ~vecx2, z = ~coonx1y,
				type = "scatter3d",
				mode = "lines",
				line = list(width = 3, dash = "solid", color = "#3366ff"),
				opacity = .8)
		}
	}
	
	# save plot and send to viewer
	htmlwidgets::saveWidget(widget = jn.plot, file = paste0(plot.name, ".html"), selfcontained = TRUE)
	jn.plot
}

