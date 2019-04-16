# 3D plot function --------------------------------------------------------
plot3dNoInt <- function(file.loc.name, y.var, y.name, x1.var, x1.name, x2.var, x2.name, 
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
	mod <- lm(y ~ . + x1 + x2, data = df)
	mod.summary <- summary(mod)
	print(mod.summary)

	# identify min and max of x1 and x2 and round for later printing
	minx1 <- min(df$x1, na.rm = TRUE)
	maxx1 <- max(df$x1, na.rm = TRUE)
	minx2 <- min(df$x2, na.rm = TRUE)
	maxx2 <- max(df$x2, na.rm = TRUE)
	
	r.minx1 <- format(round(minx1, digits = 2), nsmall = 2)
	r.maxx1 <- format(round(maxx1, digits = 2), nsmall = 2)
	r.minx2 <- format(round(minx2, digits = 2), nsmall = 2)
	r.maxx2 <- format(round(maxx2, digits = 2), nsmall = 2)

	# vectors of x1 and x2 -----------------------------------
	vecx1 <- seq(minx1, maxx1, length.out = 100)
	vecx2 <- seq(minx2, maxx2, length.out = 100)
	
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
	
	# matrix for regression plane color --------------------------------
	pmcolor <- array(1L, dim(pm))
	
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
					titlefont = list(family = "Arial, sans-serif", color = "#3366ff", size = 13),
					tickfont = list(family = "Arial, sans-serif", size = 11)),
				yaxis = list(title = x2.name, 
					titlefont = list(family = "Arial, sans-serif", color = "#567001", size = 13),
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
			opacity = 0.6,
			showscale = FALSE) 
	}
	
	htmlwidgets::saveWidget(widget = jn.plot, file = paste0(plot.name, ".html"), selfcontained = TRUE)
	jn.plot
}

plot3dNoInt(file.loc.name = "simulateddata.final.testing.csv", y.var = "Depression", y.name = "Depression",
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
	scatter = "on",
	regplane = "on",
	jnROSx1mod = "off",
	jnROSx2mod = "off",
	plot.name = "Life Stress + Neuroticism Predicting Depression <br> (No Interaction Term)")
