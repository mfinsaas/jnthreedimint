require(shiny)
require(plotly)
require(foreign)
require(plyr)
require(dplyr)
require(htmlwidgets)
require(jtools)
require(psych)
require(stats)
require(DT)
library(shinyjs)
library(htmlTable)

# helper functions --------------------------------------------------------
# modified JN, also returns whether ROS falues between JN values
johnson_neyman_new <- function (pred, modx, vmat = NULL, alpha = 0.05, plot = TRUE, 
	control.fdr = FALSE, line.thickness = 0.5, df = "residual", 
	digits = getOption("jtools-digits", 2), critical.t = NULL, 
	sig.color = "#00BFC4", insig.color = "#F8766D", mod.range = NULL, 
	title = "Johnson-Neyman plot",
	df.residual.mod,
	coef.model,
	names.coef.model) # added arguments 
{
	if (df == "residual") {
		# df <- df.residual(model)
		df <- df.residual.mod
	}
	else if (df == "normal") {
		df <- Inf
	}
	else if (!is.numeric(df)) {
		stop("df argument must be 'residual', 'normal', or a number.")
	}
	out <- list()
	out <- structure(out, pred = pred, modx = modx, 
		alpha = alpha, plot = plot, digits = digits, control.fdr = control.fdr)
	intterm1 <- paste(pred, ":", modx, sep = "")
	intterm1tf <- any(intterm1 %in% names.coef.model)
	intterm2 <- paste(modx, ":", pred, sep = "")
	intterm2tf <- any(intterm2 %in% names.coef.model)
	coefs <- coef.model
	inttermstf <- c(intterm1tf, intterm2tf)
	intterms <- c(intterm1, intterm2)
	intterm <- intterms[which(inttermstf)]
	modrange <- mod.range
	modrangeo <- mod.range
	no_range_line <- TRUE
	alpha <- alpha/2
	if (control.fdr == FALSE) {
		tcrit <- qt(alpha, df = df)
		tcrit <- abs(tcrit)
	}
	else if (control.fdr == TRUE) {
		predb <- coefs[pred]
		intb <- coefs[intterm]
		vcovs <- vmat
		vcov_pred <- vcovs[pred, pred]
		vcov_int <- vcovs[intterm, intterm]
		vcov_pred_int <- vcovs[pred, intterm]
		range_sequence <- seq(from = modrangeo[1], to = modrangeo[2], 
			by = (modrangeo[2] - modrangeo[1])/1000)
		marginal_effects <- predb + intb * range_sequence
		me_ses <- sqrt(vcov_pred + (range_sequence^2) * vcov_int + 
				2 * range_sequence * vcov_pred_int)
		ts <- marginal_effects/me_ses
		df <- df.residual.mod
		ps <- 2 * pmin(pt(ts, df = df), (1 - pt(ts, df = df))) # not this
		multipliers <- seq_along(marginal_effects)/length(marginal_effects) # THIS IS GIVING THE PROBLEMS. 
		ps_o <- order(ps)
		test <- 0
		i <- 1 + length(marginal_effects)
		while (test == 0 & i > 1) {
			i <- i - 1
			test <- min(ps[ps_o][1:i] <= multipliers[i] * (alpha * 
					2))
		}
		tcrit <- abs(qt(multipliers[i] * alpha, df))
	}
	covy3 <- vmat[intterm, intterm]
	covy1 <- vmat[pred, pred]
	covy1y3 <- vmat[intterm, pred]
	y3 <- coef.model[intterm]
	y1 <- coef.model[pred]
	if (!is.null(critical.t)) {
		tcrit <- critical.t
	}
	a <- tcrit^2 * covy3 - y3^2
	b <- 2 * (tcrit^2 * covy1y3 - y1 * y3)
	c <- tcrit^2 * covy1 - y1^2
	discriminant <- function(a, b, c) {
		disc <- b^2 - 4 * a * c
		if (disc > 0) {
			out <- disc
		}
		else if (disc == 0) {
			return(NULL)
		}
		else {
			return(NULL)
		}
		return(out)
	}
	disc <- discriminant(a, b, c)
	if (is.null(disc)) {
		failed <- TRUE
	}
	else {
		failed <- FALSE
	}
	quadsolve <- function(a, b, c, disc) {
		x1 <- (-b + sqrt(disc))/(2 * a)
		x2 <- (-b - sqrt(disc))/(2 * a)
		result <- c(x1, x2)
		result <- sort(result, decreasing = FALSE)
		return(result)
	}
	if (!is.null(disc)) {
		bounds <- quadsolve(a, b, c, disc)
	}
	else {
		bounds <- c(-Inf, Inf)
	}
	names(bounds) <- c("Lower", "Higher")
	cbands <- function(x2, y1, y3, covy1, covy3, covy1y3, tcrit, 
		predl, modx) {
		upper <- c()
		slopes <- c()
		lower <- c()
		slopesf <- function(i) {
			s <- y1 + y3 * i
			return(s)
		}
		upperf <- function(i, s) {
			u <- s + tcrit * sqrt((covy1 + 2 * i * covy1y3 + 
					i^2 * covy3))
			return(u)
		}
		lowerf <- function(i, s) {
			l <- s - tcrit * sqrt((covy1 + 2 * i * covy1y3 + 
					i^2 * covy3))
			return(l)
		}
		slopes <- sapply(x2, slopesf, simplify = "vector", USE.NAMES = FALSE)
		upper <- mapply(upperf, x2, slopes)
		lower <- mapply(lowerf, x2, slopes)
		out <- matrix(c(x2, slopes, lower, upper), ncol = 4)
		colnames(out) <- c(modx, predl, "Lower", "Upper")
		out <- as.data.frame(out)
		return(out)
	}
	x2 <- seq(from = modrange[1], to = modrange[2], length.out = 1000)
	predl <- paste("Slope of", pred)
	cbs <- cbands(x2, y1, y3, covy1, covy3, covy1y3, tcrit, predl, 
		modx)
	out$bounds <- bounds
	out <- structure(out, modrange = modrangeo)
	sigs <- which((cbs$Lower < 0 & cbs$Upper < 0) | (cbs$Lower > 
			0 & cbs$Upper > 0))
	insigs <- setdiff(1:1000, sigs)
	cbs$Significance <- rep(NA, nrow(cbs))
	cbs$Significance <- factor(cbs$Significance, levels = c("Insignificant", 
		"Significant"))
	index <- 1:1000 %in% insigs
	cbs$Significance[index] <- "Insignificant"
	index <- 1:1000 %in% sigs
	cbs$Significance[index] <- "Significant"
	out$cbands <- cbs
	index <- which(cbs$Significance == "Significant")[1]
	if (!is.na(index) & index != 0) {
		inside <- (cbs[index, modx] > bounds[1] && cbs[index, 
			modx] < bounds[2])
		all_sig <- NULL
		if (is.na(which(cbs$Significance == "Insignificant")[1])) {
			all_sig <- TRUE
		}
		else {
			all_sig <- FALSE
		}
	}
	else {
		inside <- FALSE
		all_sig <- TRUE
	}
	out$inside <- inside # added this to return JN region falls inside or outside JN mod values
	out <- structure(out, inside = inside, failed = failed, all_sig = all_sig)
	cbso1 <- cbs[cbs[, modx] < bounds[1], ]
	cbso2 <- cbs[cbs[, modx] > bounds[2], ]
	cbsi <- cbs[(cbs[, modx] > bounds[1] & cbs[, modx] < bounds[2]), 
		]
	alpha <- alpha * 2
	alpha <- gsub("0\\.", "\\.", as.character(alpha))
	pmsg <- paste("p <", alpha)
	plot <- ggplot2::ggplot() + ggplot2::geom_path(data = cbso1, 
		ggplot2::aes(x = cbso1[, modx], y = cbso1[, predl], color = cbso1[, 
			"Significance"]), size = line.thickness) + ggplot2::geom_path(data = cbsi, 
				ggplot2::aes(x = cbsi[, modx], y = cbsi[, predl], color = cbsi[, 
					"Significance"]), size = line.thickness) + ggplot2::geom_path(data = cbso2, 
						ggplot2::aes(x = cbso2[, modx], y = cbso2[, predl], color = cbso2[, 
							"Significance"]), size = line.thickness) + ggplot2::geom_ribbon(data = cbso1, 
								ggplot2::aes(x = cbso1[, modx], ymin = cbso1[, "Lower"], 
									ymax = cbso1[, "Upper"], fill = cbso1[, "Significance"]), 
								alpha = 0.2) + ggplot2::geom_ribbon(data = cbsi, ggplot2::aes(x = cbsi[, 
									modx], ymin = cbsi[, "Lower"], ymax = cbsi[, "Upper"], 
									fill = cbsi[, "Significance"]), alpha = 0.2) + ggplot2::geom_ribbon(data = cbso2, 
										ggplot2::aes(x = cbso2[, modx], ymin = cbso2[, "Lower"], 
											ymax = cbso2[, "Upper"], fill = cbso2[, "Significance"]), 
										alpha = 0.2) + ggplot2::scale_fill_manual(values = c(Significant = sig.color, 
											Insignificant = insig.color), labels = c("n.s.", pmsg), 
											breaks = c("Insignificant", "Significant"), drop = FALSE, 
											guide = ggplot2::guide_legend(title = NULL, order = 2)) + ggplot2::geom_hline(ggplot2::aes(yintercept = 0))
	if (no_range_line == FALSE) {
		plot <- plot + ggplot2::geom_segment(ggplot2::aes(x = modrangeo[1], 
			xend = modrangeo[2], y = 0, yend = 0, linetype = "Range of\nobserved\ndata"), 
			lineend = "square", size = 1.25)
	}
	plot <- plot + ggplot2::scale_linetype_discrete(name = " ", 
		guide = ggplot2::guide_legend(order = 1))
	if (out$bounds[1] < modrange[1]) {
	}
	else if (all_sig == FALSE) {
		plot <- plot + ggplot2::geom_vline(ggplot2::aes(xintercept = out$bounds[1]), 
			linetype = 2, color = sig.color)
	}
	if (out$bounds[2] > modrange[2]) {
	}
	else if (all_sig == FALSE) {
		plot <- plot + ggplot2::geom_vline(ggplot2::aes(xintercept = out$bounds[2]), 
			linetype = 2, color = sig.color)
	}
	plot <- plot + ggplot2::xlim(range(cbs[, modx])) + ggplot2::labs(title = title, 
		x = modx, y = predl) + ggplot2::scale_color_manual(name = "", 
			values = c(Significant = sig.color, Insignificant = insig.color), 
			# guide = "none") + theme_apa(legend.pos = "right", legend.font.size = 10) + 
			guide = "none") + 
		ggplot2::theme(legend.key.size = ggplot2::unit(1, "lines"))
	out$plot <- plot
	if (control.fdr == TRUE) {
		out$t_value <- tcrit
	}
	class(out) <- "johnson_neyman"
	return(out)
}


# adapts johnson-neyman from jtools by using vec as the moderator values and 
# testing slopes/cbs across these values
johnson_neyman_veccb2 <- function (pred, modx, vmat = NULL, alpha = 0.05, plot = TRUE, 
	control.fdr = FALSE, line.thickness = 0.5, df = "residual", 
	digits = getOption("jtools-digits", 2), critical.t = NULL, 
	sig.color = "#00BFC4", insig.color = "#F8766D", mod.range = NULL, 
	title = "Johnson-Neyman plot",
	df.residual.mod,
	coef.model,
	names.coef.model,
	vec) # added arguments 
{
	if (df == "residual") {
		# df <- df.residual(model)
		df <- df.residual.mod
	}
	else if (df == "normal") {
		df <- Inf
	}
	else if (!is.numeric(df)) {
		stop("df argument must be 'residual', 'normal', or a number.")
	}
	out <- list()
	out <- structure(out, pred = pred, modx = modx, 
		alpha = alpha, plot = plot, digits = digits, control.fdr = control.fdr)
	intterm1 <- paste(pred, ":", modx, sep = "")
	intterm1tf <- any(intterm1 %in% names.coef.model)
	intterm2 <- paste(modx, ":", pred, sep = "")
	intterm2tf <- any(intterm2 %in% names.coef.model)
	coefs <- coef.model
	inttermstf <- c(intterm1tf, intterm2tf)
	intterms <- c(intterm1, intterm2)
	intterm <- intterms[which(inttermstf)]
	modrange <- mod.range
	modrangeo <- mod.range
	no_range_line <- TRUE
	alpha <- alpha/2
	if (control.fdr == FALSE) {
		tcrit <- qt(alpha, df = df)
		tcrit <- abs(tcrit)
	}
	else if (control.fdr == TRUE) {
		predb <- coefs[pred]
		intb <- coefs[intterm]
		vcovs <- vmat
		vcov_pred <- vcovs[pred, pred]
		vcov_int <- vcovs[intterm, intterm]
		vcov_pred_int <- vcovs[pred, intterm]
		# range_sequence <- seq(from = modrangeo[1], to = modrangeo[2], 
		# 	by = (modrangeo[2] - modrangeo[1])/1000)
		range_sequence <- vec
		marginal_effects <- predb + intb * range_sequence
		me_ses <- sqrt(vcov_pred + (range_sequence^2) * vcov_int + 
				2 * range_sequence * vcov_pred_int)
		ts <- marginal_effects/me_ses
		df <- df.residual.mod
		ps <- 2 * pmin(pt(ts, df = df), (1 - pt(ts, df = df)))
		multipliers <- seq_along(marginal_effects)/length(marginal_effects)
		ps_o <- order(ps)
		test <- 0
		i <- 1 + length(marginal_effects)
		while (test == 0 & i > 1) {
			i <- i - 1
			test <- min(ps[ps_o][1:i] <= multipliers[i] * (alpha * 
					2))
		}
		tcrit <- abs(qt(multipliers[i] * alpha, df))
	}
	covy3 <- vmat[intterm, intterm]
	covy1 <- vmat[pred, pred]
	covy1y3 <- vmat[intterm, pred]
	y3 <- coef.model[intterm]
	y1 <- coef.model[pred]
	if (!is.null(critical.t)) {
		tcrit <- critical.t
	}
	a <- tcrit^2 * covy3 - y3^2
	b <- 2 * (tcrit^2 * covy1y3 - y1 * y3)
	c <- tcrit^2 * covy1 - y1^2
	discriminant <- function(a, b, c) {
		disc <- b^2 - 4 * a * c
		if (disc > 0) {
			out <- disc
		}
		else if (disc == 0) {
			return(NULL)
		}
		else {
			return(NULL)
		}
		return(out)
	}
	disc <- discriminant(a, b, c)
	if (is.null(disc)) {
		failed <- TRUE
	}
	else {
		failed <- FALSE
	}
	quadsolve <- function(a, b, c, disc) {
		x1 <- (-b + sqrt(disc))/(2 * a)
		x2 <- (-b - sqrt(disc))/(2 * a)
		result <- c(x1, x2)
		result <- sort(result, decreasing = FALSE)
		return(result)
	}
	if (!is.null(disc)) {
		bounds <- quadsolve(a, b, c, disc)
	}
	else {
		bounds <- c(-Inf, Inf)
	}
	names(bounds) <- c("Lower", "Higher")
	cbands <- function(x2, y1, y3, covy1, covy3, covy1y3, tcrit, 
		predl, modx) {
		upper <- c()
		slopes <- c()
		lower <- c()
		slopesf <- function(i) {
			s <- y1 + y3 * i
			return(s)
		}
		upperf <- function(i, s) {
			u <- s + tcrit * sqrt((covy1 + 2 * i * covy1y3 + 
					i^2 * covy3))
			return(u)
		}
		lowerf <- function(i, s) {
			l <- s - tcrit * sqrt((covy1 + 2 * i * covy1y3 + 
					i^2 * covy3))
			return(l)
		}
		slopes <- sapply(x2, slopesf, simplify = "vector", USE.NAMES = FALSE)
		upper <- mapply(upperf, x2, slopes)
		lower <- mapply(lowerf, x2, slopes)
		out <- matrix(c(x2, slopes, lower, upper), ncol = 4)
		colnames(out) <- c(modx, predl, "Lower", "Upper")
		out <- as.data.frame(out)
		return(out)
	}
	# x2 <- seq(from = modrange[1], to = modrange[2], length.out = 1000)
	x2 <- vec
	predl <- paste("Slope of", pred)
	cbs <- cbands(x2, y1, y3, covy1, covy3, covy1y3, tcrit, predl, 
		modx)
	out$bounds <- bounds
	out <- structure(out, modrange = modrangeo)
	sigs <- which((cbs$Lower < 0 & cbs$Upper < 0) | (cbs$Lower > 
			0 & cbs$Upper > 0))
	# insigs <- setdiff(1:1000, sigs)
	insigs <- setdiff(1:length(vec), sigs)
	cbs$Significance <- rep(NA, nrow(cbs))
	cbs$Significance <- factor(cbs$Significance, levels = c("Insignificant", 
		"Significant"))
	# index <- 1:1000 %in% insigs
	index <- 1:length(vec) %in% insigs
	cbs$Significance[index] <- "Insignificant"
	# index <- 1:1000 %in% sigs
	index <- 1:length(vec) %in% sigs
	cbs$Significance[index] <- "Significant"
	out$cbands <- cbs
	index <- which(cbs$Significance == "Significant")[1]
	if (!is.na(index) & index != 0) {
		inside <- (cbs[index, modx] > bounds[1] && cbs[index, 
			modx] < bounds[2])
		all_sig <- NULL
		if (is.na(which(cbs$Significance == "Insignificant")[1])) {
			all_sig <- TRUE
		}
		else {
			all_sig <- FALSE
		}
	}
	else {
		inside <- FALSE
		all_sig <- TRUE
	}
	out$inside <- inside # added this to return JN region falls inside or outside JN mod values
	out <- structure(out, inside = inside, failed = failed, all_sig = all_sig)
	cbso1 <- cbs[cbs[, modx] < bounds[1], ]
	cbso2 <- cbs[cbs[, modx] > bounds[2], ]
	cbsi <- cbs[(cbs[, modx] > bounds[1] & cbs[, modx] < bounds[2]), 
		]
	alpha <- alpha * 2
	alpha <- gsub("0\\.", "\\.", as.character(alpha))
	pmsg <- paste("p <", alpha)
	plot <- ggplot2::ggplot() + ggplot2::geom_path(data = cbso1, 
		ggplot2::aes(x = cbso1[, modx], y = cbso1[, predl], color = cbso1[, 
			"Significance"]), size = line.thickness) + ggplot2::geom_path(data = cbsi, 
				ggplot2::aes(x = cbsi[, modx], y = cbsi[, predl], color = cbsi[, 
					"Significance"]), size = line.thickness) + ggplot2::geom_path(data = cbso2, 
						ggplot2::aes(x = cbso2[, modx], y = cbso2[, predl], color = cbso2[, 
							"Significance"]), size = line.thickness) + ggplot2::geom_ribbon(data = cbso1, 
								ggplot2::aes(x = cbso1[, modx], ymin = cbso1[, "Lower"], 
									ymax = cbso1[, "Upper"], fill = cbso1[, "Significance"]), 
								alpha = 0.2) + ggplot2::geom_ribbon(data = cbsi, ggplot2::aes(x = cbsi[, 
									modx], ymin = cbsi[, "Lower"], ymax = cbsi[, "Upper"], 
									fill = cbsi[, "Significance"]), alpha = 0.2) + ggplot2::geom_ribbon(data = cbso2, 
										ggplot2::aes(x = cbso2[, modx], ymin = cbso2[, "Lower"], 
											ymax = cbso2[, "Upper"], fill = cbso2[, "Significance"]), 
										alpha = 0.2) + ggplot2::scale_fill_manual(values = c(Significant = sig.color, 
											Insignificant = insig.color), labels = c("n.s.", pmsg), 
											breaks = c("Insignificant", "Significant"), drop = FALSE, 
											guide = ggplot2::guide_legend(order = 2)) + ggplot2::geom_hline(ggplot2::aes(yintercept = 0))
	if (no_range_line == FALSE) {
		plot <- plot + ggplot2::geom_segment(ggplot2::aes(x = modrangeo[1], 
			xend = modrangeo[2], y = 0, yend = 0, linetype = "Range of\nobserved\ndata"), 
			lineend = "square", size = 1.25)
	}
	plot <- plot + ggplot2::scale_linetype_discrete(name = " ", 
		guide = ggplot2::guide_legend(order = 1))
	if (out$bounds[1] < modrange[1]) {
	}
	else if (all_sig == FALSE) {
		plot <- plot + ggplot2::geom_vline(ggplot2::aes(xintercept = out$bounds[1]), 
			linetype = 2, color = sig.color)
	}
	if (out$bounds[2] > modrange[2]) {
	}
	else if (all_sig == FALSE) {
		plot <- plot + ggplot2::geom_vline(ggplot2::aes(xintercept = out$bounds[2]), 
			linetype = 2, color = sig.color)
	}
	plot <- plot + ggplot2::xlim(range(cbs[, modx])) + ggplot2::labs(title = title, 
		x = modx, y = predl) + ggplot2::scale_color_manual(name = "", 
			values = c(Significant = sig.color, Insignificant = insig.color), 
			# guide = "none") + theme_apa(legend.pos = "right", legend.font.size = 10) + 
			guide = "none") + 
		ggplot2::theme(legend.key.size = ggplot2::unit(1, "lines"))
	out$plot <- plot
	if (control.fdr == TRUE) {
		out$t_value <- tcrit
	}
	class(out) <- "johnson_neyman"
	return(out)
}

# UI ----------------------------------------------------------------------
ui <- fluidPage(
	useShinyjs(),
	titlePanel(
		h4("Visualizing Continuous x Continuous Interactions in 3D", 
			h6("M. C. Finsaas & B. L. Goldstein, Department of Psychology, Stony Brook University")),
		"3D Interactions"
	),
	sidebarLayout(
		# sidebar panel -----------------------------------------------------------
		# manual input
		sidebarPanel(width = 3,
			radioButtons(inputId = "inputType",
				label = "Select Input Method",
				choices = c("Upload raw data" = "rawData", 
					"Manual entry (beta)" = "manual"),
				inline = TRUE),
			hr(),
			fluidRow(
				column(6, uiOutput(outputId = "yManNameUI")),
				column(6, uiOutput(outputId = "x1ManNameUI")),
				column(6, uiOutput(outputId = "x2ManNameUI"))),
			htmlOutput(outputId = "line1"), 
			fluidRow(
				column(6, uiOutput(outputId = "minx1UI")),
				column(6, uiOutput(outputId = "maxx1UI")),
				column(6, uiOutput(outputId = "minx2UI")),
				column(6, uiOutput(outputId = "maxx2UI"))),
			htmlOutput(outputId = "line2"), 
			#		htmlOutput(outputId = "coefTitle"),
			fluidRow(
				column(6, uiOutput(outputId = "b0UI")),
				column(6, uiOutput(outputId = "bx1UI")),
				column(6, uiOutput(outputId = "bx2UI")),
				column(6, uiOutput(outputId = "bIntUI"))),
			htmlOutput(outputId = "line3"), 
			#		htmlOutput(outputId = "coefVarTitle"),
			fluidRow(
				column(6, uiOutput(outputId = "b0varUI")),
				column(6, uiOutput(outputId = "bx1varUI")),
				column(6, uiOutput(outputId = "bx2varUI")),
				column(6, uiOutput(outputId = "bIntvarUI"))),
			htmlOutput(outputId = "line4"), 
			#	htmlOutput(outputId = "coefCovarTitle"),
			fluidRow(
				column(6, uiOutput(outputId = "bx1IntcovUI")),
				column(6, uiOutput(outputId = "bx2IntcovUI"))),
			uiOutput(outputId = "dfsUI"),
			uiOutput(outputId = "fdrManualUI"),
			htmlOutput(outputId = "messageMinGEMaxX1"),
			htmlOutput(outputId = "messageMinGEMaxX2"),
			htmlOutput(outputId = "messageNegVar"),
			htmlOutput(outputId = "messageNotRight"),
			htmlOutput(outputId = "messageOneNA"),
			htmlOutput(outputId = "messageDf"),
			htmlOutput(outputId = "messageInt0"),
			# raw data upload
			uiOutput(outputId = "uploadFileUI"),
			htmlOutput(outputId = "line6"), 
			uiOutput(outputId = "x1UI"),
			uiOutput(outputId = "std.x1UI"),
			uiOutput(outputId = "x1PrintNameUI"),
			htmlOutput(outputId = "line7"), 
			uiOutput(outputId = "x2UI"),
			uiOutput(outputId = "std.x2UI"),
			uiOutput(outputId = "x2PrintNameUI"),
			htmlOutput(outputId = "line8"), 
			uiOutput(outputId = "yUI"),
			uiOutput(outputId = "std.yUI"),
			uiOutput(outputId = "yPrintNameUI"),
			htmlOutput(outputId = "line9"), 
			uiOutput(outputId = "covUI"),
			htmlOutput(outputId = "line10"), 
			uiOutput(outputId = "fdrRawDataUI"),
			checkboxInput(inputId = "grayscale",
				label = "Grayscale?"),
			column(12, align="center", 
				actionButton("run",label = "Run / Update")),
			br(),
			br(),
			htmlOutput(outputId = "messageNonNum"),
			htmlOutput(outputId = "messageSameVar"),
			htmlOutput(outputId = "messageInt"),
			tags$style("body{font-size: 10px}"),
			tags$style(".btn{font-size: 10px}"),
			tags$style(".form-control{font-size: 10px; height: 28px}"),
			tags$style(".well{padding: 14px}"),
			tags$style(".form-group{margin-bottom: 10px}"),
			tags$style("hr{margin-bottom: 10px; margin-top: 10px; border-top: 1px solid #cccccc}"),
			tags$style(".dropdown-menu{font-size: 10px}"),
			tags$style(".progress-bar{font-size: 10px")),
		# main panel --------------------------------------------------------------
		mainPanel(width = 9,
			navbarPage(title = "", id = "inTabset",
				tabPanel(title = "How To",
					htmlOutput(outputId = "guide")),
				tabPanel(title = "3D Plot", value = "3DPlotPanel",
					fluidRow(
						column(7,
							plotlyOutput(outputId = "plot3dPlot")),
						column(5,
							fluidRow(
								column(4, uiOutput(outputId = "regplaneUI")),
								column(4, uiOutput(outputId = "scatterUI")),
								column(4, uiOutput(outputId = "predCIUI"))),
							wellPanel(
								htmlOutput(outputId = "mx1CheckboxTitle"),
								uiOutput(outputId = "pm.mx1.wROSUI"), 
								uiOutput(outputId = "coonx1UI"),
								htmlOutput("jnmx1.sum")),
							wellPanel(
								htmlOutput(outputId = "mx2CheckboxTitle"),
								uiOutput(outputId = "pm.mx2.wROSUI"), 
								uiOutput(outputId = "coonx2UI"),
								htmlOutput("jnmx2.sum")))),
					fluidRow(
						column(4, 
							htmlOutput(outputId = "mx1cbscale.title"),
							htmlOutput(outputId = "mx1cbscale"))),
					fluidRow(
						column(1,
							htmlOutput(outputId = "mx1cbscale.text1")),
						column(3,
							htmlOutput(outputId = "mx1cbscale.text2"))),
					fluidRow(
						column(4, 
							htmlOutput(outputId = "mx2cbscale.title"),
							htmlOutput(outputId = "mx2cbscale"))),
					fluidRow(
						column(1,
							htmlOutput(outputId = "mx2cbscale.text1")),
						column(3,
							htmlOutput(outputId = "mx2cbscale.text2")))),
				navbarMenu(title = "Johnson-Neyman",
					tabPanel(title = "Marginal Effects Plots",
						htmlOutput(outputId = "jnmx1Title"), 
						br(),
						fluidRow(
							column(4, htmlOutput(outputId = "jnmx1Info")),
							column(8,	plotOutput(outputId = "jnmx1Plot"))),
						br(),
						htmlOutput(outputId = "jnmx2Title"), 
						br(),
						fluidRow(
							column(4, htmlOutput(outputId = "jnmx2Info")),
							column(8, plotOutput(outputId = "jnmx2Plot"))),
						br()),
					tabPanel(title = "Slope Confidence Band Tables",
						htmlOutput(outputId = "mx1cbTitle"),
						htmlOutput(outputId = "jnmx1Text"),
						DT::dataTableOutput("jnmx1cbInfo"),
						htmlOutput(outputId = "mx2cbTitle"),
						htmlOutput(outputId = "jnmx2Text"),
						DT::dataTableOutput("jnmx2cbInfo"))),
				tabPanel(title = "Model Output",
					htmlOutput(outputId = "modFit"),
					inlineCSS(list("table" = "font-size: 10px"))),
				tabPanel(title = "Descriptives", 
					htmlOutput(outputId = "desc.stdTitle"),
					tableOutput(outputId = "desc.std"),
					htmlOutput(outputId = "descTitle"),
					tableOutput(outputId = "desc")),
				tabPanel(title = "Raw Data", 
					DT::dataTableOutput(outputId = "rawdata"),
					htmlOutput(outputId = "rawDataUnavail")),
				tabPanel(title = "Code for Saving Plots", 
					htmlOutput(outputId = "codetitle1")),
				tabPanel(title = "Educational Materials", 
					htmlOutput(outputId = "edu"))))
	)
)

server <- function(input, output, session) {
	output$guide <- renderUI(
		HTML(paste(
			"<div style='font-size: 12px'> 
			<h5> Welcome </h5>
					This visualization tool represents continuous by continuous interactions as a regression plane in 
					3D space in combination with the regions of significance from the Johnson-Neyman technique. You can see an example of the type of plot produced <a href=http://rpubs.com/sbu_mfinsaas/Example target=_blank> here</a>. 
					Please contact Megan Finsaas at megan.finsaas@gmail.com if you have questions or comments about the app. <br> <br>
					For the best performance, please run the app in <b> Google Chrome</b>. <br> <br>
					The plot is powered by Plotly (Sievert, 2018) and the Johnson-Neyman analysis is conducted using the <i> interactions </i> package (Long, 2019). The marginal effects plots are also produced using this package. <br> <br> 
					
					
			<h5> Uploading a Datafile (Raw Data Upload Method) </h5>
			<ul>	
					<li> File requirements: </li>
				    <ul>
							<li>Must have unique variable names at top of the file. Since these names will appear on the plot and the output, we recommend keeping them short (< 8 characters). Variable names that are too long can push the model summary off screen. </li>
							<li>No missing data indicators allowed, but blanks are okay. The program handles missing data using listwise deletion on the variables in the model. </li>
							<li>Must be a comma-separated (.csv) or SPSS (.sav) file. Sometimes these file types may appear grayed out/unavailable in the file upload window. 
								To access them, change the setting in the file upload window to view all file types. </li>
						<li> Variables must be numeric. If using an SPSS file, the program may read variables with few integer values or with value labels as factors. Try removing the value labels and/or using a .csv file instead. </li>
					</ul>
						<li> Once the file has uploaded, the dropdowns will populate with the variable names from the file. </li>
						<li> If the file contains a large number of variables, not all of their names will appear in the dropdown list, but they can still be accessed by beginning to type the variable names. You may wish to subset your dataset to include only the variables in your model prior to uploading the file to mitigate this (and to speed up processing time in general). </li>
</ul>					
<h5> Inputting Variables for Regression Model (Raw Data Upload Method) </h5>
					<ul>
						<li> Click the dropdown buttons for variable names or delete “--” and type the variable name. </li>
						<li> Once the predictors and an outcome are input, click “Run/Update”.
						<li> Predictors and the outcome are centered by default. You can alternatively standardize them or keep them in the raw form. The model will automatically re-estimate whenever the standardization option is changed. </li>
						<li> To add covariates, click the text box for a dropdown list or type the variable names. Covariates are always centered at their means. They can be entered at any time and the model will automatically update to include them.</li>
						<li> At this point, the program does not accept squared terms. </li>			
					</ul>
<h5> Manual Entry Method </h5>
					<ul>
						<li> Ensure that all variables are centered as desired (all covariates must be centered or standardized). Using this input method assumes the interaction is significant.  </li>
						<li> Enter the variable names and minimums and maximums of both predictors. </li>
						<li> Enter the coefficient estimates for the intercept (b0), both predictors (bX1 & bX2), and interaction (bX1*bX2). </li>
						<li> Next enter the asymptotic variances for each coefficient, the covariances for each predictor with the interaction, and the degrees of freedom. Click “Run/Update” to see the plot and other output.  </li>
						<li> This input method is in beta as bug testing is ongoing. Please email with any questions or concerns. </li>
					</ul>
			<h5> 3D Plot Tab </h5>
					<ul>
						<li> Along with the  <b>3D plot</b>, when	a region of significance falls within the observed data, the <b>Johnson-Neyman output</b> will include the values
							on the moderator that mark the bounds of the region of significance, the range of slope estimates for the 
							predictor-outcome relationship and the percentage of cases within this region, and the observed range of the moderator variable. </li>
						<li> When a <b>crossover point</b> falls within the observed data, it will also print on this tab. </li>
					</ul>
			<h5> False Discovery Rate Adjustment for Johnson-Neyman </h5>
					<ul>
						<li> Checking this box adjusts critical t-value for the Johnson-Neyman test according to the paper by Esarey and Sumner (2018), titled, “Marginal Effects in Interaction Models: Determining and Controlling the False Positive Rate.”  </li>
						<li> The adjusted critical t-value is printed next to the marginal effects plots under the <b> Johnson-Neyman </b> tab. </li>
					</ul>
			<h5> Plot Features </h5>
					<ul> 
						<li> The Plotly buttons at the top right of the plot allow you to <b> zoom, pan, rotate, and save a .png </b> of the plot.  </li>
						<li> For all models, checkboxes for <b> adding the scatterplot, regression plane, </b>and <b>95% confidence interval</b> around the predicted values (for raw data upload option only) will appear to the right of the plot. </li>
						<li>	Which other plot feature options appear below these depend on whether there are regions of significance or crossover points within the observed range of data. </li>
							<ul>
								<li> If a region of significance falls within the range of data, the option to <b> shade the region of significance </b> on the regression plane will appear. Selecting the “Solid ROS” button will shade the region of significance on the plot in one color; selecting the “Gradient ROS” button will shade the region using a gradient that reflects the width of the 95% confidence band around the slope estimate.
								Finally, selecting the “Gradient All” option will shade the entire regression plane according to the confidence bands. You can also view the confidence band estimates under the Johnson-Neyman tab. </li>
								<li> If a crossover point falls within the range of data, the “Crossover” checkbox will appear, which <b> adds the crossover point </b>to the figure.</li>
							</ul>
						</ul>
			<h5> Additional Output </h5>
					<ul>
						<li> <b> Johnson-Neyman Tab: </b> Marginal effects plots and Johnson-Neyman bounds regardless of whether they fall within the range of observed data; tables with slope estimates, 95% confidence bands and <i>p</i>-values at all values of the moderator.</li> 
						<li> <b> Model Output Tab: </b> Regression model output. </li>
						<li> <b> Descriptives Tab: </b> Table of descriptives for variables in model with standardization options applied; table for all variables in the dataset in raw form </li>
						<li> <b> Raw Data Tab: </b> Sortable data frame of raw data. </li>
					</ul>
			<h5> Saving and Sharing Plots </h5>
					There are three ways to save and share your plots:
					<ol>
						<li> Take screen/snapshots of the figure in various rotations (i.e., PrtScn on Windows; Shift + Command + 4 on Macs) or click the Plotly camera button to save a .png file. </li>
						<li> Record the plot while manually rotating the plot. To keep the mouse outside the plot while recording, click and drag above it. </li>
						<li> To create a “live” rotatable plot with a permanent link (this option is only available for the raw data upload input method):
						<ul> 
							<li>Download <a href=https://www.r-project.org target=_blank>R</a> and <a href=https://www.rstudio.com target=_blank>RStudio</a>, which are both free programs. </li>
							<li> Run the function code linked in the <b> Code for Saving Plots </b> tab in RStudio. </li>
							<li> Run the input code linked in the <b> Code for Saving Plots </b> tab in RStudio after setting the arguments to your preference. </li>
							<li> Click the blue “Publish” icon in the RStudio viewer and follow the steps to publish to the free service <a href=http://rpubs.com target=_blank>Rpubs</a>. </li>
							<li> This approach works because RStudio supports online publication directly from its viewer window, whereas this action is not supported directly from Shiny apps. </li>
						</ul>
					</ol>
			<h5> Trouble Seeing Plot? </h5>
			If you run into trouble viewing the plot (e.g., the text boxes overlap; can't see the plot fully), please: 
			<ul>
				<li> Double check that you're using Google Chrome. </li>
				<li> Expand your browser window fully. </li>
				<li> Zoom out in your browser window (this option is usually located under View; on Macs, the shortcut is Command -). </li>
				<li> Use the Plotly options to zoom, pan, and rotate the plot. These options appear at the top right side of the plot when you hover your mouse over the plot. Once you have clicked on zoom, pan, or rotate, click on the plot and hold while moving the mouse side to side or up and down to change the view. </li>
				<li> Make sure your variable names are < 8 characters. </li>
			</ul>
			If you continue to have problems, please email me at megan.finsaas@gmail.com.<br> <br> <br> <br>"))
	)
	
	# render UI depending on input method selection ---------------------------
	
	# for manual entry
	
	output$yManNameUI <- renderUI({
		req(input$inputType == "manual")
		textInput(inputId = "yManName", label = "Y (Outcome) Name")
	})
	
	# names of variables
	output$x1ManNameUI <- renderUI({
		req(input$inputType == "manual")
		textInput(inputId = "x1ManName", label = "X1 (Predictor) Name")
	})
	
	output$x2ManNameUI <- renderUI({
		req(input$inputType == "manual")
		textInput(inputId = "x2ManName", label = "X2 (Predictor) Name")
	})
	
	# min/max input
	output$minx1UI <- renderUI({
		req(input$inputType == "manual")
		numericInput(inputId = "minx1", value = 0.00, label = "X1 Min")#
	})
	output$maxx1UI <- renderUI({
		req(input$inputType == "manual")
		numericInput(inputId = "maxx1", value = 0.00, label = "X1 Max") #paste0("X1 Max (",strtrim(x1.name(), 3), ")"))
	})
	
	output$line1 <- renderText({
		req(input$inputType == "manual")
		HTML("<hr>")
	})
	
	output$minx2UI <- renderUI({
		req(input$inputType == "manual")
		numericInput(inputId = "minx2", value = 0.00, label = "X2 Min")
	})
	output$maxx2UI <- renderUI({
		req(input$inputType == "manual")
		numericInput(inputId = "maxx2", value = 0.00, label = "X2 Max")
	})
	
	output$line2 <- renderText({
		req(input$inputType == "manual")
		HTML("<hr>")
	})
	
	# coefficient input
	output$coefTitle <- renderText({
		req(input$inputType == "manual")
		HTML(paste0("<div style='font-weight: 700' >Coefficients</div>"))
	})
	output$b0UI <- renderUI({
		req(input$inputType == "manual")
		numericInput(inputId = "b0", value = 0.00, label = "b0")
	})
	output$bx1UI <- renderUI({
		req(input$inputType == "manual")
		numericInput(inputId = "bx1", value = 0.00, label = "bX1") #paste0("bX1 (",strtrim(x1.name(), 3), ")"))
	})
	output$bx2UI <- renderUI({
		req(input$inputType == "manual")
		numericInput(inputId = "bx2", value = 0.00, label = "bX2") #paste0("bX2 (",strtrim(x2.name(), 3), ")"))
	})
	output$bIntUI <- renderUI({
		req(input$inputType == "manual")
		numericInput(inputId = "bInt", value = 0.00, label = "bX1*X2")
	})
	output$line3 <- renderText({
		req(input$inputType == "manual")
		HTML("<hr>")
	})
	# coefficient variance input
	output$coefVarTitle <- renderText({
		req(input$inputType == "manual")
		HTML(paste0("<div style='font-weight: 700' >Coefficient Variances</div>"))
	})
	output$b0varUI <- renderUI({
		req(input$inputType == "manual")
		numericInput(inputId = "b0var", value = 0.00, label = "b0 Variance")
	})
	output$bx1varUI <- renderUI({
		req(input$inputType == "manual")
		numericInput(inputId = "bx1var", value = 0.00, label = "bX1 Variance") # paste0("b1 (",strtrim(x1.name(), 3), ")"))
	})
	output$bx2varUI <- renderUI({
		req(input$inputType == "manual")
		numericInput(inputId = "bx2var", value = 0.00, label = "bX2 Variance") #paste0("b2 (",strtrim(x2.name(), 3), ")"))
	})
	output$bIntvarUI <- renderUI({
		req(input$inputType == "manual")
		numericInput(inputId = "bIntvar", value = 0.00, label = "bX1*X2 Variance")
	})
	output$line4 <- renderText({
		req(input$inputType == "manual")
		HTML("<hr>")
	})
	# coefficient covariance input
	output$coefCovarTitle <- renderText({
		req(input$inputType == "manual")
		HTML(paste0("<div style='font-weight: 700' >Coefficient Covariances</div>"))
	})
	output$bx1IntcovUI <- renderUI({
		req(input$inputType == "manual")
		numericInput(inputId = "bx1Intcov", value = 0.00, label = "bX1, bX1*X2 Covariance")
	})
	output$bx2IntcovUI <- renderUI({
		req(input$inputType == "manual")
		numericInput(inputId = "bx2Intcov", value = 0.00, label = "bX2, bX1*X2 Covariance")
	})
	
	# df
	output$dfsUI <- renderUI({
		req(input$inputType == "manual") 
		numericInput(inputId = "dfs", value = 0, label = "df (N - # predictors - 1)")
	})
	output$fdrManualUI <- renderUI({
		req(input$inputType == "manual")
		checkboxInput(inputId = "fdrManualCheckbox",
			label = "Limit false discovery rate for Johnson-Neyman?",
			value = checkboxes$fdrManualCB)
	})
	# for raw data
	
	# upload file
	output$uploadFileUI <- renderUI({
		req(input$inputType == "rawData")
		fileInput(inputId = "uploadFile", # input: select a file
			label = "Choose File", multiple = FALSE,
			accept = c(".csv", ".sav"))
	})
	output$line6 <- renderText({
		req(input$inputType == "rawData")
		HTML("<hr>")
	})
	output$fdrRawDataUI <- renderUI({
		req(input$inputType == "rawData")
		checkboxInput(inputId = "fdrRawDataCheckbox",
			label = "Limit false discovery rate for Johnson-Neyman?",
			value = checkboxes$fdrRawDataCB)
	})
	output$line7 <- renderText({
		req(input$inputType == "rawData")
		HTML("<hr>")
	})
	output$yUI <- renderUI({
		req(input$inputType == "rawData")
		selectInput(inputId = "y", 
			label = "Select Y (Dependent) Variable", choices = "--")
	})
	output$std.yUI <- renderUI({
		req(input$inputType == "rawData")
		radioButtons(inputId = "std.y",
			label = NULL,
			choices = c(Center = "center", Standardize = "standardize", Raw = "raw"),
			selected = "center", inline = TRUE)
	})
	output$line8 <- renderText({
		req(input$inputType == "rawData")
		HTML("<hr>")
	})
	output$x1UI <- renderUI({
		req(input$inputType == "rawData")
		selectInput(inputId = "x1", 
			label = "Select X1 (Predictor) Variable", choices = "--")
	})
	output$std.x1UI <- renderUI({
		req(input$inputType == "rawData")
		radioButtons(inputId = "std.x1",
			label = NULL,
			choices = c(Center = "center", Standardize = "standardize", Raw = "raw"),
			selected = "center", inline = TRUE)
	})
	output$line9 <- renderText({
		req(input$inputType == "rawData")
		HTML("<hr>")
	})
	output$x2UI <- renderUI({
		req(input$inputType == "rawData")
		selectInput(inputId = "x2", 
			label = "Select X2 (Predictor) Variable", choices = "--")
	})
	output$std.x2UI <- renderUI({
		req(input$inputType == "rawData")
		radioButtons(inputId = "std.x2",
			label = NULL,
			choices = c(Center = "center", Standardize = "standardize", Raw = "raw"),
			selected = "center", inline = TRUE)
	})
	output$line10 <- renderText({
		req(input$inputType == "rawData")
		HTML("<hr>")
	})
	output$covUI <- renderUI({
		req(input$inputType == "rawData")
		selectInput(inputId = "cov", 
			label = "Select Covariate(s) (If none, leave blank)", choices = "--", multiple = TRUE)
	})
	output$line11 <- renderText({
		req(input$inputType == "rawData")
		HTML("<hr>")
	})
	output$yPrintNameUI <- renderUI({
		req(input$inputType == "rawData")
		textInput(inputId = "yPrintName", label = "Name for Display")
	})
	output$x1PrintNameUI <- renderUI({
		req(input$inputType == "rawData")
		textInput(inputId = "x1PrintName", label = "Name for Display")
	})
	output$x2PrintNameUI <- renderUI({
		req(input$inputType == "rawData")
		textInput(inputId = "x2PrintName", label = "Name for Display")
	})
	
	#load file
	df <- reactive({
		req(input$uploadFile)
		if(tools::file_ext(input$uploadFile) == "sav") {
			foreign::read.spss(input$uploadFile$datapath, to.data.frame = TRUE)
		} else {
			if(tools::file_ext(input$uploadFile) == "csv") {
				read.csv(input$uploadFile$datapath, header = TRUE)
			}
		}
	})
	
	# update variable drop downs in ui
	observeEvent(df(), {
		updateSelectInput(session, "y", choices = c("--", names(df())))
		updateSelectInput(session, "x1", choices = c("--", names(df())))
		updateSelectInput(session, "x2", choices = c("--", names(df())))
		updateSelectInput(session, "cov", choices = c(names(df())))
	})
	
	# display raw data
	output$rawdata <- renderDataTable({
		if (input$inputType == "rawData") {
			req(input$uploadFile)
			DT::datatable(data = df(),
				options = list(pageLength = 10, lengthMenu = c(10, 25, 40)),
				rownames = FALSE)
		}
	})
	
	output$rawDataUnavail <- renderUI({
		if (input$inputType == "manual") {
			"This feature is only available when the 'upload raw data' input option is used."
		}
	})
	
	# test conditions ---------------------------------------------------------
	
	# manual input
	# all0 <- eventReactive(input$run, {
	# 	req(input$inputType == "manual")
	# 	if (input$minx1 == 0 & input$maxx1 == 0 & input$minx2 == 0 & input$maxx2 == 0 &
	# 			input$b0 == 0 & input$bx1 == 0 & input$bx2 == 0 & input$bInt == 0 &
	# 			input$b0var == 0 & input$bx1var == 0 & input$bx2var == 0 & input$bIntvar == 0 &
	# 			input$bx1Intcov == 0 & input$bx2Intcov == 0 & input$dfs == 0) {
	# 		all0 <- TRUE
	# 	} else {
	# 		all0 <- FALSE
	# 	}
	# })
	# 
	
	int0 <- eventReactive(input$run, {
		req(input$inputType == "manual")
		if (input$bInt == 0) {
			int0 <- TRUE
		} else {
			int0 <- FALSE
		}
	})
	
	anyNA <- eventReactive(input$run, {
		req(input$inputType == "manual")
		if ((input$minx1 == "" | input$maxx1 == "" | input$minx2 == "" | input$maxx2 == "" |
				input$b0 == "" | input$bx1 == "" | input$bx2 == "" | input$bInt == "" &
				input$b0var == "" | input$bx1var == "" | input$bx2var == "" | input$bIntvar == "" |
				input$bx1Intcov == "" | input$bx2Intcov == "" | input$dfs == "") |
				!is.numeric(input$minx1) | !is.numeric(input$maxx1) | !is.numeric(input$minx2) | !is.numeric(input$maxx2) |
				!is.numeric(input$b0) | !is.numeric(input$bx1) | !is.numeric(input$bx2) | !is.numeric(input$bInt) &
				!is.numeric(input$b0var) | !is.numeric(input$bx1var) | !is.numeric(input$bx2var) | !is.numeric(input$bIntvar) |
				!is.numeric(input$bx1Intcov) | !is.numeric(input$bx2Intcov) | !is.numeric(input$dfs)) {
			anyNA <- TRUE
		} else {
			anyNA <- FALSE
		}
	})
	
	minGEMaxX1 <- eventReactive(input$run, {
		req(input$inputType == "manual") 
		if(input$minx1 >= input$maxx1) {
			minGEMaxX1 <- TRUE
		} else {
			minGEMaxX1 <- FALSE
		}
	})
	
	minGEMaxX2 <- eventReactive(input$run, {
		req(input$inputType == "manual") 
		if(input$minx2 >= input$maxx2) {
			minGEMaxX2 <- TRUE
		} else {
			minGEMaxX2 <- FALSE
		}
	})
	
	negVar <- eventReactive(input$run, {
		req(input$inputType == "manual") 
		if(input$b0var < 0 | input$bx1var < 0 | input$bx2var < 0 | input$bIntvar < 0) {
			negVar <- TRUE
		} else {
			negVar <- FALSE
		}
	})
	
	notRight <- eventReactive(input$run, {
		req(input$inputType == "manual")
		bx1var <- input$bx1var
		bIntvar <- input$bIntvar
		bx2var <- input$bx2var
		bx2Intcov <- input$bx2Intcov
		bx1Intcov <- input$bx1Intcov
		vecx1w <- seq(minx1(), maxx1(), length.out = 100)
		vecx2w <- seq(minx2(), maxx2(), length.out = 100)
		if (sum(sign(bx2var + 2 * vecx1w * bx2Intcov + 	vecx1w^2 * bIntvar)) < length(vecx1w) |  # x1 as mod
				sum(sign(bx1var + 2 * vecx2w * bx1Intcov + 	vecx2w^2 * bIntvar)) < length(vecx2w)) {  #x1 as mod
			notRight <- TRUE
		} else {
			notRight <- FALSE
		}
	})
	
	# raw data input
	# same variable input twice 
	sameVar <- reactive({
		req(input$inputType == "rawData", input$x1 != "--", input$x2 != "--", input$y != "--")
		if(is.null(input$cov)) {
			if(input$x1 == input$x2 | input$x1 == input$y | input$x2 == input$y) {
				sameVar <- TRUE
			} else {
				sameVar <- FALSE
			}
		} else {
			if(length(input$cov) == 1) {
				if(input$x1 == input$x2 | input$x1 == input$y | input$x2 == input$y | input$cov == input$x1 |
						input$cov == input$x2 | input$cov == input$y) {
					sameVar <- TRUE
				} else {
					sameVar <- FALSE
				}
			} else {
				if(length(input$cov) > 1) {
					if(input$x1 == input$x2 | input$x1 == input$y | input$x2 == input$y | any(input$cov == input$x1) |
							any(input$cov == input$x2) | any(input$cov == input$y)) {
						sameVar <- TRUE
					} else {
						sameVar <- FALSE
					}
				}
			}
		}
	})
	
	# variables not numeric
	non.numeric <- reactive({
		req(input$inputType == "rawData", input$x1 != "--", input$x2 != "--", input$y != "--")
		mat <- as.matrix(df()[, c(input$y, input$x1, input$x2, input$cov)])
		if(is.numeric(mat) == TRUE) {
			non.numeric <- FALSE
		} else {
			non.numeric <- TRUE
		}
	})
	

# move to 3D plot ---------------------------------------------------------
	
	observeEvent(input$run, {
		if (input$inputType == "rawData") {
			req(input$x1 != "--", input$x2 != "--", input$y != "--", non.numeric() == FALSE, sameVar() == FALSE)
		} else {
			if(input$inputType == "manual") {			
				req(input$dfs > 1, anyNA() == FALSE, int0() == FALSE, notRight() == FALSE, minGEMaxX1() == FALSE, minGEMaxX2 () == FALSE, negVar() == FALSE)
			}
		}
		updateTabsetPanel(session, "inTabset",
			selected = "3DPlotPanel")
	})
	
	
	# render warnings/messages for sidebar ------------
	messageSameVar <- eventReactive(input$run, {
		req(input$inputType == "rawData", input$x1 != "--", input$x2 != "--", input$y != "--")
		if (sameVar() == TRUE) {
			HTML(paste0("<br/> <div style='background-color: #ff9999; border-radius: 5px;
				text-align: center; border-width: 1px; border-style: solid; border-color: #808080; padding:5px;
				font-size: 10px; font-color: #808080'> Check your model. The same variable was input twice. </div>"))
		} 
	})
	
	output$messageSameVar <- renderText({
		messageSameVar()
	})
	
	messageNonNum <- eventReactive(input$run, {
		req(input$inputType == "rawData", input$x1 != "--", input$x2 != "--", input$y != "--")
		if (non.numeric() == TRUE) {
			HTML(paste0("<br/> <div style='background-color: #ff9999; border-radius: 5px;
			text-align: center; border-width: 1px; border-style: solid; border-color: #808080; padding:5px;
			font-size: 12px; font-color: #808080'> Check your variables. At least one seems to be non-numeric. If
		  they are all numeric, it may be that one variable has only a few integer values or (if using an SPSS file) has value labels and is being read as a factor. Try 
			using a different file type and/or removing the value labels. </div>"))
		} 
	})
	
	output$messageNonNum <- renderText({
		messageNonNum()
	})
	
	messageInt <- eventReactive(input$run, {
		if (input$inputType == "rawData") {
			req(input$x1 != "--", input$x2 != "--", input$y != "--", non.numeric() == FALSE, sameVar() == FALSE)
			if (intsig() == FALSE) {
				HTML(paste0("<br/> <div style='background-color: #ff9999; border-radius: 5px;
					text-align: center; border-width: 1px; border-style: solid; border-color:#808080;
					font-size: 10px; font-color: #808080; padding:10px;' >", 
					paste0("The interaction between ", x1.name(), " and ", x2.name(), " does not significantly predict ",
						y.name(), ". To view all 3D plot features and the Johnson-Neyman test results, enter a model with a significant ",
						"interaction term. You can still see the regression plane and scatterplot on the <b> 3D Plot </b> tab."), "</div>"))
			} else {
				HTML(paste0("<br/> <div style='background-color: #d9d9d9; border-radius: 5px;
					text-align: center; border-width: 1px; border-style: solid; border-color: #808080; padding:8px;
					font-size: 10px; font-color: #808080'> Click the <b> 3D Plot </b> tab to see the 3D plot and 
					Johnson-Neyman output for your model. </div>"))
			}
		} else { 
			req(input$dfs > 1, anyNA() == FALSE, int0() == FALSE, notRight() == FALSE, minGEMaxX1() == FALSE, minGEMaxX2 () == FALSE, negVar() == FALSE)
			HTML(paste0("<br/> <div style='background-color: #d9d9d9; border-radius: 5px;
					text-align: center; border-width: 1px; border-style: solid; border-color: #808080; padding:8px;
					font-size: 10px; font-color: #808080'> Click the <b> 3D Plot </b> tab to see the 3D plot and 
					Johnson-Neyman output for your model. </div>"))
		}
	})
	
	output$messageInt <- renderText({
		messageInt()
	})
	
	# manual input 
	
	messageMinGEMaxX1 <- eventReactive(input$run, {
		req(input$inputType == "manual", input$minx1 != "", input$maxx1 != "") 
		if(minGEMaxX1() == TRUE) {
			HTML(paste0("<div style='background-color: #ff9999; border-radius: 5px;
					text-align: center; border-width: 1px; border-style: solid; border-color:#808080;
					font-size: 10px; font-color: #808080; padding:10px;' >", 	paste0("The minimum of ", x1.name(), 
						" cannot be equal to or greater than its maximum. Please check your entries."), "</div>"))
		}
	})
	
	output$messageMinGEMaxX1 <- renderText({
		messageMinGEMaxX1()
	})
	
	
	messageMinGEMaxX2 <- eventReactive(input$run, {
		req(input$inputType == "manual", input$minx2 != "", input$maxx2 != "") 
		if(minGEMaxX2() == TRUE) {
			HTML(paste0("<br/> <div style='background-color: #ff9999; border-radius: 5px;
					text-align: center; border-width: 1px; border-style: solid; border-color:#808080;
					font-size: 10px; font-color: #808080; padding:10px;' >", 
				paste0("The minimum of ", x2.name(), " cannot be equal to or greater than its maximum. Please check your entries."), "</div>"))
		}
	})
	
	output$messageMinGEMaxX2 <- renderText({
		messageMinGEMaxX2()
	})
	
	
	messageNegVar <- eventReactive(input$run, {
		req(input$inputType == "manual", input$b0var != "", input$bx1var != "", input$bx2var != "", input$bIntvar !="")
		if(negVar() == TRUE) {
			HTML(paste0("<br/> <div style='background-color: #ff9999; border-radius: 5px;
					text-align: center; border-width: 1px; border-style: solid; border-color:#808080;
					font-size: 10px; font-color: #808080; padding:10px;' >", 
				paste0("The asymptotic variances of the coefficients cannot be negative, and the asymptotic variances of the x1, x2, and interaction coefficents cannot equal zero. Please check your entries."), "</div>"))
		}
	})
	
	output$messageNegVar <- renderText({
		messageNegVar()
	})
	
	messageNotRight <- eventReactive(input$run, {
		req(input$inputType == "manual", anyNA() == FALSE, input$dfs > 1, minGEMaxX1() == FALSE, minGEMaxX1() == FALSE, negVar() == FALSE) 
		if (notRight() == TRUE) { 
			HTML(paste0("<br/> <div style='background-color: #ff9999; border-radius: 5px;
					text-align: center; border-width: 1px; border-style: solid; border-color:#808080;
					font-size: 10px; font-color: #808080; padding:10px;' >", 
				paste0("Something is not right. Please check your entries."), "</div>"))
		}
	})
	
	output$messageNotRight <- renderText({
		messageNotRight()
	})
	
	messageOneNA <- eventReactive(input$run, {
		req(input$inputType == "manual") 
		if (anyNA() == TRUE) {  
			HTML(paste0("<br/> <div style='background-color: #ff9999; border-radius: 5px;
					text-align: center; border-width: 1px; border-style: solid; border-color:#808080;
					font-size: 10px; font-color: #808080; padding:10px;' >", 
				paste0("At least one input is blank or a non-plausible value. Please check your entries."), "</div>"))
			# }
		}
	})
	
	output$messageOneNA <- renderText({
		messageOneNA()
	})
	
	messageDf <- eventReactive(input$run, {
		req(input$inputType == "manual") 
		if (input$dfs <= 1) {  
			HTML(paste0("<br/> <div style='background-color: #ff9999; border-radius: 5px;
					text-align: center; border-width: 1px; border-style: solid; border-color:#808080;
					font-size: 10px; font-color: #808080; padding:10px;' >", 
				paste0("Degrees of freedom must be larger than 1. Please check your entries."), "</div>"))
		}
	})
	
	output$messageDf <- renderText({
		messageDf()
	})
	
	messageInt0 <- eventReactive(input$run, {
		req(input$inputType == "manual", !is.na(input$bInt)) 
		if (int0() == TRUE) {  
			HTML(paste0("<br/> <div style='background-color: #ff9999; border-radius: 5px;
					text-align: center; border-width: 1px; border-style: solid; border-color:#808080;
					font-size: 10px; font-color: #808080; padding:10px;' >", 
				paste0("The interaction coefficient cannot be zero. Please check your entries."), "</div>"))
		}
	})
	
	output$messageInt0 <- renderText({
		messageInt0()
	})
	
	# standardize variables and subset data --------
	df.std <- eventReactive(input$run, {
		req(input$inputType == "rawData", input$x1 != '--', input$x2 != '--', input$y != '--', sameVar() == FALSE, non.numeric() == FALSE)
		data <- df() 
		
		# keep only cases with values on all variables
		data	<-	data[ , c(input$y, input$x1, input$x2, input$cov)]
		data <- data[complete.cases(data), ]
		
		
		if (input$std.x1 == "standardize") {
			data[,input$x1] <- c(scale(data[,input$x1]))
		}  else {
			if (input$std.x1 == "raw") {
				data[,input$x1] <- data[,input$x1]
			} else {
				if (input$std.x1 == "center") {
					data[,input$x1] <- c(scale(data[,input$x1], scale = FALSE))
				}
			}
		}
		if (input$std.x2 == "standardize") {
			data[,input$x2] <- c(scale(data[,input$x2]))
		}  else {
			if (input$std.x2 == "raw") {
				data[,input$x2] <- data[,input$x2]
			} else {
				if (input$std.x2 == "center") {
					data[,input$x2] <- c(scale(data[,input$x2], scale = FALSE))
				}
			}
		}
		if (input$std.y == "standardize") {
			data[,input$y] <- c(scale(data[,input$y]))
		}  else {
			if (input$std.y == "raw") {
				data[,input$y] <- data[,input$y]
			} else {
				if (input$std.y == "center") {
					data[,input$y] <- c(scale(data[,input$y], scale = FALSE))
				}
			}
		}
		data[,input$cov] <- c(scale(data[,input$cov], scale = FALSE)) # covariates always centered
		data[, c(input$y, input$x1, input$x2, input$cov)] 
	})
	
	# descriptives ------------------------------------------------------------
	output$desc <- renderTable({ # descriptives on raw data
		req(input$inputType == "rawData", input$uploadFile)
		descdf <- describe(df())[c("n","mean","sd","median","min","max","skew","kurtosis")]
		descdf <- as.data.frame(descdf)
		Variable <- colnames(df())
		cbind(Variable, descdf)
	})
	
	output$desc.std <- renderTable({ # descriptives on variables in model with standardization option applied
		req(input$inputType == "rawData")
		descdf <- describe(df.std())[c("n","mean","sd","median","min","max","skew","kurtosis")]
		descdf <- as.data.frame(descdf)
		Variable <- colnames(df.std())
		if (input$x1PrintName != "" | input$x2PrintName != "" | input$yPrintName != "") {
			Display <- c(y.name(), x1.name(), x2.name(),  
				rep("", length.out = length(coef(mod())) - 4))
			descdf <- cbind(Variable, Display, descdf)
			names(descdf) <- c("Variable", "Display Name", "n","mean","sd","median","min","max","skew","kurtosis")
			descdf
		} else {
			cbind(Variable, descdf)
		}
	})
	
	# render descriptives output
	output$descTitle <- renderUI({
		if (input$inputType == "rawData") {
			req(input$uploadFile)
			h6("Descriptive Statistics on Raw Data")
		} else {
			"This feature is only available when the 'upload raw data' input option is used."
		}
	})
	
	output$desc.stdTitle <- renderUI({
		if (input$inputType == "rawData") {
			req(input$uploadFile, input$x1 != '--', input$x2 != '--', input$y != '--', sameVar() == FALSE, non.numeric() == FALSE)
			h6("Descriptive Statistics on Variables in Model with Standardization Options Applied")
		} else {
			NULL
		}
	})
	
	
	# estimate linear model ---------------------------------------------------
	mod <- eventReactive(input$run, {
		lm(paste(input$y, " ~ . + ", input$x1, " * ", input$x2), data = df.std())
	})
	
	modsum <- eventReactive(input$run, {
		summary(mod())
	})
	
	# model output ------------------------------------------------------------
	# output$modInfo <- renderUI({
	# 	req(input$x1 != '--', input$x2 != '--', input$y != '--', sameVar() == FALSE, non.numeric() == FALSE)
	# 	HTML(paste("<div style='font-weight: 700' >Model Summary</div>"))
	# })
	# 
	modFit <- eventReactive(input$run, {
		if (input$inputType == "rawData") {
			req(input$x1 != '--', input$x2 != '--', input$y != '--', sameVar() == FALSE, non.numeric() == FALSE)
			if (input$yPrintName != "") {
				dv <- paste0("Dependent variable: ", input$y, " (", y.name(), ")")
			} else {
				dv <- paste0("Dependent variable: ", input$y)
			}
			
			rsq <- modsum()$r.squared
			rsq <- format(round(rsq, digits = 3), nsmall = 3)
			rsq <- paste0("R-squared = ", rsq)
			arsq <- modsum()$adj.r.squared
			arsq <- format(round(arsq, digits = 3), nsmall = 3)
			arsq <- paste0("Adjusted R-squared = ", arsq)
			
			f <- modsum()$fstatistic[1]
			fdf1 <- modsum()$fstatistic[2]
			fdf2 <- modsum()$fstatistic[3]
			fp <- 1 - stats::pf(f, fdf1, fdf2)
			f <- format(round(f, digits = 3), nsmall = 3)
			fp <- format(round(fp, digits = 3), nsmall = 3)
			fsum <- paste0("<i>F</i>(",fdf1, ", ",fdf2, ") = ", f, ", <i>p</i> = ", fp)
			
			nobs <- length(coef(mod())) + modsum()$fstatistic[3]
			nobs <- paste0("Number of observations: ", nobs)
			# 
			HTML(paste("<div class='well' > <label> Model Summary </label> <br/>", modTable(), dv, "<br/>", nobs, "<br/>", fsum, "<br/>", rsq, "; ", arsq, "</div>"))
		} else {
			"This feature is only available when the 'upload raw data' input option is used."
		}
	})
	
	modTable <- renderTable({
		req(input$x1 != '--', input$x2 != '--', input$y != '--', sameVar() == FALSE, non.numeric() == FALSE)
		modtable <- unlist(modsum()[4])
		modtable <- data.frame(matrix(modtable, ncol = 4))
		modtable <- format(round(modtable, digits = 3), nsmall = 3)
		varNames <- names(coef(mod()))
		if (input$x1PrintName != "" | input$x2PrintName != "") {
			dispName <- c("", paste0(" (", x1.name(), ")"), paste0(" (", x2.name(), ")"), 
				rep("", length.out = length(coef(mod())) - 4), paste0(" (", x1.name(), " X ", x2.name(), ")"))
			varNames <- paste0(varNames, dispName)
			#		modtable <- cbind(, dispName, modtable)
			#	colnames(modtable) <- c(" ", "Est.","SE", "t-value","p-value")
			# colnames(modtable) <- c(" ", "Display Name", "Est.","SE", "t-value","p-value")
		} 
		modtable <- cbind(varNames, modtable)
		colnames(modtable) <- c(" ", "Est.","SE", "t-value","p-value")
		
		return(modtable)
	})
	
	output$modFit <- renderUI({
		modFit()
	})
	# test whether interaction term is significant ----------------------------
	intsig <- eventReactive(input$run, {
		if (input$inputType == "rawData") {
			modsum <- summary(mod())
			modsum$coef[nrow(modsum$coef), 4] < .05 # access significance column for interaction term 
		} else {
			req(int0() == FALSE)
			intsig <- TRUE
		}
	})
	
	
	# identify crossover points -----------------------------------------------
	
	# x1 as mod: -Bx1/Bx1*x2, plot on x2
	coonx2 <- eventReactive(input$run, {
		if (input$inputType == "manual") {
			req(input$dfs > 1, negVar() == FALSE, minGEMaxX1() == FALSE, minGEMaxX2() == FALSE, anyNA() == FALSE)
		}
		coonx2 <- -bx1()/bInt()
		if(coonx2 >= minx2() & coonx2 <= maxx2()) {
			coonx2 <- as.numeric(coonx2)
		} else {
			coonx2 <- NA
		}
	})
	
	# predicted y at crossover
	coonx2y <- eventReactive(input$run, {
		if (!is.na(coonx2())) {
			b0() + bx1()*vecx1() + bx2()*coonx2() + bInt()*vecx1()*coonx2()
		} else {
			coonx2y <- NA
		}
	})
	
	# x2 as mod: -Bx2/Bx1*x2, plot on x1
	coonx1 <- eventReactive(input$run, {
		if (input$inputType == "manual") {
			req(input$dfs > 1, negVar() == FALSE, minGEMaxX1() == FALSE, minGEMaxX2() == FALSE, anyNA() == FALSE)
		}
		coonx1 <- -bx2()/bInt()
		if(coonx1 >= minx1() & coonx1 <= maxx1()) {
			coonx1 <- as.numeric(coonx1)
		} else {
			coonx1 <- NA
		}
	})
	
	# predicted y at crossover
	coonx1y <- eventReactive(input$run, {
		if (!is.na(coonx1())) {
			b0() + bx1()*coonx1() + bx2()*vecx2() + bInt()*coonx1()*vecx2()
		} else {
			coonx1y <- NA
		}
	})
	
	# save coefficients -------------------------------------------------------
	
	b0 <- eventReactive(input$run, {
		if (input$inputType == "rawData") {
			b0 <- coef(mod())[1]
		} else {
			b0 <- input$b0
		}
	})
	
	bx1 <- eventReactive(input$run, {
		if (input$inputType == "rawData") {
			mod <- mod()
			bx1 <- coef(mod)[2]
		} else {
			bx1 <- input$bx1
		}
	})
	
	bx2 <- eventReactive(input$run, {
		if (input$inputType == "rawData") {
			bx2 <- coef(mod())[3]
		} else {
			bx2 <- input$bx2
		}
	})
	
	bInt <- eventReactive(input$run, {
		if (input$inputType == "rawData") {
			bInt <- coef(mod())[length(coef(mod()))]
		} else {
			req(input$bInt != 0)
			bInt <- input$bInt
		}
	})
	
	
	# min/max of x1/x2 --------------------------------------------------------
	minx1 <- eventReactive(input$run, {
		if (input$inputType == "rawData") {
			data <- df.std()
			min(data[,input$x1], na.rm = TRUE)
		} else {
			req(anyNA() == FALSE, input$dfs > 1, minGEMaxX1() == FALSE, minGEMaxX2() == FALSE, negVar() == FALSE) 
			minx1 <- input$minx1
		}
	})
	
	maxx1 <- eventReactive(input$run, {
		if (input$inputType == "rawData") {
			data <- df.std()
			max(data[,input$x1], na.rm = TRUE)
		} else {
			req(anyNA() == FALSE, input$dfs > 1, minGEMaxX1() == FALSE, minGEMaxX2() == FALSE, negVar() == FALSE) 
			maxx1 <- input$maxx1
		}
	})
	
	minx2 <- eventReactive(input$run, {
		if (input$inputType == "rawData") {
			data <- df.std()
			min(data[,input$x2], na.rm = TRUE)
		} else {
			req(anyNA() == FALSE, input$dfs > 1, minGEMaxX1() == FALSE, minGEMaxX2() == FALSE, negVar() == FALSE) 
			minx2 <- input$minx2
		}
	})
	
	maxx2 <- eventReactive(input$run, {
		if (input$inputType == "rawData") {
			data <- df.std()
			max(data[,input$x2], na.rm = TRUE)
		} else {
			req(anyNA() == FALSE, input$dfs > 1, minGEMaxX1() == FALSE, minGEMaxX2() == FALSE, negVar() == FALSE) 
			maxx2 <- input$maxx2
		}
	})
	
	
	# grayscale option --------------------------------------------------------
	
	x1modcol <- eventReactive(input$run, {
		if (input$grayscale == FALSE) {
			x1modcol <- "#567001"
		} else {
			x1modcol <- "#404040"
		}
	})	
	x2modcol <- eventReactive(input$run, {
		if (input$grayscale == FALSE) {
			x2modcol <- "#3366ff"
		} else {
			x2modcol <- "#404040"
		}
	})	
	
	# save inputs for JN functions -------------------------------------------------------
	
	fdr <- eventReactive(input$run, {
		if (input$inputType == "rawData") {
			fdr <- input$fdrRawDataCheckbox
		} else {
			fdr <- input$fdrManualCheckbox
		}
	})
	
	names.coef.manual <- eventReactive(input$run, {
		if (input$inputType == "rawData") {
			names(coef(mod()))
		} else {
			names.coef.manual <- c("(Intercept)", x1.name(), x2.name(), 
				paste0(x1.name(),":", x2.name()))
		}
	})
	
	vcov.manual <- eventReactive(input$run, {
		if (input$inputType == "rawData") {
			vcov(mod())	
		} else {
			b0var <- input$b0var
			bx1var <- input$bx1var
			bIntvar <- input$bIntvar
			bx2var <- input$bx2var
			bx2Intcov <- input$bx2Intcov
			bx1Intcov <- input$bx1Intcov
			bIntvar <- input$bIntvar
			vcov.manual <- matrix(c(
				b0var, 0,         0,         0, 
				0,     bx1var,    0,         bx1Intcov,
				0,     0,         bx2var,    bx2Intcov,
				0,     bx1Intcov, bx2Intcov, bIntvar), ncol=4)
			rownames(vcov.manual) <- names.coef.manual()
			colnames(vcov.manual) <- names.coef.manual()
			vcov.manual
		}
	})
	
	df.residual.manual <- eventReactive(input$run, {
		if (input$inputType == "rawData") {
			df.residual(mod())	
		} else {
			df.residual.manual <- input$dfs
			df.residual.manual
		}
	})
	
	coef.model.manual <- eventReactive(input$run, {
		if (input$inputType == "rawData") {
			coef(mod())
		} else {
			coef.model.manual <- c(b0(), bx1(), bx2(), bInt())
			names(coef.model.manual) <- names.coef.manual()
			coef.model.manual
		}
	})
	
	
	# save variable names -----------------------------------------------------
	
	y.name <- eventReactive(input$run, {
		if (input$inputType == "rawData") {
			if(input$yPrintName == "") {
				input$y
			} else {
				input$yPrintName
			}
		} else {
			if(input$yManName == "") {
				"Y"
			} else {
				input$yManName
			}
		}
	})
	
	x1.name <- eventReactive(input$run, {
		if (input$inputType == "rawData") {
			if(input$x1PrintName == "") {
				input$x1
			} else {
				input$x1PrintName
			}
		} else {
			if(input$x1ManName == "") {
				"X1"
			} else {
				input$x1ManName
			}
		}
	})
	
	x2.name <- eventReactive(input$run, {
		if (input$inputType == "rawData") {
			if(input$x2PrintName == "") {
				input$x2
			} else {
				input$x2PrintName
			}
		} else {
			if(input$x2ManName == "") {
				"X2" 
			} else {
				input$x2ManName
			}
		}
	})
	
	jnx1.name <- eventReactive(input$run, {
		if (input$inputType == "rawData") {
			input$x1
		} else {
			x1.name()
		}
	})
	
	jnx2.name <- eventReactive(input$run, {
		if (input$inputType == "rawData") {
			input$x2
		} else {
			x2.name()
		}
	})
	
	
	# johnson neyman analysis -------------------------------------------------
	# when x1 is mod
	jnmx1 <- eventReactive(input$run, {
		if (input$inputType == "manual") {
			req(anyNA() == FALSE, notRight() == FALSE, input$dfs > 1, negVar() == FALSE, minGEMaxX1() == FALSE, minGEMaxX2() == FALSE)
		}
		johnson_neyman_new(pred = jnx2.name(), modx = jnx1.name(), vmat = vcov.manual(), alpha = 0.05, plot = TRUE, 
			control.fdr = fdr(), line.thickness = 0.5, df = "residual", 
			digits = getOption("jtools-digits", 2), critical.t = NULL, 
			sig.color = x1modcol(), insig.color = "#b3b2b2", mod.range = c(minx1(), maxx1()), #"#749702"
			title = "Marginal Effects Plot",
			df.residual.mod = df.residual.manual(),
			coef.model = coef.model.manual(),
			names.coef.model = names.coef.manual()) 
	})
	
	# save jn values if within mod range
	jnmx1.1 <- eventReactive(input$run, {
		if(jnmx1()$bounds[1] >= minx1() & jnmx1()$bounds[1] <= maxx1()) {
			jnmx1.1 <- as.numeric(jnmx1()$bounds[1])
		} else {
			jnmx1.1 <- NA
		}
		jnmx1.1
	})
	
	jnmx1.2 <- eventReactive(input$run, {
		if(jnmx1()$bounds[2] >= minx1() & jnmx1()$bounds[2] <= maxx1()) {
			jnmx1.2 <- as.numeric(jnmx1()$bounds[2])
		} else {
			jnmx1.2 <- NA
		}
	})
	
	# when x2 is mod
	jnmx2 <- eventReactive(input$run, {
		if (input$inputType == "manual") {
			req(anyNA() == FALSE, notRight() == FALSE, input$dfs > 1, negVar() == FALSE, minGEMaxX1() == FALSE, minGEMaxX2() == FALSE)
		}
		johnson_neyman_new(pred = jnx1.name(), modx = jnx2.name(), vmat = vcov.manual(), alpha = 0.05, plot = TRUE, 
			control.fdr = fdr(), line.thickness = 0.5, df = "residual", 
			digits = getOption("jtools-digits", 2), critical.t = NULL, 
			sig.color = x2modcol(), insig.color = "#b3b2b2", mod.range = c(minx2(), maxx2()), 
			title = "Marginal Effects Plot",
			df.residual.mod = df.residual.manual(),
			coef.model = coef.model.manual(),
			names.coef.model = names.coef.manual()) # added arguments 
	})
	
	# save jn values if within mod range
	jnmx2.1 <- eventReactive(input$run, {
		if(jnmx2()$bounds[1] >= minx2() & jnmx2()$bounds[1] <= maxx2()) {
			jnmx2.1 <- as.numeric(jnmx2()$bounds[1])
		} else {
			jnmx2.1 <- NA
		}
	})
	
	jnmx2.2 <- eventReactive(input$run, {
		if(jnmx2()$bounds[2] >= minx2() & jnmx2()$bounds[2] <= maxx2()) {
			jnmx2.2 <- as.numeric(jnmx2()$bounds[2])
		} else {
			jnmx2.2 <- NA
		}
	})
	
	
	
	# jn confidence band tables -----------------------------------------------
	jnmx1cb <- eventReactive(input$run, {
		johnson_neyman_veccb2(pred = jnx2.name(), modx = jnx1.name(), vmat = vcov.manual(), alpha = 0.05, plot = FALSE, 
			control.fdr = fdr(), line.thickness = 0.5, df = "residual", 
			digits = getOption("jtools-digits", 2), critical.t = NULL, 
			mod.range = c(minx1(), maxx1()), 
			title = "Johnson-Neyman plot",
			df.residual.mod = df.residual.manual(),
			coef.model = coef.model.manual(),
			names.coef.model = names.coef.manual(), vec = vecx1())$cbands
	})
	
	jnmx2cb <- eventReactive(input$run, {
		johnson_neyman_veccb2(pred = jnx1.name(), modx = jnx2.name(), vmat = vcov.manual(), alpha = 0.05, plot = FALSE, 
			control.fdr = fdr(), line.thickness = 0.5, df = "residual", 
			digits = getOption("jtools-digits", 2), critical.t = NULL, 
			mod.range = c(minx2(), maxx2()), 
			title = "Johnson-Neyman plot",
			df.residual.mod = df.residual.manual(),
			coef.model = coef.model.manual(),
			names.coef.model = names.coef.manual(), vec = vecx2())$cbands
	})
	
	# johnson neyman output ---------------------------------------------------
	# marginal effects plots
	jnmx1Title <- eventReactive(input$run, {
		if (input$inputType == "rawData") {
			req(input$x1 != '--', input$x2 != '--', input$y != '--', sameVar() == FALSE, non.numeric() == FALSE)
		} else {
			req(anyNA() == FALSE, notRight() == FALSE, input$dfs > 1, negVar() == FALSE, minGEMaxX1() == FALSE, minGEMaxX2() == FALSE)
		}
		h5(paste0(x1.name(), " as Moderator"))
	})
	
	output$jnmx1Title <- renderUI({
		jnmx1Title()
	})
	
	jnmx1Plot <- eventReactive(input$run, {
		if (input$inputType == "rawData") {
			req(input$x1 != '--', input$x2 != '--', input$y != '--', sameVar() == FALSE, non.numeric() == FALSE)
		} else {
			req(anyNA() == FALSE, notRight() == FALSE, input$dfs > 1, negVar() == FALSE, minGEMaxX1() == FALSE, minGEMaxX2() == FALSE)
		}
		if(intsig() == TRUE & mx1.ROS.none() == FALSE) {
			jnmx1()$plot
		} 
	})
	output$jnmx1Plot <- renderPlot({
		jnmx1Plot()
	})
	
	jnmx1Info <- eventReactive(input$run, {
		if (input$inputType == "rawData") {
			req(input$x1 != '--', input$x2 != '--', input$y != '--', sameVar() == FALSE, non.numeric() == FALSE)
		} else {
			req(anyNA() == FALSE, notRight() == FALSE, input$dfs > 1, negVar() == FALSE, minGEMaxX1() == FALSE, minGEMaxX2() == FALSE)
		}
		if(intsig() == TRUE & mx1.ROS.none() == FALSE) {
			jnmx1()
		} else {
			HTML(paste0("The interaction term is not significant and/or there is no region of significance within the range of ", x1.name(),"."))
		}
	})
	
	output$jnmx1Info <- renderPrint({
		jnmx1Info()
	})
	
	jnmx2Title <- eventReactive(input$run, {
		if (input$inputType == "rawData") {
			req(input$x1 != '--', input$x2 != '--', input$y != '--', sameVar() == FALSE, non.numeric() == FALSE)
		} else {
			req(anyNA() == FALSE, notRight() == FALSE, input$dfs > 1, negVar() == FALSE, minGEMaxX1() == FALSE, minGEMaxX2() == FALSE)
		}
		h5(paste0(x2.name(), " as Moderator"))
	})
	
	output$jnmx2Title <- renderUI({
		jnmx2Title()
	})
	
	jnmx2Plot <- eventReactive(input$run, {
		if (input$inputType == "rawData") {
			req(input$x1 != '--', input$x2 != '--', input$y != '--', sameVar() == FALSE, non.numeric() == FALSE)
		} else {
			req(anyNA() == FALSE, notRight() == FALSE, input$dfs > 1, negVar() == FALSE, minGEMaxX1() == FALSE, minGEMaxX2() == FALSE)
		}
		if(intsig() == TRUE & mx2.ROS.none() == FALSE) {
			jnmx2()$plot
		}
	})
	
	output$jnmx2Plot <- renderPlot({
		jnmx2Plot()
	})
	
	jnmx2Info <- eventReactive(input$run, {
		if (input$inputType == "rawData") {
			req(input$x1 != '--', input$x2 != '--', input$y != '--', sameVar() == FALSE, non.numeric() == FALSE)
		} else {
			req(anyNA() == FALSE, notRight() == FALSE, input$dfs > 1, negVar() == FALSE, minGEMaxX1() == FALSE, minGEMaxX2() == FALSE)
		}
		if(intsig() == TRUE & mx2.ROS.none() == FALSE) {
			jnmx2()
		} else {
			HTML(paste0("The interaction term is not significant and/or there is no region of significance within the range of ", x2.name(),"."))
		}
	})
	
	output$jnmx2Info <- renderPrint({
		jnmx2Info()
	})
	
	
	# confidence band tables
	mx1cbTitle <- eventReactive(input$run, {
		if (input$inputType == "rawData") {
			req(input$x1 != '--', input$x2 != '--', input$y != '--', sameVar() == FALSE, non.numeric() == FALSE)
		} else {
			req(anyNA() == FALSE, notRight() == FALSE, input$dfs > 1, negVar() == FALSE, minGEMaxX1() == FALSE, minGEMaxX2() == FALSE)
		}
		h5(paste0(x1.name(), " as Moderator"))
	})
	
	output$mx1cbTitle <- renderUI({
		mx1cbTitle()
	})
	
	jnmx1cbInfo <- eventReactive(input$run, {
		if (input$inputType == "rawData") {
			req(input$x1 != '--', input$x2 != '--', input$y != '--', sameVar() == FALSE, non.numeric() == FALSE)
		} else {
			req(anyNA() == FALSE, notRight() == FALSE, input$dfs > 1, negVar() == FALSE, minGEMaxX1() == FALSE, minGEMaxX2() == FALSE)
		}
		if(intsig() == TRUE & mx1.ROS.none() == FALSE) {
			jnmx1cb <- jnmx1cb()
			jnmx1cb[,1:4] <- format(round(jnmx1cb[,1:4], digits = 3), nsmall = 3)
			DT::datatable(jnmx1cb, options = list(lengthMenu = c(10,20,50, length(vecx1()))))
		}
	})
	
	output$jnmx1cbInfo <- DT::renderDataTable({
		jnmx1cbInfo()
	})
	
	jnmx1Text <- eventReactive(input$run, {
		if (input$inputType == "rawData") {
			req(input$x1 != '--', input$x2 != '--', input$y != '--', sameVar() == FALSE, non.numeric() == FALSE)
		} else {
			req(anyNA() == FALSE, notRight() == FALSE, input$dfs > 1, negVar() == FALSE, minGEMaxX1() == FALSE, minGEMaxX2() == FALSE)
		}
		if(intsig() == FALSE | mx1.ROS.none() == TRUE) {
			HTML(paste0("The interaction term is not significant and/or there is no region of significance within the range of ", x1.name(),"."))
		}
	})
	
	output$jnmx1Text <- renderPrint({
		jnmx1Text()
	})
	
	mx2cbTitle <- eventReactive(input$run, {
		if (input$inputType == "rawData") {
			req(input$x1 != '--', input$x2 != '--', input$y != '--', sameVar() == FALSE, non.numeric() == FALSE)
		} else {
			req(anyNA() == FALSE, notRight() == FALSE, input$dfs > 1, negVar() == FALSE, minGEMaxX1() == FALSE, minGEMaxX2() == FALSE)
		}
		h5(paste0(x2.name(), " as Moderator"))
	})
	
	output$mx2cbTitle <- renderUI({
		mx2cbTitle()
	})
	
	jnmx2cbInfo <- eventReactive(input$run, {
		if (input$inputType == "rawData") {
			req(input$x1 != '--', input$x2 != '--', input$y != '--', sameVar() == FALSE, non.numeric() == FALSE)
		} else {
			req(anyNA() == FALSE, notRight() == FALSE, input$dfs > 1, negVar() == FALSE, minGEMaxX1() == FALSE, minGEMaxX2() == FALSE)
		}
		if(intsig() == TRUE & mx2.ROS.none() == FALSE) {
			jnmx2cb <- jnmx2cb()
			jnmx2cb[,1:4] <- format(round(jnmx2cb[,1:4], digits = 3), nsmall = 3)
			DT::datatable(jnmx2cb, options = list(lengthMenu = c(10,20,50, length(vecx2()))))
		}
	})
	
	output$jnmx2cbInfo <- DT::renderDataTable({
		jnmx2cbInfo()
	})
	
	jnmx2Text <- eventReactive(input$run, {
		if (input$inputType == "rawData") {
			req(input$x1 != '--', input$x2 != '--', input$y != '--', sameVar() == FALSE, non.numeric() == FALSE)
		} else {
			req(anyNA() == FALSE, notRight() == FALSE, input$dfs > 1, negVar() == FALSE, minGEMaxX1() == FALSE, minGEMaxX2() == FALSE)
		}
		if(intsig() == FALSE | mx2.ROS.none() == TRUE) {
			HTML(paste0("The interaction term is not significant and/or there is no region of significance within the range of ", x2.name(),"."))
		}
	})
	
	output$jnmx2Text <- renderPrint({
		jnmx2Text()
	})
	
	# johnson neyman and crossover text output for main tab ---------------------------------------------------------
	jnmx1.sum <- eventReactive(input$run, {
		if (input$inputType == "rawData") {
			req(input$x1 != '--', input$x2 != '--', input$y != '--', sameVar() == FALSE, non.numeric() == FALSE)
		} else {
			req(anyNA() == FALSE, notRight() == FALSE, input$dfs > 1, negVar() == FALSE, minGEMaxX1() == FALSE, minGEMaxX2() == FALSE)
		}
		if(intsig() == TRUE & mx1.ROS.none() == FALSE) {
			# get sign of slope at jn bounds
			if (mx1.ROS.all() == TRUE) {
				jnmx1.1sl <- sign(jnmx1cb()[1, 2])
				jnmx1.2sl <- sign(jnmx1cb()[nrow(jnmx1cb()), 2])
			} else {
				jnmx1.1sl <- sign(jnmx1cb()[jnmx1.1rn(), 2])
				jnmx1.2sl <- sign(jnmx1cb()[jnmx1.2rn(), 2])
			}
			jnmx1.1sl <- ifelse(jnmx1.1sl > 0, "positive", "negative")
			jnmx1.2sl <- ifelse(jnmx1.2sl > 0, "positive", "negative")
			
			# set values for text printing
			jnmx1.sl <- c(jnmx1.1sl, jnmx1.2sl) # pos or neg slope
			jnmx1.sign <- c(" >= ",  " <= ")
			minmaxx1 <- c(jnmx1cb()[1, 2], jnmx1cb()[nrow(jnmx1cb()), 2]) # min max of mod
			jnmx1.val <- c(jnmx1cb()[jnmx1.1rn(), 1], jnmx1cb()[jnmx1.2rn(), 1]) # jn bound values
			jnmx1.slest <- c(jnmx1cb()[jnmx1.1rn(), 2], jnmx1cb()[jnmx1.2rn(), 2]) # slope at jn bounds
			minmaxx1 <- c(jnmx1cb()[1, 1], jnmx1cb()[nrow(jnmx1cb()), 1]) # min max of mod
			jnmx1.slestminmax <- c(jnmx1cb()[1, 2], jnmx1cb()[nrow(jnmx1cb()), 2]) # slope at min max of mod
			obsorjnmx1 <- c(" (minimum)", " (maximum)", " (JN bound)")
			
			# formatting
			jnmx1.val <- format(round(jnmx1.val, digits = 2), nsmall = 2)
			jnmx1.slest <- format(round(jnmx1.slest, digits = 2), nsmall = 2)
			minmaxx1 <- format(round(minmaxx1, digits = 2), nsmall = 2)
			jnmx1.slestminmax <- format(round(jnmx1.slestminmax, digits = 2), nsmall = 2)
			
			# jn text summary
			if(mx1.ROS.all() == TRUE) { # ROS covers whole range of mod
				jnmx1.sum <- list(jnmx1.sl[1], jnmx1.slestminmax[1], jnmx1.slestminmax[2], minmaxx1[1], minmaxx1[2],
					obsorjnmx1[1], obsorjnmx1[2])
			} else {
				if(jnmx1.inside() == FALSE) { # ROS falls outside JN values
					if(is.na(jnmx1.1()) & !is.na(jnmx1.2())) { # high band / outside
						if (input$inputType == "rawData") {
							jnmx1.prop <- prop.table(table(df.std()[,input$x1] >= jnmx1.2()))
							jnmx1.prop <- jnmx1.prop[dim(jnmx1.prop)]
						} else {
							jnmx1.prop <- NULL
						}
						jnmx1.sum <- list(jnmx1.sl[2], jnmx1.slest[2], jnmx1.slestminmax[2], jnmx1.sign[1], jnmx1.val[2],
							jnmx1.sign[2], minmaxx1[2], obsorjnmx1[3], obsorjnmx1[2])
					} else {
						if(!is.na(jnmx1.1()) & is.na(jnmx1.2())) { # low band / outside
							if (input$inputType == "rawData") {
								jnmx1.prop <- prop.table(table(df.std()[,input$x1] <= jnmx1.1()))
								jnmx1.prop <- jnmx1.prop[dim(jnmx1.prop)]
							} else {
								jnmx1.prop <- NULL
							}
							jnmx1.sum <- list(jnmx1.sl[1], jnmx1.slestminmax[1], jnmx1.slest[1], jnmx1.sign[1], minmaxx1[1],
								jnmx1.sign[2], jnmx1.val[1], obsorjnmx1[1], obsorjnmx1[3])
						} else { # high and low band / outside
							if(!is.na(jnmx1.1()) & !is.na(jnmx1.2())) {
								if (input$inputType == "rawData") {
									jnmx1.prop <- prop.table(table(df.std()[,input$x1] >= jnmx1.1() & 
											df.std()[,input$x1] <= jnmx1.2()))
									jnmx1.prop <- jnmx1.prop[dim(jnmx1.prop)]
								} else {
									jnmx1.prop <- NULL
								}
								jnmx1.sum <- list(jnmx1.sl[1], jnmx1.sign[1], minmaxx1[1], jnmx1.sign[2], jnmx1.val[1], 
									jnmx1.slestminmax[1], jnmx1.slest[1], 
									jnmx1.sl[2], jnmx1.sign[1], jnmx1.val[2], jnmx1.sign[2], minmaxx1[2], 
									jnmx1.slest[2], jnmx1.slestminmax[2], obsorjnmx1[1], obsorjnmx1[3], obsorjnmx1[3], obsorjnmx1[2])
							}
						}
					}
				} else { # ROS falls inside JN values
					if(!is.na(jnmx1.1()) & !is.na(jnmx1.2())) { 
						if (input$inputType == "rawData") {
							jnmx1.prop <- prop.table(table(df.std()[,input$x1] >= jnmx1.1() & 
									df.std()[,input$x1] <= jnmx1.2()))
							jnmx1.prop <- jnmx1.prop[dim(jnmx1.prop)]
						} else {
							jnmx1.prop <- NULL
						}
						jnmx1.sum <- list(jnmx1.sl, jnmx1.slest[1], jnmx1.slest[2], 
							jnmx1.sign[1], jnmx1.val[1], jnmx1.sign[2], minmaxx1[2], obsorjnmx1[3], obsorjnmx1[3])
					} else {
						if(is.na(jnmx1.1()) & !is.na(jnmx1.2())) { # low band / inside
							if (input$inputType == "rawData") {
								jnmx1.prop <- prop.table(table(df.std()[,input$x1] <= jnmx1.2()))
								jnmx1.prop <- jnmx1.prop[dim(jnmx1.prop)]
							} else {
								jnmx1.prop <- NULL
							}
							jnmx1.sum <- list(jnmx1.sl[2], jnmx1.slestminmax[1], jnmx1.slest[2], jnmx1.sign[1], minmaxx1[1],
								jnmx1.sign[2], jnmx1.val[2], obsorjnmx1[1], obsorjnmx1[3])
						} else {
							if (!is.na(jnmx1.1()) & is.na(jnmx1.2())) { # high band / inside
								if (input$inputType == "rawData") {
									jnmx1.prop <- prop.table(table(df.std()[,input$x1] >= jnmx1.1()))
									jnmx1.prop <- jnmx1.prop[dim(jnmx1.prop)]
								} else {
									jnmx1.prop <- NULL
								}
								jnmx1.sum <- list(jnmx1.sl[1], jnmx1.slest[1], jnmx1.slestminmax[2], jnmx1.sign[1], jnmx1.val[1], 
									jnmx1.sign[2], minmaxx1[2], obsorjnmx1[3], obsorjnmx1[2])
							} 
						}
					}
				}
			}
			
			if(mx1.ROS.all() == TRUE) {
				if (input$inputType == "rawData") {
					jnmx1.perc <- paste0("This region of significance includes 100.00% of the sample. ")
				} else {
					jnmx1.perc <- NULL
				}
				jnmx1.text <- paste0("The relationship between ", x2.name(), " and ", y.name(),
					" is significantly ", jnmx1.sum[1], " at all values of ", x1.name(), " (minimum: ",
					jnmx1.sum[4], "; maximum: ", jnmx1.sum[5], "; slope estimate range: ", jnmx1.sum[2], " to " , jnmx1.sum[3],
					"). ", jnmx1.perc)
			} else {
				if (jnmx1.inside() == FALSE & !is.na(jnmx1.1()) & !is.na(jnmx1.2())) { # two ROS (high and low) 
					if (input$inputType == "rawData") {
						jnmx1.perc <- jnmx1.prop*100
						jnmx1.perc <- format(round(jnmx1.perc, digits = 2), nsmall = 2)
						jnmx1.perc <- paste0("These regions of significance include ", jnmx1.perc, "% of the sample. ")
					} else {
						jnmx1.perc <- NULL
					}
					jnmx1.text <- paste0("The relationship between ", x2.name(), " and ", y.name(),
						" is significantly ", jnmx1.sum[1], " when ", x1.name(), " is",
						jnmx1.sum[2], jnmx1.sum[3], jnmx1.sum[15], " and", jnmx1.sum[4], jnmx1.sum[5], jnmx1.sum[16],
						" (slope estimate range: ", jnmx1.sum[6], " to ", jnmx1.sum[7], ") and significantly ", jnmx1.sum[8], " when ",
						x1.name(), " is", jnmx1.sum[9], jnmx1.sum[10],jnmx1.sum[17], " and ", jnmx1.sum[11], jnmx1.sum[12], jnmx1.sum[18],
						" (slope estimate range: ", jnmx1.sum[13], " to ", jnmx1.sum[14], 
						"). ", jnmx1.perc)
				} else {
					if (jnmx1.inside() == TRUE & !is.na(jnmx1.1()) & !is.na(jnmx1.2())) { # one center ROS 
						if (input$inputType == "rawData") {
							jnmx1.perc <- jnmx1.prop*100
							jnmx1.perc <- format(round(jnmx1.perc, digits = 2), nsmall = 2)
							jnmx1.perc <- paste0("These regions of significance include ", jnmx1.perc, "% of the sample. ")
						} else {
							jnmx1.perc <- NULL
						}
						jnmx1.text <- paste0("When ", x1.name(), " is between ",	jnmx1.val[1], " (JN bound) and ", jnmx1.val[2], " (JN bound), ", 
							x2.name(), " and ", y.name(), " are significantly related (slope estimate range: ", jnmx1.slest[1], 
							" to ", jnmx1.slest[2], "). ", jnmx1.perc)
					} else { # one high or low ROS
						if (input$inputType == "rawData") {
							jnmx1.perc <- jnmx1.prop*100
							jnmx1.perc <- format(round(jnmx1.perc, digits = 2), nsmall = 2)
							jnmx1.perc <- paste0("This region of significance includes ", jnmx1.perc, "% of the sample. ")
						} else {
							jnmx1.perc <- NULL
						}
						jnmx1.text <- paste0("When ", x1.name(), " is",	jnmx1.sum[4], jnmx1.sum[5], jnmx1.sum[8], " and", jnmx1.sum[6], jnmx1.sum[7], jnmx1.sum[9],
							", the relationship between ", x2.name(), " and ", y.name(),
							" is significantly ", jnmx1.sum[1], " (slope estimate range: ", jnmx1.sum[2], 
							" to ", jnmx1.sum[3], "). ", jnmx1.perc)
					}
				}
			}
		} else {
			jnmx1.text <- paste0("The relationship between ", x2.name(), " and ", y.name(), " is not significantly different from zero ",
				"at any value of ", x1.name(), ".")
		}
		
		if (!is.na(coonx1())) {
			coonx1 <- format(round(coonx1(), digits = 2), nsmall = 2)
			coonx1.text <- paste0("Crossover point: ", coonx1)
		} else {
			coonx1.text <- paste0("The crossover point falls outside the range of ", x1.name(), ".")
		}
		
		HTML(paste0( "<br/>",
			jnmx1.text, "<br/> <br/>", coonx1.text, "</div>"))
	})
	
	output$jnmx1.sum <- renderUI({
		jnmx1.sum()
	})
	
	jnmx2.sum <- eventReactive(input$run, {
		if (input$inputType == "rawData") {
			req(input$x1 != '--', input$x2 != '--', input$y != '--', sameVar() == FALSE, non.numeric() == FALSE)
		} else {
			req(anyNA() == FALSE, notRight() == FALSE, input$dfs > 1, negVar() == FALSE, minGEMaxX1() == FALSE, minGEMaxX2() == FALSE)
		}
		if(intsig() == TRUE & mx2.ROS.none() == FALSE) {
			
			# get sign of slope at jn bounds
			if (mx2.ROS.all() == TRUE) {
				jnmx2.1sl <- sign(jnmx2cb()[1, 2])
				jnmx2.2sl <- sign(jnmx2cb()[nrow(jnmx2cb()), 2])
			} else {
				jnmx2.1sl <- sign(jnmx2cb()[jnmx2.1cn(), 2])
				jnmx2.2sl <- sign(jnmx2cb()[jnmx2.2cn(), 2])
			}
			jnmx2.1sl <- ifelse(jnmx2.1sl > 0, "positive", "negative")
			jnmx2.2sl <- ifelse(jnmx2.2sl > 0, "positive", "negative")
			
			# set values for text printing 
			jnmx2.sl <- c(jnmx2.1sl, jnmx2.2sl) # pos or neg slope at jn values
			jnmx2.sign <- c(" >= ",  " <= ")
			minmaxx2 <- c(jnmx2cb()[1, 2], jnmx2cb()[nrow(jnmx2cb()), 2]) # min max of mod
			jnmx2.val <- c(jnmx2cb()[jnmx2.1cn(), 1], jnmx2cb()[jnmx2.2cn(), 1]) # jn bound values
			jnmx2.slest <- c(jnmx2cb()[jnmx2.1cn(), 2], jnmx2cb()[jnmx2.2cn(), 2]) # slope at jn bounds
			minmaxx2 <- c(jnmx2cb()[1, 1], jnmx2cb()[nrow(jnmx2cb()), 1]) # min max of mod
			jnmx2.slestminmax <- c(jnmx2cb()[1, 2], jnmx2cb()[nrow(jnmx2cb()), 2]) # slope at min max of mod
			obsorjnmx2 <- c(" (minimum)", " (maximum)", " (JN bound)")
			
			# formatting
			jnmx2.val <- format(round(jnmx2.val, digits = 2), nsmall = 2)
			jnmx2.slest <- format(round(jnmx2.slest, digits = 2), nsmall = 2)
			minmaxx2 <- format(round(minmaxx2, digits = 2), nsmall = 2)
			jnmx2.slestminmax <- format(round(jnmx2.slestminmax, digits = 2), nsmall = 2)
			
			# jn text summary
			if(mx2.ROS.all() == TRUE) { # ROS covers whole range of mod
				jnmx2.sum <- list(jnmx2.sl[1], jnmx2.slestminmax[1], jnmx2.slestminmax[2], minmaxx2[1], minmaxx2[2],
					obsorjnmx2[1], obsorjnmx2[2])
			} else {
				if(jnmx2.inside() == FALSE) { # ROS falls outside JN values
					if(is.na(jnmx2.1()) & !is.na(jnmx2.2())) { # high band / outside
						if(input$inputType == "rawData") {
							jnmx2.prop <- prop.table(table(df.std()[,input$x2] >= jnmx2.2()))
							jnmx2.prop <- jnmx2.prop[dim(jnmx2.prop)]
						} else {
							jnmx2.prop <- NULL
						}
						jnmx2.sum <- list(jnmx2.sl[2], jnmx2.slest[2], jnmx2.slestminmax[2], jnmx2.sign[1], jnmx2.val[2],
							jnmx2.sign[2], minmaxx2[2], obsorjnmx2[3], obsorjnmx2[2])
					} else {
						if(!is.na(jnmx2.1()) & is.na(jnmx2.2())) { # low band / outside
							if(input$inputType == "rawData") {
								jnmx2.prop <- prop.table(table(df.std()[,input$x2] <= jnmx2.1()))
								jnmx2.prop <- jnmx2.prop[dim(jnmx2.prop)]
							} else {
								jnmx2.prop <- NULL
							}
							jnmx2.sum <- list(jnmx2.sl[1], jnmx2.slestminmax[1], jnmx2.slest[1], jnmx2.sign[1], minmaxx2[1],
								jnmx2.sign[2], jnmx2.val[1], obsorjnmx2[1], obsorjnmx2[3])
						} else { # high and low band / outside
							if(!is.na(jnmx2.1()) & !is.na(jnmx2.2())) {
								if(input$inputType == "rawData") {
									jnmx2.prop <- prop.table(table(df.std()[,input$x2] >= jnmx2.1() & 
											df.std()[,input$x2] <= jnmx2.2()))
									jnmx2.prop <- jnmx2.prop[dim(jnmx2.prop)]
								} else {
									jnmx2.prop <- NULL
								}
								jnmx2.sum <- list(jnmx2.sl[1], jnmx2.sign[1], minmaxx2[1], jnmx2.sign[2], jnmx2.val[1], 
									jnmx2.slestminmax[1], jnmx2.slest[1], 
									jnmx2.sl[2], jnmx2.sign[1], jnmx2.val[2], jnmx2.sign[2], minmaxx2[2], 
									jnmx2.slest[2], jnmx2.slestminmax[2], obsorjnmx2[1], obsorjnmx2[3], obsorjnmx2[3], obsorjnmx2[2])
							}
						}
					}
				} else { # ROS falls inside JN values
					if(!is.na(jnmx2.1()) & !is.na(jnmx2.2())) { # (no squared terms allowed, center band not run)
						if (input$inputType == "rawData") {
							jnmx2.prop <- prop.table(table(df.std()[,input$x2] >= jnmx2.1() & 
									df.std()[,input$x2] <= jnmx2.2()))
							jnmx2.prop <- jnmx2.prop[dim(jnmx2.prop)]
						} else {
							jnmx2.prop <- NULL
						}
						jnmx2.sum <- list(jnmx2.sl[1], jnmx2.slest[1], jnmx2.slest[2], 
							jnmx2.sign[1], jnmx2.val[1], jnmx2.sign[2], minmaxx2[2], obsorjnmx2[3], obsorjnmx2[3])
					} else {
						if(is.na(jnmx2.1()) & !is.na(jnmx2.2())) { # low band / inside
							if (input$inputType == "rawData") {
								jnmx2.prop <- prop.table(table(df.std()[,input$x2] <= jnmx2.2()))
								jnmx2.prop <- jnmx2.prop[dim(jnmx2.prop)]
							} else {
								jnmx2.prop <- NULL
							}
							jnmx2.sum <- list(jnmx2.sl[2], jnmx2.slestminmax[1], jnmx2.slest[2], jnmx2.sign[1], minmaxx2[1],
								jnmx2.sign[2], jnmx2.val[2], obsorjnmx2[1], obsorjnmx2[3])
						} else {
							if (!is.na(jnmx2.1()) & is.na(jnmx2.2())) { # high band / inside
								if (input$inputType == "rawData") {
									jnmx2.prop <- prop.table(table(df.std()[,input$x2] >= jnmx2.1()))
									jnmx2.prop <- jnmx2.prop[dim(jnmx2.prop)]
								} else {
									jnmx2.prop <- NULL
								}
								jnmx2.sum <- list(jnmx2.sl[1], jnmx2.slest[1], jnmx2.slestminmax[2], jnmx2.sign[1], jnmx2.val[1], 
									jnmx2.sign[2], minmaxx2[2], obsorjnmx2[3], obsorjnmx2[2])
							}
						}
					}
				}
			}
			
			if(mx2.ROS.all() == TRUE) {
				if (input$inputType == "rawData") {
					jnmx2.perc <- paste0("This region of significance includes 100.00% of the sample. ")
				} else {
					jnmx2.perc <- NULL
				}
				jnmx2.text <- paste0("The relationship between ", x1.name(), " and ", y.name(),
					" is significantly ", jnmx2.sum[1], " at all values of ", x2.name(), " (minimum: ",
					jnmx2.sum[4], "; maximum: ", jnmx2.sum[5], "; slope estimate range: ", jnmx2.sum[2], " to " , jnmx2.sum[3],
					"). ", jnmx2.perc)
			} else {
				if (jnmx2.inside() == FALSE & !is.na(jnmx2.1()) & !is.na(jnmx2.2())) { # two ROS (high and low) 
					if (input$inputType == "rawData") {
						jnmx2.perc <- jnmx2.prop*100
						jnmx2.perc <- format(round(jnmx2.perc, digits = 2), nsmall = 2)
						jnmx2.perc <- paste0("These regions of significance include ", jnmx2.perc, "% of the sample. ")
					} else {
						jnmx2.perc <- NULL
					}
					jnmx2.text <- paste0("The relationship between ", x1.name(), " and ", y.name(),
						" is significantly ", jnmx2.sum[1], " when ", x2.name(), " is", 
						jnmx2.sum[2], jnmx2.sum[3], jnmx2.sum[15], " and", jnmx2.sum[4], jnmx2.sum[5], jnmx2.sum[16],
						" (slope estimate range: ", jnmx2.sum[6], " to ", jnmx2.sum[7], ") and significantly ", jnmx2.sum[8], " when ",
						x2.name(), " is", jnmx2.sum[9], jnmx2.sum[10],jnmx2.sum[17], " and ", jnmx2.sum[11], jnmx2.sum[12], jnmx2.sum[18],
						" (slope estimate range: ", jnmx2.sum[13], " to ", jnmx2.sum[14], 
						"). ", jnmx2.perc)
				} else {
					if (jnmx2.inside() == TRUE & !is.na(jnmx2.1()) & !is.na(jnmx2.2())) { # one center ROS 
						if (input$inputType == "rawData") {
							jnmx2.perc <- jnmx2.prop*100
							jnmx2.perc <- format(round(jnmx2.perc, digits = 2), nsmall = 2)
							jnmx2.perc <- paste0("These regions of significance include ", jnmx2.perc, "% of the sample. ")
						} else {
							jnmx2.perc <- NULL
						}
						jnmx2.text <- paste0("When ", x2.name(), " is between ",	jnmx2.val[1], " (JN bound) and ", jnmx2.val[2], " (JN bound), ", 
							x1.name(), " and ", y.name(), " are significantly related (slope estimate range: ", jnmx2.slest[1], 
							" to ", jnmx2.slest[2], "). ", jnmx2.perc)
					} else { # one high or low ROS
						if (input$inputType == "rawData") {
							jnmx2.perc <- jnmx2.prop*100
							jnmx2.perc <- format(round(jnmx2.perc, digits = 2), nsmall = 2)
							jnmx2.perc <- paste0("This region of significance includes ", jnmx2.perc, "% of the sample. ")
						} else {
							jnmx2.perc <- NULL
						}
						jnmx2.text <- paste0("When ", x2.name(), " is",	jnmx2.sum[4], jnmx2.sum[5], jnmx2.sum[8], " and", jnmx2.sum[6], jnmx2.sum[7], jnmx2.sum[9],
							", the relationship between ", x1.name(), " and ", y.name(),
							" is significantly ", jnmx2.sum[1], " (slope estimate range: ", jnmx2.sum[2], 
							" to ", jnmx2.sum[3], "). ", jnmx2.perc)
					}
				}
			}
		} else {
			jnmx2.text <- paste0("The relationship between ", x1.name(), " and ", y.name(), " is not significantly different from zero ",
				"at any value of ", x2.name(), ".")
		}
		
		if (!is.na(coonx2())) {
			coonx2 <- format(round(coonx2(), digits = 2), nsmall = 2)
			coonx2.text <- paste0("Crossover point: ", coonx2)
		} else {
			coonx2.text <- paste0("The crossover point falls outside the range of ", x2.name(), ".")
		}
		
		HTML(paste0("<br/>",
			jnmx2.text, "<br/> <br/>", coonx2.text, "</div>"))
	})
	
	output$jnmx2.sum <- renderUI({
		jnmx2.sum()
	})
	
	# vectors of x1 and x2 -----------------------------------
	vecx1 <- eventReactive(input$run, {
		vecx1 <- seq(minx1(), maxx1(), length.out = 100)
		vecx1 <- c(vecx1, jnmx1.1(), jnmx1.2()) # add in JN values 
		vecx1 <- sort(vecx1)
	})
	
	vecx2 <- eventReactive(input$run, {
		vecx2 <- seq(minx2(), maxx2(), length.out = 100)
		vecx2 <- c(vecx2, jnmx2.1(), jnmx2.2())	# add in JN values 
		vecx2 <- sort(vecx2)
	})
	
	# predicted values matrix -------------------------------------------------
	pm <- eventReactive(input$run, {
		pm <- outer(vecx1(), vecx2(), function(x1, x2) {
			b0() + bx1()*x1 + bx2()*x2 + bInt()*x1*x2
		})
		rownames(pm) <- vecx1()
		colnames(pm) <- vecx2()
		pm
	})
	# pm <- eventReactive(input$run, {
	# 	pm <- outer(vecx1(), vecx2(), function(pmx1, pmx2) {
	# 		data <- df.std()
	# 		if (is.null(input$cov)) {
	# 			newdatapm <- data.frame(pmx1, pmx2)
	# 			names(newdatapm) <- c(input$x1, input$x2)
	# 		} else {
	# 			if (length(input$cov) == 1) {
	# 				x <- mean(data[,4])
	# 				pmcov.l <- list(x) 
	# 			} else {
	# 				pmcov.l <- sapply(c(data[,4:ncol(data)]), mean, simplify = FALSE)
	# 			}
	# 			newdatapm <- data.frame(pmx1, pmx2, pmcov.l)
	# 			names(newdatapm) <- c(input$x1, input$x2, input$cov)
	# 		}
	# 		predict.lm(mod(), newdatapm)
	# 	})
	# 	rownames(pm) <- vecx1()
	# 	colnames(pm) <- vecx2()
	# 	pm
	# })
	# 
	# lower ci of predicted values matrix -------------------------------------------------
	pmcil <- eventReactive(input$run, {
		pmcil <- outer(vecx1(), vecx2(), function(pmx1, pmx2) {
			data <- df.std()
			if (is.null(input$cov)) {
				newdatapm <- data.frame(pmx1, pmx2)
				names(newdatapm) <- c(input$x1, input$x2)
			} else {
				if (length(input$cov) == 1) {
					x <- mean(data[,4])
					pmcov.l <- list(x) 
				} else {
					pmcov.l <- sapply(c(data[,4:ncol(data)]), mean, simplify = FALSE)
				}
				newdatapm <- data.frame(pmx1, pmx2, pmcov.l)
				names(newdatapm) <- c(input$x1, input$x2, input$cov)
			}
			predict.lm(mod(), newdatapm, se.fit = TRUE, interval = "confidence",
				level = 0.95)$fit[,"lwr"]
		})
		rownames(pmcil) <- vecx1()
		colnames(pmcil) <- vecx2()
		pmcil
	})
	
	# upper ci of predicted values matrix -------------------------------------------------
	pmciu <- eventReactive(input$run, {
		pmciu <- outer(vecx1(), vecx2(), function(pmx1, pmx2) {
			data <- df.std()
			if (is.null(input$cov)) {
				newdatapm <- data.frame(pmx1, pmx2)
				names(newdatapm) <- c(input$x1, input$x2)
			} else {
				if (length(input$cov) == 1) {
					x <- mean(data[,4])
					pmcov.l <- list(x) 
				} else {
					pmcov.l <- sapply(c(data[,4:ncol(data)]), mean, simplify = FALSE)
				}
				newdatapm <- data.frame(pmx1, pmx2, pmcov.l)
				names(newdatapm) <- c(input$x1, input$x2, input$cov)
			}
			predict.lm(mod(), newdatapm, se.fit = TRUE, interval = "confidence",
				level = 0.95)$fit[,"upr"]
		})
		rownames(pmciu) <- vecx1()
		colnames(pmciu) <- vecx2()
		pmciu
	})
	
	
	# create jn ROS planes ----------------------------------------------------
	
	# for x1 as mod
	# identify whether ROS is inside JN values
	jnmx1.inside <- eventReactive(input$run, {
		jnmx1()$inside
	})
	
	# identify whether ROS covers whole range of x1
	mx1.ROS.all <- eventReactive(input$run, {
		if((jnmx1.inside() == FALSE & (jnmx1()$bounds[2] < minx1() | jnmx1()$bounds[1] > maxx1())) |
				(jnmx1.inside() == TRUE & jnmx1()$bounds[1] < minx1() & jnmx1()$bounds[2] > maxx1())) {
			mx1.ROS.all <- TRUE
			mx1.ROS.all
		} else {
			mx1.ROS.all <- FALSE
			mx1.ROS.all
		}
	})
	
	# identify whether ROS within range of x1
	mx1.ROS.none <- eventReactive(input$run, {
		if((is.na(jnmx1.1()) & is.na(jnmx1.2())) & mx1.ROS.all() == FALSE) {
			mx1.ROS.none <- TRUE
			mx1.ROS.none
		} else {
			mx1.ROS.none <- FALSE
			mx1.ROS.none
		}
	})
	
	# identify row number of jn values
	jnmx1.1rn <- eventReactive(input$run, {
		jnmx1.1rn <- match(jnmx1.1(), vecx1())
	})
	
	jnmx1.2rn <- eventReactive(input$run, {
		jnmx1.2rn <- match(jnmx1.2(), vecx1())
	})
	
	# make x1 ROS planes
	pm.mx1.wROS <- eventReactive(input$run, {
		if (mx1.ROS.none() == FALSE & intsig() == TRUE) {
			
			pm <- pm()
			jnmx1.1 <- jnmx1.1()
			jnmx1.2 <- jnmx1.2()
			vecx1 <- vecx1()
			vecx2 <- vecx2()
			
			# ROS covers whole range
			if(mx1.ROS.all()==TRUE ) {
				pm.mx1.wROS <- pm() 
			}
			
			# ROS inside JN bounds
			if(jnmx1.inside() == TRUE) {
				# inside both JN values (center band) 
				if(!is.na(jnmx1.1) & !is.na(jnmx1.2)) {
					pm.mx1.ROS1i <- pm[jnmx1.1rn():jnmx1.2rn(), ]
					mx1.pad1 <- matrix(nrow = jnmx1.1rn() - 1, ncol = ncol(pm)) # lower pad made of NAs
					mx1.pad2 <- matrix(nrow = nrow(pm) - jnmx1.2rn(), ncol = ncol(pm)) # upper pad  made of NAs
					pm.mx1.wROS <- rbind(mx1.pad1, pm.mx1.ROS1i, mx1.pad2)
				}
				# lower JN bound within range (high band)
				if(!is.na(jnmx1.1) & is.na(jnmx1.2)) {
					pm.mx1.ROS1i <- pm[vecx1 >= jnmx1.1,]  # ROS covers first JN bound through maximum 
					mx1.pad <- matrix(nrow = nrow(pm) - nrow(pm.mx1.ROS1i), ncol = ncol(pm))
					pm.mx1.wROS <- rbind(mx1.pad, pm.mx1.ROS1i)
				}
				# upper JN bound within range (low band)
				if(!is.na(jnmx1.2) & is.na(jnmx1.1)) {
					pm.mx1.ROS2i <- pm[vecx1 <= jnmx1.2,] # ROS covers minimum through second JN bound
					mx1.pad <- matrix(nrow = nrow(pm) - nrow(pm.mx1.ROS2i), ncol = ncol(pm))
					pm.mx1.wROS <- rbind(pm.mx1.ROS2i, mx1.pad)
				}
			}
			
			# ROS outside JN bounds
			if(jnmx1.inside() == FALSE) {
				# first JN bound within range (low band)
				if(!is.na(jnmx1.1)) { 
					pm.mx1.ROS1o <- pm[vecx1 <= jnmx1.1,] # ROS covers minimum through first JN bound
				}
				# second JN bound within range (high band)
				if(!is.na(jnmx1.2)) { 
					pm.mx1.ROS2o <- pm[vecx1 >= jnmx1.2,] # ROS covers JN bound through maximum
				}
				
				# combining with pads
				# both within range (low and high)
				if(!is.na(jnmx1.1) & !is.na(jnmx1.2)) { 
					mx1.pad <- matrix(nrow = nrow(pm) - nrow(pm.mx1.ROS1o) - nrow(pm.mx1.ROS2o), ncol = ncol(pm)) # pad made of NAs
					pm.mx1.wROS <- rbind(pm.mx1.ROS1o, mx1.pad, pm.mx1.ROS2o) # pad between ROS with extra rows
				} 
				# first within range only (low band only)
				if (!is.na(jnmx1.1) & is.na(jnmx1.2)) { 
					mx1.pad <- matrix(nrow = nrow(pm) - nrow(pm.mx1.ROS1o), ncol = ncol(pm))
					pm.mx1.wROS <- rbind(pm.mx1.ROS1o, mx1.pad)
				}
				# second within range only (high band only)
				if (is.na(jnmx1.1) & !is.na(jnmx1.2)) { 
					mx1.pad <- matrix(nrow = nrow(pm) - nrow(pm.mx1.ROS2o), ncol = ncol(pm))
					pm.mx1.wROS <- rbind(mx1.pad, pm.mx1.ROS2o)
				}
			}
			# add vector values to row and column names of matrix
			rownames(pm.mx1.wROS) <- vecx1
			colnames(pm.mx1.wROS) <- vecx2
			pm.mx1.wROS
		} else {
			NULL
		}
	})
	
	
	# for x2 as mod
	# identify whether ROS is inside JN values
	jnmx2.inside <- eventReactive(input$run, {
		jnmx2()$inside
	})
	
	# identify whether ROS covers whole range of x2
	mx2.ROS.all <- eventReactive(input$run, {
		if((jnmx2.inside() == FALSE & (jnmx2()$bounds[2] < minx2() | jnmx2()$bounds[1] > maxx2())) |
				(jnmx2.inside() == TRUE & jnmx2()$bounds[1] < minx2() & jnmx2()$bounds[2] > maxx2())) {
			mx2.ROS.all <- TRUE
			mx2.ROS.all
		} else {
			mx2.ROS.all <- FALSE
			mx2.ROS.all
		}
	})
	
	# identify whether ROS within range of x2
	mx2.ROS.none <- eventReactive(input$run, {
		if((is.na(jnmx2.1()) & is.na(jnmx2.2())) & mx2.ROS.all() == FALSE) {
			mx2.ROS.none <- TRUE
			mx2.ROS.none
		} else {
			mx2.ROS.none <- FALSE
			mx2.ROS.none
		}
	})
	
	# identify column number of jn bounds
	jnmx2.1cn <- eventReactive(input$run, {
		jnmx2.1cn <- match(jnmx2.1(), vecx2())
	})
	
	jnmx2.2cn <- eventReactive(input$run, {
		jnmx2.2cn <- match(jnmx2.2(), vecx2())
	})
	
	
	# make x2 ROS planes
	pm.mx2.wROS <- eventReactive(input$run, {
		if (mx2.ROS.none() == FALSE & intsig() == TRUE) {
			
			pm <- pm()
			jnmx2.1 <- jnmx2.1()
			jnmx2.2 <- jnmx2.2()
			vecx1 <- vecx1()
			vecx2 <- vecx2()
			
			# ROS covers whole range
			if(mx2.ROS.all()==TRUE) {
				pm.mx2.wROS <- pm() 
			}
			
			
			# ROS inside JN bounds
			if(jnmx2.inside() == TRUE) {
				# inside both JN values (center band) 
				if(!is.na(jnmx2.1) & !is.na(jnmx2.2)) {
					pm.mx2.ROS1i <- pm[,jnmx2.1cn():jnmx2.2cn()]
					mx2.pad1 <- matrix(nrow = nrow(pm), ncol = jnmx2.1cn() - 1) # lower pad made of NAs
					mx2.pad2 <- matrix(nrow = nrow(pm), ncol = ncol(pm) - jnmx2.2cn()) # upper pad  made of NAs
					pm.mx2.wROS <- cbind(mx2.pad1, pm.mx2.ROS1i, mx2.pad2)
				}
				# lower JN bound within range (high band)
				if(!is.na(jnmx2.1) & is.na(jnmx2.2)) {
					pm.mx2.ROS1i <- pm[,vecx2 >= jnmx2.1]  # ROS covers first JN bound through maximum 
					mx2.pad <- matrix(nrow = nrow(pm), ncol = ncol(pm) - ncol(pm.mx2.ROS1i))
					pm.mx2.wROS <- cbind(mx2.pad, pm.mx2.ROS1i)
				}
				# upper JN bound within range (low band)
				if(!is.na(jnmx2.2) & is.na(jnmx2.1)) {
					pm.mx2.ROS2i <-	pm[,vecx2 <= jnmx2.2] # ROS covers minimum through second JN bound
					mx2.pad <- matrix(nrow = nrow(pm), ncol = ncol(pm) - ncol(pm.mx2.ROS2i))
					pm.mx2.wROS <- cbind(pm.mx2.ROS2i, mx2.pad)
				}
			}
			
			# ROS outside JN bounds
			if(jnmx2.inside() == FALSE) {
				
				# first JN bound within range (low band)
				if(!is.na(jnmx2.1)) { 
					pm.mx2.ROS1o <- pm[,vecx2 <= jnmx2.1] # ROS covers minimum through first JN bound
				}
				# second JN bound within range (high band)
				if(!is.na(jnmx2.2)) { 
					pm.mx2.ROS2o <- pm[,vecx2 >= jnmx2.2] # ROS covers JN bound through maximum
				}
				
				# combining with pads
				# both within range (low and high)
				if(!is.na(jnmx2.1) & !is.na(jnmx2.2)) { 
					mx2.pad <- matrix(nrow = nrow(pm), ncol = ncol(pm) - ncol(pm.mx2.ROS1o) - ncol(pm.mx2.ROS2o))
					pm.mx2.wROS <- cbind(pm.mx2.ROS1o, mx2.pad, pm.mx2.ROS2o)
				} 
				# first within range only (low band only)
				if (!is.na(jnmx2.1) & is.na(jnmx2.2)) { 
					mx2.pad <- matrix(nrow = nrow(pm), ncol = ncol(pm) - ncol(pm.mx2.ROS1o))
					pm.mx2.wROS <- cbind(pm.mx2.ROS1o, mx2.pad)
				}
				# second within range only (high band only)
				if (is.na(jnmx2.1) & !is.na(jnmx2.2)) { 
					mx2.pad <- matrix(nrow = nrow(pm), ncol = ncol(pm) - ncol(pm.mx2.ROS2o))
					pm.mx2.wROS <- cbind(mx2.pad, pm.mx2.ROS2o)
				}
			}
			
			# add vector values to row and column names of matrix
			rownames(pm.mx2.wROS) <- vecx1
			colnames(pm.mx2.wROS) <- vecx2
			pm.mx2.wROS
		} else {
			NULL
		}
	})
	
	
	# set reactive values for plot features for continuity  --------
	checkboxes <- reactiveValues(scatterCB = TRUE, regplaneCB = TRUE,
		predCICB = FALSE, pm.mx1.wROSCB = "none", pm.mx2.wROSCB = "none",
		jn.gradient.mx1CB = FALSE, jn.gradient.mx2CB = FALSE,
		coonx1CB = FALSE, coonx2 = FALSE, fdrRawDataCB = TRUE, fdrManualCB = TRUE)
	
	# fdr
	observeEvent(input$fdrRawDataCheckbox, {
		if(input$fdrRawDataCheckbox == TRUE) {
			checkboxes$fdrRawDataCB <- TRUE
		}
		if(input$fdrRawDataCheckbox == FALSE) {
			checkboxes$fdrRawDataCB <- FALSE
		}
	})
	
	observeEvent(input$fdrManualCheckbox, {
		if(input$fdrManualCheckbox == TRUE) {
			checkboxes$fdrManualCB <- TRUE
		}
		if(input$fdrManualCheckbox == FALSE) {
			checkboxes$fdrManualCB <- FALSE
		}
	})
	
	# scatter
	observeEvent(input$scatterCheckbox, {
		if(input$scatterCheckbox == TRUE) {
			checkboxes$scatterCB <- TRUE
		}
		if(input$scatterCheckbox == FALSE) {
			checkboxes$scatterCB <- FALSE
		}
	})
	
	# reg plane
	observeEvent(input$regplaneCheckbox, {
		if(input$regplaneCheckbox == TRUE) {
			checkboxes$regplaneCB <- TRUE
		}
		if(input$regplaneCheckbox == FALSE) {
			checkboxes$regplaneCB <- FALSE
		}
	})
	
	# pred CI
	observeEvent(input$predCICheckbox, {
		if(input$predCICheckbox == TRUE) {
			checkboxes$predCICB <- TRUE
		}
		if(input$predCICheckbox == FALSE) {
			checkboxes$predCICB <- FALSE
		}
	})
	
	# ROS mod x1
	observeEvent(input$pm.mx1.wROSCheckbox, {
		if(input$pm.mx1.wROSCheckbox == "gradient") {
			checkboxes$pm.mx1.wROSCB <- "gradient"
		}
		if(input$pm.mx1.wROSCheckbox == "gradientAll") {
			checkboxes$pm.mx1.wROSCB <- "gradientAll"
		}
		if(input$pm.mx1.wROSCheckbox == "none") {
			checkboxes$pm.mx1.wROSCB <- "none"
		}
		if(input$pm.mx1.wROSCheckbox == "solid") {
			checkboxes$pm.mx1.wROSCB <- "solid"
		}
	})
	
	# ROS mod x2
	observeEvent(input$pm.mx2.wROSCheckbox, {
		if(input$pm.mx2.wROSCheckbox == "gradient") {
			checkboxes$pm.mx2.wROSCB <- "gradient"
		}
		if(input$pm.mx2.wROSCheckbox == "gradientAll") {
			checkboxes$pm.mx2.wROSCB <- "gradientAll"
		}
		if(input$pm.mx2.wROSCheckbox == "none") {
			checkboxes$pm.mx2.wROSCB <- "none"
		}
		if(input$pm.mx2.wROSCheckbox == "solid") {
			checkboxes$pm.mx2.wROSCB <- "solid"
		}
	})
	
	# crossover on x1
	observeEvent(input$coonx1Checkbox, {
		if(input$coonx1Checkbox == TRUE) {
			checkboxes$coonx1CB <- TRUE
		}
		if(input$coonx1Checkbox == FALSE) {
			checkboxes$coonx1CB <- FALSE
		}
	})
	
	# crossover on x2
	observeEvent(input$coonx2Checkbox, {
		if(input$coonx2Checkbox == TRUE) {
			checkboxes$coonx2CB <- TRUE
		}
		if(input$coonx2Checkbox == FALSE) {
			checkboxes$coonx2CB <- FALSE
		}
	})
	
	x1.name.plot <- eventReactive(input$run, {
		x1.name.plot <- input$x1
	})
	
	x2.name.plot <- eventReactive(input$run, {
		x2.name.plot <- input$x2
	})
	
	y.name.plot <- eventReactive(input$run, {
		y.name.plot <- input$y
	})
	
	
	# create 3D plot ---------------------------------------------------------
	plot3d <- eventReactive(c(input$run, input$regplaneCheckbox, input$scatterCheckbox, input$predCICheckbox, input$pm.mx1.wROSCheckbox, input$pm.mx2.wROSCheckbox, input$coonx1Checkbox, input$coonx2Checkbox), {
		pmcolor <- array(1L, dim(pm())) # matrix for regression plane color
		if (input$inputType == "rawData") {
			# req(input$x1 != "--", input$x2 != "--",input$y != "--", sameVar() == FALSE, non.numeric() == FALSE)
			data <- df.std()
			
			# scatter plot toggle
			scattersize <- 1.5
			scatterop <- 0
			if (!is.null(input$scatterCheckbox)) {
				if (input$scatterCheckbox == TRUE) {
					scattersize <- 1.5
					scatterop <- 1
				} else {
					scattersize <- 0.01
					scatterop <- 0
				}
			} else {
				NULL
			}
			
			p <- plot_ly(width = 575, height = 450) %>% 
				add_trace(data = data,
					x = ~data[,x1.name.plot()],
					y = ~data[,x2.name.plot()],
					z = ~data[,y.name.plot()],
					type = "scatter3d",
					mode = "markers",
					marker = list(size = scattersize,
						color = "grey",
						opacity = scatterop)) %>%
				# plot formatting
				layout(
					scene = list(
						camera = list(eye = list(x = 1.6, y = -1.6, z = 1.25)),
						xaxis = list(title = x1.name(),
							titlefont = list(family = "Arial, sans-serif", color = x2modcol(), size = 10),
							tickfont = list(family = "Arial, sans-serif", size = 10)),
						yaxis = list(title = x2.name(),
							titlefont = list(family = "Arial, sans-serif", color = x1modcol(), size = 10),
							tickfont = list(family = "Arial, sans-serif", size = 10)),
						zaxis = list(title = y.name(),
							titlefont = list(family = "Arial, sans-serif", color = "#484747", size = 10),
							tickfont = list(family = "Arial, sans-serif",	size = 10)))) %>%
				hide_legend()
			
		} else { 
			req(anyNA() == FALSE, notRight() == FALSE, input$dfs > 1, negVar() == FALSE, minGEMaxX1() == FALSE, minGEMaxX2() == FALSE)
			p <- plot_ly(width = 575, height = 450) %>% 
				add_trace(x = c(minx1(), maxx1()),
					y = c(minx2(), maxx2()),
					z = c(min(pm()), max(pm())),
					type = "scatter3d",
					mode = "markers",
					marker = list(size = .01,
						color = "grey",
						opacity = 0)) %>%
				# plot formatting
				layout(
					scene = list(
						camera = list(eye = list(x = 1.6, y = -1.6, z = 1.25)),
						xaxis = list(title = x1.name(),
							titlefont = list(family = "Arial, sans-serif", color = x2modcol(), size = 10),
							tickfont = list(family = "Arial, sans-serif", size = 10)),
						yaxis = list(title = x2.name(),
							titlefont = list(family = "Arial, sans-serif", color = x1modcol(), size = 10),
							tickfont = list(family = "Arial, sans-serif", size = 10)),
						zaxis = list(title = y.name(),
							titlefont = list(family = "Arial, sans-serif", color = "#484747", size = 10),
							tickfont = list(family = "Arial, sans-serif",	size = 10)))) %>%
				hide_legend()
			
		} 
		
		# add regression plane if checked
		if (!is.null(input$regplaneCheckbox)) {
			if (input$regplaneCheckbox == TRUE) {
				p <- add_surface(p, x = ~vecx1(), y = ~vecx2(), z = ~t(pm()),
					type = "surface",
					# text = customhover,
					# hoverinfo = "text",
					surfacecolor = t(pmcolor),
					colorscale = list(c(0, 1), c('gray', 'gray')),
					opacity = 0.55,
					showscale = FALSE)
			}
		} else {
			NULL
		}
		
		# confidence interval for predicted values
		if(!is.null(input$predCICheckbox)) {
			if (input$predCICheckbox == TRUE) {
				# lower confidence band
				p <- add_surface(p, x = ~vecx1(), y = ~vecx2(), z = ~t(pmcil()),
					type = "surface",
					surfacecolor=t(pmcolor),
					colorscale = list(c(0, 1), c('gray', 'gray')),
					opacity = 0.4,
					showscale = FALSE) %>%
					# upper confidence band
					add_surface(x = ~vecx1(), y = ~vecx2(), z = ~t(pmciu()),
						type = "surface",
						surfacecolor= t(pmcolor),
						colorscale = list(c(0, 1), c('gray', 'gray')),
						opacity = 0.4,
						showscale = FALSE)
			}
		} else {
			NULL
		}
		
		# ROS for x2 as mod
		if (mx2.ROS.none() == FALSE & intsig() == TRUE & !is.null(input$pm.mx2.wROSCheckbox)) {
			if (input$pm.mx2.wROSCheckbox == "gradient") {
				pm.mx2.wROSg <- t(pm.mx2.wROS())
				pm.mx2.wROSg[!is.na(pm.mx2.wROSg)] <- jnmx2cbd()
				pm.mx2.wROSg <- t(pm.mx2.wROSg)
				p <- add_surface(p, x = ~vecx1(), y = ~vecx2(), z = ~t(pm.mx2.wROS()),
					type = "surface",
					surfacecolor=t(pm.mx2.wROSg),
					colorscale = list(c(0, 1), c(x2modcol(), 'white')),
					opacity = 0.7,
					showscale = FALSE)
			}
			if(input$pm.mx2.wROSCheckbox == "solid") {
				pm.mx2.wROSg <- pm.mx2.wROS()
				pm.mx2.wROSg[!is.na(pm.mx2.wROSg)] <- 1
				p <- add_surface(p, x = ~vecx1(), y = ~vecx2(), z = ~t(pm.mx2.wROS()),
					type = "surface",
					surfacecolor=t(pm.mx2.wROSg),
					colorscale = list(c(0, 1), c(x2modcol(), 'white')),
					opacity = 0.7,
					showscale = FALSE)
			}
			if (input$pm.mx2.wROSCheckbox == "gradientAll") {
				pm.mx2.wROS <- pm()
				pm.mx2.wROSg <- t(pm())
				# pm.mx2.wROSg <- t(pm.mx2.wROS())
				pm.mx2.wROSg[!is.na(pm.mx2.wROSg)] <- jnmx2cbd()
				# pm.mx2.wROSg[!is.na(pm.mx2.wROSg)] <- jnmx2cbd()
				pm.mx2.wROSg <- t(pm.mx2.wROSg)
				p <- add_surface(p, x = ~vecx1(), y = ~vecx2(), z = ~t(pm.mx2.wROS),
					type = "surface",
					surfacecolor=t(pm.mx2.wROSg),
					colorscale = list(c(0, 1), c(x2modcol(), 'white')),
					opacity = 0.7,
					showscale = FALSE)
			}
		} else {
			NULL
		}
		
		# ROS for x1 as mod
		if (mx1.ROS.none() == FALSE & intsig() == TRUE & !is.null(input$pm.mx1.wROSCheckbox)) {
			if(input$pm.mx1.wROSCheckbox == "gradient") {
				pm.mx1.wROSg <- pm.mx1.wROS()
				pm.mx1.wROSg[!is.na(pm.mx1.wROSg)] <- jnmx1cbd()
				p <- add_surface(p, x = ~vecx1(), y = ~vecx2(), z = ~t(pm.mx1.wROS()),
					type = "surface",
					surfacecolor=t(pm.mx1.wROSg),
					colorscale = list(c(0, 1), c(x1modcol(), 'white')),
					opacity = 0.7,
					showscale = FALSE)
			}
			if(input$pm.mx1.wROSCheckbox == "solid") {
				pm.mx1.wROSg <- pm.mx1.wROS()
				pm.mx1.wROSg[!is.na(pm.mx1.wROSg)] <- 1
				p <- add_surface(p, x = ~vecx1(), y = ~vecx2(), z = ~t(pm.mx1.wROS()),
					type = "surface",
					surfacecolor=t(pm.mx1.wROSg),
					colorscale = list(c(0, 1), c(x1modcol(), 'white')),
					opacity = 0.7,
					showscale = FALSE)
			}
			if (input$pm.mx1.wROSCheckbox == "gradientAll") {
				pm.mx1.wROS <- pm()
				pm.mx1.wROSg <- pm()
				pm.mx1.wROSg[!is.na(pm.mx1.wROSg)] <- jnmx1cbd()
				p <- add_surface(p, x = ~vecx1(), y = ~vecx2(), z = ~t(pm.mx1.wROS),
					type = "surface",
					surfacecolor=t(pm.mx1.wROSg),
					colorscale = list(c(0, 1), c(x1modcol(), 'white')),
					opacity = 0.7,
					showscale = FALSE)
			}
		} else {
			NULL
		}
		
		# crossover line for x2 as mod
		if (!is.null(input$coonx1Checkbox)) {
			if(input$coonx1Checkbox == TRUE) {
				p <- add_trace(p, x = coonx1(), y = vecx2(), z = ~coonx1y(),
					type = "scatter3d",
					mode = "lines",
					line = list(width = 2, dash = "solid", color = x1modcol()),
					opacity = 0.8)
			}
		} else {
			NULL
		}
		
		# crossover line for x1 as mod
		if (!is.null(input$coonx2Checkbox)) {
			if (input$coonx2Checkbox == TRUE) {
				p <- add_trace(p, x = vecx1(), y = coonx2(), z = ~coonx2y(),
					type = "scatter3d",
					mode = "lines",
					line = list(width = 2, dash = "solid", color = x2modcol()),
					opacity = 0.7)
			}
		} else {
			NULL
		}
		p
	})
	
	
	# scale for cb gradient ---------------------------------------------------
	# for x1 as mod	
	jnmx1cbd <- eventReactive(c(input$run, input$pm.mx1.wROSCheckbox), {
		if(intsig() == TRUE & !is.null(input$pm.mx1.wROSCheckbox) & mx1.ROS.none() == FALSE) {
			if(input$pm.mx1.wROSCheckbox == "gradient") {
				jnmx1cb <- jnmx1cb()
				if(mx1.ROS.all() == TRUE) {
					jnmx1cb <- jnmx1cb[!is.na(jnmx1cb[,1]),]
				} else {
					if(jnmx1.inside() == FALSE) {
						jnmx1cb <- jnmx1cb[vecx1() <= jnmx1.1() | vecx1() >= jnmx1.2(),]
					}
					if(jnmx1.inside() == TRUE) {
						jnmx1cb <- jnmx1cb[vecx1() >= jnmx1.1() | vecx1() <= jnmx1.2(),]
					}
					jnmx1cb <- jnmx1cb[!is.na(jnmx1cb[,1]),]
				}
				jnmx1cbd <- (jnmx1cb[,4] - jnmx1cb[,3]) # upper CI minus lower CI
			} else {
				if(input$pm.mx1.wROSCheckbox == "gradientAll") {
					jnmx1cb <- jnmx1cb()
					jnmx1cb <- jnmx1cb[!is.na(jnmx1cb[,1]),]
					jnmx1cbd <- (jnmx1cb[,4] - jnmx1cb[,3]) # upper CI minus lower CI
				}
			} 
		} else {
			NULL
		}
	})
	
	jnmx1cbd.print <- eventReactive(c(input$run, input$pm.mx1.wROSCheckbox), {
		jnmx1cbd <- seq(jnmx1cbd()[1], jnmx1cbd()[length(jnmx1cbd())], length.out = 10)
		jnmx1cbd <- sort(jnmx1cbd, decreasing = FALSE)
		jnmx1cbd <- format(round(jnmx1cbd, digits = 2), nsmall = 2)
		jnmx1cbd
	})
	
	mx1cbscale <- eventReactive(c(input$run, input$pm.mx1.wROSCheckbox), {
		if (input$inputType == "rawData") {
			req(non.numeric() == FALSE, sameVar() == FALSE)
		} else {
			req(anyNA() == FALSE, notRight() == FALSE, input$dfs > 1, negVar() == FALSE, minGEMaxX1() == FALSE, minGEMaxX2() == FALSE) 
		}
		if(intsig() == TRUE & !is.null(input$pm.mx1.wROSCheckbox) & mx1.ROS.none() == FALSE) {
			if(input$pm.mx1.wROSCheckbox == "gradient" | input$pm.mx1.wROSCheckbox == "gradientAll") {
				HTML("<div style='background-image: linear-gradient(to right, ", x1modcol(), ", white);  
							border-width: 1px; border-style: solid; border-color:gray;
							height: 20px;   text-align: justify' </div>")
			}
		} else {
			NULL
		}
	})
	
	output$mx1cbscale <- renderUI({
		mx1cbscale()
	})
	
	mx1cbscale.title <- eventReactive(c(input$run, input$pm.mx1.wROSCheckbox), {
		if (input$inputType == "rawData") {
			req(non.numeric() == FALSE, sameVar() == FALSE)
		} else {
			req(anyNA() == FALSE, notRight() == FALSE, input$dfs > 1, negVar() == FALSE, minGEMaxX1() == FALSE, minGEMaxX2() == FALSE) 
		}
		if(intsig() == TRUE & !is.null(input$pm.mx1.wROSCheckbox) & mx1.ROS.none() == FALSE) {
			if(input$pm.mx1.wROSCheckbox == "gradient" | input$pm.mx1.wROSCheckbox == "gradientAll") {
				HTML("<div style='height: 20px;   text-align: center' >", paste0("Width of 95% CB around ", x2.name(), " Slope"), "</div>")
			}
		} else {
			NULL
		}
	})
	
	output$mx1cbscale.title <- renderUI({
		mx1cbscale.title()
	})
	
	mx1cbscale.text1 <- eventReactive(c(input$run, input$pm.mx1.wROSCheckbox), {
		if (input$inputType == "rawData") {
			req(non.numeric() == FALSE, sameVar() == FALSE)
		} else {
			req(anyNA() == FALSE, notRight() == FALSE, input$dfs > 1, negVar() == FALSE, minGEMaxX1() == FALSE, minGEMaxX2() == FALSE) 
		}
		if(intsig() == TRUE & !is.null(input$pm.mx1.wROSCheckbox) & mx1.ROS.none() == FALSE) {
			if(input$pm.mx1.wROSCheckbox == "gradient" | input$pm.mx1.wROSCheckbox == "gradientAll") {
				HTML("<div style='height: 20px;   text-align: left' >",  paste0("Min: ", min(jnmx1cbd.print())), "</div>")
			}
		} else {
			NULL
		}
	})
	
	output$mx1cbscale.text1 <- renderUI({
		mx1cbscale.text1()
	})
	
	mx1cbscale.text2 <- eventReactive(c(input$run, input$pm.mx1.wROSCheckbox), {
		if (input$inputType == "rawData") {
			req(non.numeric() == FALSE, sameVar() == FALSE)
		} else {
			req(anyNA() == FALSE, notRight() == FALSE, input$dfs > 1, negVar() == FALSE, minGEMaxX1() == FALSE, minGEMaxX2() == FALSE) 
		}
		if(intsig() == TRUE & !is.null(input$pm.mx1.wROSCheckbox) & mx1.ROS.none() == FALSE) {
			if(input$pm.mx1.wROSCheckbox == "gradient" | input$pm.mx1.wROSCheckbox == "gradientAll") {
				HTML("<div style='height: 20px;   text-align: right' >",  paste0("Max: ", max(jnmx1cbd.print())), "</div>") 
			}
		} else {
			NULL
		}
	})
	
	output$mx1cbscale.text2 <- renderUI({
		mx1cbscale.text2()
	})
	
	# for x2 as mod
	jnmx2cbd <- eventReactive(c(input$run, input$pm.mx2.wROSCheckbox), {
		if(intsig() == TRUE & !is.null(input$pm.mx2.wROSCheckbox) & mx2.ROS.none() == FALSE) {
			if(input$pm.mx2.wROSCheckbox == "gradient") {
				jnmx2cb <- jnmx2cb()
				if(mx2.ROS.all() == TRUE) {
					jnmx2cb <- jnmx2cb[!is.na(jnmx2cb[,1]),]
				} else {
					if(jnmx2.inside() == FALSE) {
						jnmx2cb <- jnmx2cb[vecx2() <= jnmx2.1() | vecx2() >= jnmx2.2(),]
					} 
					if (jnmx2.inside() == TRUE) {
						jnmx2cb <- jnmx2cb[vecx2() >= jnmx2.1() | vecx2() <= jnmx2.2(),]
					}
					jnmx2cb <- jnmx2cb[!is.na(jnmx2cb[,1]),]
				}
				jnmx2cbd <- (jnmx2cb[,4] - jnmx2cb[,3]) # upper CI minus lower CI
			} else {
				if(input$pm.mx2.wROSCheckbox == "gradientAll") {
					jnmx2cb <- jnmx2cb()
					jnmx2cb <- jnmx2cb[!is.na(jnmx2cb[,1]),]
					jnmx2cbd <- (jnmx2cb[,4] - jnmx2cb[,3]) # upper CI minus lower CI
				}
			} 
		} else {
			NULL
		}
	})
	
	jnmx2cbd.print <- eventReactive(c(input$run, input$pm.mx2.wROSCheckbox), {	
		jnmx2cbd <- seq(jnmx2cbd()[1], jnmx2cbd()[length(jnmx2cbd())], length.out = 10)
		jnmx2cbd <- sort(jnmx2cbd, decreasing = FALSE)
		jnmx2cbd <- format(round(jnmx2cbd, digits = 2), nsmall = 2)
		jnmx2cbd
	})
	
	mx2cbscale <- eventReactive(c(input$run, input$pm.mx2.wROSCheckbox), {
		if (input$inputType == "rawData") {
			req(non.numeric() == FALSE, sameVar() == FALSE)
		} else {
			req(anyNA() == FALSE, notRight() == FALSE, input$dfs > 1, negVar() == FALSE, minGEMaxX1() == FALSE, minGEMaxX2() == FALSE) 
		}
		if(intsig() == TRUE & !is.null(input$pm.mx2.wROSCheckbox) & mx2.ROS.none() == FALSE) {
			if(input$pm.mx2.wROSCheckbox == "gradient" | input$pm.mx2.wROSCheckbox == "gradientAll") {
				HTML("<div style='background-image: linear-gradient(to right, ", x2modcol(), ", white);  
							border-width: 1px; border-style: solid; border-color:gray;
							height: 20px;   text-align: justify' </div>")
			}
		} else {
			NULL
		}
	})
	
	output$mx2cbscale <- renderUI({
		mx2cbscale()
	})
	
	mx2cbscale.title <- eventReactive(c(input$run, input$pm.mx2.wROSCheckbox), {
		if (input$inputType == "rawData") {
			req(non.numeric() == FALSE, sameVar() == FALSE)
		} else {
			req(anyNA() == FALSE, notRight() == FALSE, input$dfs > 1, negVar() == FALSE, minGEMaxX1() == FALSE, minGEMaxX2() == FALSE) 
		}
		if(intsig() == TRUE & !is.null(input$pm.mx2.wROSCheckbox) & mx2.ROS.none() == FALSE) {
			if(input$pm.mx2.wROSCheckbox == "gradient" | input$pm.mx2.wROSCheckbox == "gradientAll") {
				HTML("<div style='height: 20px;   text-align: center' >", paste0("Width of 95% CB around ", x1.name(), " Slope"), "</div>")
			}
		} else {
			NULL
		}
	})
	
	output$mx2cbscale.title <- renderUI({
		mx2cbscale.title()
	})
	
	mx2cbscale.text1 <- eventReactive(c(input$run, input$pm.mx2.wROSCheckbox), {
		if (input$inputType == "rawData") {
			req(non.numeric() == FALSE, sameVar() == FALSE)
		} else {
			req(anyNA() == FALSE, notRight() == FALSE, input$dfs > 1, negVar() == FALSE, minGEMaxX1() == FALSE, minGEMaxX2() == FALSE) 
		}
		if(intsig() == TRUE & !is.null(input$pm.mx2.wROSCheckbox) & mx2.ROS.none() == FALSE) {
			if(input$pm.mx2.wROSCheckbox == "gradient" | input$pm.mx2.wROSCheckbox == "gradientAll") {
				HTML("<div style='height: 20px;   text-align: left' >",  paste0("Min: ", min(jnmx2cbd.print())), "</div>")
			}
		} else {
			NULL
		}
	})
	
	output$mx2cbscale.text1 <- renderUI({
		mx2cbscale.text1()
	})
	
	mx2cbscale.text2 <- eventReactive(c(input$run, input$pm.mx2.wROSCheckbox), {
		if (input$inputType == "rawData") {
			req(non.numeric() == FALSE, sameVar() == FALSE)
		} else {
			req(anyNA() == FALSE, notRight() == FALSE, input$dfs > 1, negVar() == FALSE, minGEMaxX1() == FALSE, minGEMaxX2() == FALSE) 
		}
		if(intsig() == TRUE & !is.null(input$pm.mx2.wROSCheckbox) & mx2.ROS.none() == FALSE) {
			if(input$pm.mx2.wROSCheckbox == "gradient" | input$pm.mx2.wROSCheckbox == "gradientAll") {
				HTML("<div style='height: 20px;   text-align: right' >",  paste0("Max: ", max(jnmx2cbd.print())), "</div>") 
			}
		} else {
			NULL
		}
	})
	
	output$mx2cbscale.text2 <- renderUI({
		mx2cbscale.text2()
	})
	
	# render plot -------------------------------------------------------------
	
	output$plot3dPlot <- renderPlotly({
		# if (input$inputType == "rawData") {
		# 	req(non.numeric() == FALSE, sameVar() == FALSE)
		# } else {
		# 	req(anyNA() == FALSE, input$dfs > 1, negVar() == FALSE, minGEMaxX1() == FALSE, minGEMaxX2() == FALSE)
		# }
		plot3d()
	})
	
	# render UI for plot features ---------------------------------------------
	
	scatterUI <- eventReactive(input$run, {
		if (input$inputType == "rawData") {
			req(input$x1 != '--', input$x2 != '--', input$y != '--', sameVar() == FALSE, non.numeric() == FALSE)
			checkboxInput(inputId = "scatterCheckbox",
				label = "Scatterplot",
				value = checkboxes$scatterCB)
		}
	})
	output$scatterUI <- renderUI({ 
		scatterUI()
	})
	regplaneUI <- eventReactive(input$run, {
		if (input$inputType == "rawData") {
			req(input$x1 != '--', input$x2 != '--', input$y != '--', sameVar() == FALSE, non.numeric() == FALSE)
		} else {
			req(anyNA() == FALSE, int0() == FALSE, notRight() == FALSE, input$dfs > 1, negVar() == FALSE, minGEMaxX1() == FALSE, minGEMaxX2() == FALSE)
		}
		checkboxInput(inputId = "regplaneCheckbox",
			label = "Regression plane",
			value = checkboxes$regplaneCB)
	}) 
	output$regplaneUI <- renderUI({ 
		regplaneUI()
	})
	predCIUI <- eventReactive(input$run, {
		req(input$inputType == "rawData", input$x1 != '--', input$x2 != '--', input$y != '--', sameVar() == FALSE, non.numeric() == FALSE)
		checkboxInput(inputId = "predCICheckbox", 
			label = "95% CI",
			value = checkboxes$predCICB)
	})
	output$predCIUI <- renderUI({ 
		predCIUI()
	})
	output$pm.mx1.wROSUI <- renderUI({
		if (input$inputType == "rawData") {
			req(mx1.ROS.none() == FALSE, intsig() == TRUE, non.numeric() == FALSE)
		} else {
			req(mx1.ROS.none() == FALSE, intsig() == TRUE)
		}
		radioButtons(inputId = "pm.mx1.wROSCheckbox",
			label = NULL,
			choices = c("No ROS" = "none", "Solid ROS" = "solid", "Gradient ROS (Slope 95% CB)" = "gradient", 
				"Gradient All (Slope 95% CB)" = "gradientAll"),
			selected = checkboxes$pm.mx1.wROSCB, inline = FALSE)
	})
	output$pm.mx2.wROSUI <- renderUI({
		if (input$inputType == "rawData") {
			req(mx2.ROS.none() == FALSE, intsig() == TRUE, non.numeric() == FALSE)
		} else {
			req(mx2.ROS.none() == FALSE, intsig() == TRUE)
		}
		radioButtons(inputId = "pm.mx2.wROSCheckbox",
			label = NULL,
			choices = c("No ROS" = "none", "Solid ROS" = "solid", "Gradient ROS (Slope 95% CB)" = "gradient",
				"Gradient All (Slope 95% CB)" = "gradientAll"),
			selected = checkboxes$pm.mx2.wROSCB, inline = FALSE)
	})
	
	coonx1UI <- eventReactive(input$run, {
		if (input$inputType == "rawData") {
			req(!is.na(coonx1()), non.numeric() == FALSE, sameVar() == FALSE)
		} else {
			req(!is.na(coonx1()), anyNA() == FALSE, int0() == FALSE, notRight() == FALSE, input$dfs > 1, negVar() == FALSE, minGEMaxX1() == FALSE, minGEMaxX2() == FALSE)
		}
		checkboxInput(inputId = "coonx1Checkbox", 
			label = "Crossover",
			value = checkboxes$coonx1CB)
	})
	output$coonx1UI <- renderUI({
		coonx1UI()
	})
	
	coonx2UI <- eventReactive(input$run, {
		if (input$inputType == "rawData") {
			req(!is.na(coonx2()), non.numeric() == FALSE, sameVar() == FALSE)
		} else {
			req(!is.na(coonx2()), anyNA() == FALSE, int0() == FALSE, notRight() == FALSE, input$dfs > 1, negVar() == FALSE, minGEMaxX1() == FALSE, minGEMaxX2() == FALSE)
		}
		coonx2 <- format(round(coonx2(), digits = 2), nsmall = 2)
		checkboxInput(inputId = "coonx2Checkbox", 
			label = "Crossover",
			value = checkboxes$coonx2CB)
	})
	
	output$coonx2UI <- renderUI({
		coonx2UI()
	})
		
	mx1CheckboxTitle <- eventReactive(input$run, {
		if (input$inputType == "rawData") {
			req(input$x1 != '--', input$x2 != '--', input$y != '--', sameVar() == FALSE, non.numeric() == FALSE)
		} else {
			req(anyNA() == FALSE, int0() == FALSE, notRight() == FALSE, input$dfs > 1, negVar() == FALSE, minGEMaxX1() == FALSE, minGEMaxX2() == FALSE)
		}
		HTML(paste0("<div style='font-weight: 700' >", x1.name(), " as Moderator </div>")) 
	})
	
	output$mx1CheckboxTitle <- renderText({
		mx1CheckboxTitle()
	})
	
	mx2CheckboxTitle <- eventReactive(input$run, {
		if (input$inputType == "rawData") {
			req(input$x1 != '--', input$x2 != '--', input$y != '--', sameVar() == FALSE, non.numeric() == FALSE)
		} else {
			req(anyNA() == FALSE, int0() == FALSE, notRight() == FALSE, input$dfs > 1, negVar() == FALSE, minGEMaxX1() == FALSE, minGEMaxX2() == FALSE)
		}
		HTML(paste0("<div style='font-weight: 700' >", x2.name(), " as Moderator </div>")) 
	})
	
	output$mx2CheckboxTitle <- renderText({
		mx2CheckboxTitle()
	})
	
	
	# code for making live figures --------------------------------------------
	output$codetitle1 <- renderUI({
		HTML(paste(
			"<div style='font-size: 13px'> 
			<h4> R Code for Making 3D Figures </h4>
					To create a “live” rotatable plot with a permanent link like <a href=http://rpubs.com/sbu_mfinsaas/Example target=_blank> this one</a>:
						<ol> 
							<li> Download <a href=https://www.r-project.org target=_blank>R</a> and <a href=https://www.rstudio.com target=_blank>RStudio</a>, which are both free programs. </li>
							<li> Click <a href=https://github.com/mfinsaas/jnthreedimint/blob/master/1-FunctionCode_LivePlot.R target=_blank>here</a> to access the function code for producing a live plot and run it in RStudio. Any time you start a new R session you will need to re-run this code before running the input code (or you can save and source this file). </li>
							<li> Click <a href=https://github.com/mfinsaas/jnthreedimint/blob/master/2-InputCode_LivePlot.R target=_blank>here</a>
			to access the input code. Then update the arguments and run it in RStudio. The plot will appear in the RStudio Viewer and a copy will also be saved locally to your computer. </li>
							<li> Click the blue “Publish” icon in the RStudio viewer and follow the steps to publish to the free service <a href=http://rpubs.com target=_blank>Rpubs</a>. The permanent Rpubs link can then be emailed to colleagues, included in papers, etc.</li>
						</ol>
							This approach works because RStudio supports online publication directly from its viewer window, whereas this action is not supported directly from Shiny apps. If you have trouble with it, please email me at megan.finsaas@gmail.com."))
	})
	
	# educational materials --------------------------------------------
	output$edu <- renderUI({
		HTML(paste(
			"<div style='font-size: 13px'> 
			<h4> Educational Materials </h4>
					Check out these slides for additional help on <a href=https://drive.google.com/file/d/1bGE3fM4yWUjWiZ5HSgNFecLlJi9BIgT1/view?usp=sharing> visualizing interactions in 3D </a> and 
					<a href=https://drive.google.com/file/d/1NsgthZlswlV3i5B9OgXxnUpXmXEUsC-w/view?usp=sharing> conducting follow-up tests of continuous by continuous interactions. </a> 
						<ol>"))
			})
	
	
	# force select output to update on actionbutton ----------------------------------
	
	lapply(c("plot3dPlot", "mx1CheckboxTitle", "jnmx1.sum", "mx2CheckboxTitle", "jnmx2.sum",
		"mx1cbscale.title", "mx1cbscale", "mx1cbscale.text1", "mx1cbscale.text2", "mx2cbscale.title","mx2cbscale",
		"mx2cbscale.text1", "mx2cbscale.text2", "jnmx1Info",  "jnmx2Info", # can't include JN plots because they have unknown size before rendering 
		"mx1cbTitle", "jnmx1Text", "jnmx1cbInfo","mx2cbTitle", "jnmx2Text", "jnmx2cbInfo", "modFit", "desc.stdTitle","descTitle",
		"rawDataUnavail","desc.std"),
		function(x) {
			outputOptions(output, x, suspendWhenHidden = FALSE)	
		})
	
}

shinyApp(ui, server)

