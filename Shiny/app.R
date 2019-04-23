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
# returns whether ROS falues between J-N values
johnson_neyman_inside <- function (model, pred, modx, vmat = NULL, alpha = 0.05, plot = FALSE, 
	control.fdr = FALSE, line.thickness = 0.5, df = "residual", 
	digits = getOption("jtools-digits", 2), critical.t = NULL, 
	sig.color = "#00BFC4", insig.color = "#F8766D", mod.range = NULL, 
	title = "Johnson-Neyman plot") 
{
	predt <- as.character(substitute(pred))
	modxt <- as.character(substitute(modx))
	if (!(predt %in% as.character(attr(terms(model), "variables")))) {
		pred <- as.character(eval(pred))
		modx <- as.character(eval(modx))
	}
	else {
		pred <- predt
		modx <- modxt
	}
	if (df == "residual") {
		df <- df.residual(model)
	}
	else if (df == "normal") {
		df <- Inf
	}
	else if (!is.numeric(df)) {
		stop("df argument must be 'residual', 'normal', or a number.")
	}
	out <- list()
	out <- structure(out, pred = pred, modx = modx, alpha = alpha, 
		plot = plot, digits = digits, control.fdr = control.fdr)
	if (inherits(model, "merMod")) {
		intterm1 <- paste(pred, ":", modx, sep = "") 
		intterm1tf <- any(intterm1 %in% names(lme4::fixef(model))) 
		intterm2 <- paste(modx, ":", pred, sep = "") 
		intterm2tf <- any(intterm2 %in% names(lme4::fixef(model)))
		coefs <- lme4::fixef(model)
	}
	else {
		intterm1 <- paste(pred, ":", modx, sep = "")
		intterm1tf <- any(intterm1 %in% names(coef(model)))
		intterm2 <- paste(modx, ":", pred, sep = "")
		intterm2tf <- any(intterm2 %in% names(coef(model)))
		coefs <- coef(model) 
	}
	inttermstf <- c(intterm1tf, intterm2tf) 
	intterms <- c(intterm1, intterm2) 
	intterm <- intterms[which(inttermstf)]
	modrange <- range(model.frame(model)[, modx])
	modrangeo <- range(model.frame(model)[, modx])
	modsd <- sd(model.frame(model)[, modx])
	if (is.null(mod.range)) {
		modrange[1] <- modrange[1] - modsd
		modrange[2] <- modrange[2] + modsd
	}
	else {
		modrange <- mod.range
	}
	if (modrange[1] >= modrangeo[1] & modrange[2] <= modrangeo[2]) {
		no_range_line <- TRUE
	}
	else {
		no_range_line <- FALSE
	}
	alpha <- alpha/2
	if (control.fdr == FALSE) {
		tcrit <- qt(alpha, df = df)
		tcrit <- abs(tcrit)
	}
	else if (control.fdr == TRUE) {
		predb <- coefs[pred]
		intb <- coefs[intterm]
		vcovs <- vcov(model)
		vcov_pred <- vcovs[pred, pred]
		vcov_int <- vcovs[intterm, intterm]
		vcov_pred_int <- vcovs[pred, intterm]
		range_sequence <- seq(from = modrangeo[1], to = modrangeo[2], 
			by = (modrangeo[2] - modrangeo[1])/1000)
		marginal_effects <- predb + intb * range_sequence
		me_ses <- sqrt(vcov_pred + (range_sequence^2) * vcov_int + 
				2 * range_sequence * vcov_pred_int)
		ts <- marginal_effects/me_ses
		df <- df.residual(model)
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
	if (is.null(vmat)) {
		vmat <- vcov(model)
		covy3 <- vmat[intterm, intterm]
		covy1 <- vmat[pred, pred]
		covy1y3 <- vmat[intterm, pred]
		y3 <- coefs[intterm]
		y1 <- coefs[pred]
	}
	else {
		covy3 <- vmat[intterm, intterm]
		covy1 <- vmat[pred, pred]
		covy1y3 <- vmat[intterm, pred]
		y3 <- coef(model)[intterm]
		y1 <- coef(model)[pred]
	}
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
		return(inside) # added this to see whether J-N region falls inside or outside J-N mod values
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
		return(inside)
	}
	out <- structure(out, inside = inside, failed = failed, all_sig = all_sig)
	cbso1 <- cbs[cbs[, modx] < bounds[1], ]
	cbso2 <- cbs[cbs[, modx] > bounds[2], ]
	cbsi <- cbs[(cbs[, modx] > bounds[1] & cbs[, modx] < bounds[2]), 
		]
	alpha <- alpha * 2
	alpha <- gsub("0\\.", "\\.", as.character(alpha))
	pmsg <- paste("p <", alpha)
	class(out) <- "johnson_neyman"
	return(out)
}



# adapts johnson-neyman from jtools by using vec as the moderator values and 
# testing slopes/cbs across these values
johnson_neyman_veccb2 <- function (model, pred, modx, vmat = NULL, alpha = 0.05, plot = TRUE, 
	control.fdr = FALSE, line.thickness = 0.5, df = "residual", 
	digits = 10, critical.t = NULL, 
	sig.color = "#00BFC4", insig.color = "#F8766D", mod.range = NULL, 
	title = "Johnson-Neyman plot",
	vec) 
{
	predt <- as.character(substitute(pred))
	modxt <- as.character(substitute(modx))
	if (!(predt %in% as.character(attr(terms(model), "variables")))) {
		pred <- as.character(eval(pred))
		modx <- as.character(eval(modx))
	}
	else {
		pred <- predt
		modx <- modxt
	}
	if (df == "residual") {
		df <- df.residual(model)
	}
	else if (df == "normal") {
		df <- Inf
	}
	else if (!is.numeric(df)) {
		stop("df argument must be 'residual', 'normal', or a number.")
	}
	out <- list()
	out <- structure(out, pred = pred, modx = modx, alpha = alpha, 
		plot = plot, digits = digits, control.fdr = control.fdr)
	if (inherits(model, "merMod")) {
		intterm1 <- paste(pred, ":", modx, sep = "")
		intterm1tf <- any(intterm1 %in% names(lme4::fixef(model)))
		intterm2 <- paste(modx, ":", pred, sep = "")
		intterm2tf <- any(intterm2 %in% names(lme4::fixef(model)))
		coefs <- lme4::fixef(model)
	}
	else {
		intterm1 <- paste(pred, ":", modx, sep = "")
		intterm1tf <- any(intterm1 %in% names(coef(model)))
		intterm2 <- paste(modx, ":", pred, sep = "")
		intterm2tf <- any(intterm2 %in% names(coef(model)))
		coefs <- coef(model)
	}
	inttermstf <- c(intterm1tf, intterm2tf)
	intterms <- c(intterm1, intterm2)
	intterm <- intterms[which(inttermstf)]
	modrange <- range(model.frame(model)[, modx])
	modrangeo <- range(model.frame(model)[, modx])
	modsd <- sd(model.frame(model)[, modx])
	if (is.null(mod.range)) {
		modrange[1] <- modrange[1] - modsd 
		modrange[2] <- modrange[2] + modsd
	}
	else {
		modrange <- mod.range
	}
	if (modrange[1] >= modrangeo[1] & modrange[2] <= modrangeo[2]) {
		no_range_line <- TRUE
	}
	else {
		no_range_line <- FALSE
	}
	alpha <- alpha/2
	if (control.fdr == FALSE) {
		tcrit <- qt(alpha, df = df)
		tcrit <- abs(tcrit)
	}
	else if (control.fdr == TRUE) {
		predb <- coefs[pred]
		intb <- coefs[intterm]
		vcovs <- vcov(model)
		vcov_pred <- vcovs[pred, pred]
		vcov_int <- vcovs[intterm, intterm]
		vcov_pred_int <- vcovs[pred, intterm]
		# range_sequence <- seq(from = modrangeo[1], to = modrangeo[2], 
		#     by = (modrangeo[2] - modrangeo[1])/length(vec))
		range_sequence <- vec
		marginal_effects <- predb + intb * range_sequence
		me_ses <- sqrt(vcov_pred + (range_sequence^2) * vcov_int + 
				2 * range_sequence * vcov_pred_int)
		ts <- marginal_effects/me_ses
		df <- df.residual(model)
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
	if (is.null(vmat)) {
		vmat <- vcov(model)
		covy3 <- vmat[intterm, intterm]
		covy1 <- vmat[pred, pred]
		covy1y3 <- vmat[intterm, pred]
		y3 <- coefs[intterm]
		y1 <- coefs[pred]
	}
	else {
		covy3 <- vmat[intterm, intterm]
		covy1 <- vmat[pred, pred]
		covy1y3 <- vmat[intterm, pred]
		y3 <- coef(model)[intterm]
		y1 <- coef(model)[pred]
	}
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
	options(digits = 10)
	# x2 <- seq(from = modrange[1], to = modrange[2], length.out = length(vec))
	x2 <- vec
	predl <- paste("Slope of", pred)
	cbs <- cbands(x2, y1, y3, covy1, covy3, covy1y3, tcrit, predl, 
		modx)
	out$bounds <- bounds
	out <- structure(out, modrange = modrangeo)
	sigs <- which((cbs$Lower < 0 & cbs$Upper < 0) | (cbs$Lower > 
			0 & cbs$Upper > 0))
	insigs <- setdiff(1:length(vec), sigs)
	cbs$Significance <- rep(NA, nrow(cbs))
	cbs$Significance <- factor(cbs$Significance, levels = c("Insignificant", 
		"Significant"))
	index <- 1:length(vec) %in% insigs
	cbs$Significance[index] <- "Insignificant"
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
	plot <- plot + ggplot2::scale_linetype_discrete(guide = ggplot2::guide_legend(order = 1))
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
		x = modx, y = predl) + ggplot2::scale_color_manual(values = c(Significant = sig.color, 
			Insignificant = insig.color), guide = "none") + theme_apa(legend.pos = "right", 
				legend.font.size = 10) + ggplot2::theme(legend.key.size = ggplot2::unit(1, 
					"lines"))
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
		sidebarPanel(width = 3,
			fileInput(inputId = "uploadFile", # input: select a file
				label = "Choose File", multiple = FALSE,
				accept = c("text/csv/text/comma-separated-values,text/solid")),
			#	hr(),
			selectInput(inputId = "y", 
				label = "Select Y (Dependent) Variable", choices = "--"),
			radioButtons(inputId = "std.y",
				label = NULL,
				choices = c(Center = "center", Standardize = "standardize", Raw = "raw"),
				selected = "center", inline = TRUE),
			hr(),
			selectInput(inputId = "x1", 
				label = "Select X1 (Predictor) Variable", choices = "--"),
			radioButtons(inputId = "std.x1",
				label = NULL,
				choices = c(Center = "center", Standardize = "standardize", Raw = "raw"),
				selected = "center", inline = TRUE),
			hr(),
			selectInput(inputId = "x2", 
				label = "Select X2 (Predictor) Variable", choices = "--"),
			radioButtons(inputId = "std.x2",
				label = NULL,
				choices = c(Center = "center", Standardize = "standardize", Raw = "raw"),
				selected = "center", inline = TRUE),
			hr(),
			selectInput(inputId = "cov", 
				label = "Select Covariate(s) (If none, leave blank)", choices = "--", multiple = TRUE),
			htmlOutput(outputId = "messageNonNum"),
			htmlOutput(outputId = "messageSameVar"),
			htmlOutput(outputId = "messageInt"),
			br(),
			htmlOutput(outputId = "plotfeatTitle"),
			fluidRow(
				column(4, uiOutput(outputId = "scatterUI")),
				column(4, uiOutput(outputId = "regplaneUI")),
				column(4, uiOutput(outputId = "predCIUI"))),
			htmlOutput(outputId = "mx1CheckboxTitle"),
			fluidRow(
				column(6, uiOutput(outputId = "pm.mx1.wROSUI")),
				column(4, uiOutput(outputId = "coonx1UI"))),
			htmlOutput(outputId = "mx2CheckboxTitle"),
			fluidRow(
				column(6, uiOutput(outputId = "pm.mx2.wROSUI")),
				column(4, uiOutput(outputId = "coonx2UI"))),
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
			navbarPage(title = "",
				tabPanel(title = "How To",
					htmlOutput(outputId = "guide")),
				tabPanel(title = "Main Output", 
					fluidRow(
						column(7,
							plotlyOutput(outputId = "plot3dPlot")),
						column(5,
							htmlOutput(outputId = "modFit"),
							inlineCSS(list("table" = "font-size: 10px")),
							htmlOutput("jnmx1.sum"),
							htmlOutput("jnmx2.sum"))),
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
					tabPanel(title = "Slope Confidence Band Tables",
						htmlOutput(outputId = "mx1cbTitle"),
						htmlOutput(outputId = "jnmx1Text"),
						DT::dataTableOutput("jnmx1cbInfo"),
						htmlOutput(outputId = "mx2cbTitle"),
						htmlOutput(outputId = "jnmx2Text"),
						DT::dataTableOutput("jnmx2cbInfo")),
					tabPanel(title = "Marginal Effect Plots",
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
						br())),
				tabPanel(title = "Descriptives", 
					htmlOutput(outputId = "desc.stdTitle"),
					tableOutput(outputId = "desc.std"),
					htmlOutput(outputId = "descTitle"),
					tableOutput(outputId = "desc")),
				tabPanel(title = "Raw Data", DT::dataTableOutput(outputId = "rawdata")),
				tabPanel(title = "Code for Saving Plots", 
					htmlOutput(outputId = "codetitle1"))))
	)
)

server <- function(input, output, session) {
	output$guide <- renderUI(
		HTML(paste(
			"<div style='font-size: 12px'> 
			<h5> Welcome </h5>
					This visualization tool represents continuous by continuous interactions as a regression plane in 
					3D space in combination with the regions of significance from the Johnon-Neyman technique. You can see an example of the type of plot produced <a href=http://rpubs.com/sbu_mfinsaas/Example target=_blank> here</a>. 
					Please contact Megan Finsaas at megan.finsaas@gmail.com if you have questions or comments about the app. <br> <br>
					For the best performance, please run the app in <b> Google Chrome</b>. 
			<h5> Uploading a Datafile </h5>
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
<h5> Inputting Variables for Regression Model </h5>
					<ul>
						<li> Click the dropdown buttons for variable names or delete “--” and type the variable name. </li>
						<li> Once the predictors and an outcome are input, the regression model will automatically estimate.
						<li> Predictors and the outcome are centered by default. You can alternatively standardize them or keep them in the raw form. The model will automatically re-estimate whenever the standardization option is changed. </li>
						<li> To add covariates, click the text box for a dropdown list or type the variable names. Covariates are always centered at their means. They can be entered at any time and the model will automatically update to include them.</li>
						<li> At this point, the program does not accept squared terms. </li>			
					</ul>

			<h5> Main Output Tab </h5>
					<ul>
						<li> Here you will find the <b>model summary</b> and <b>3D plot</b>. </li>
						<li> When	a region of significance falls within the observed data, the <b>Johnson-Neyman output</b> will include the values
							on the moderator that mark the bounds of the region of significance, the range of slope estimates for the 
							predictor-outcome relationship and the percentage of cases within this region, and the observed range of the moderator variable. </li>
						<li> When a <b>crossover point</b> falls within the observed data, it will also print on this tab. </li>
					</ul>
			<h5> Plot Features </h5>
					<ul> 
						<li> The Plotly buttons at the top right of the plot allow you to <b> zoom, pan, rotate, and save a .png </b> of the plot.  </li>
						<li> For all models, checkboxes for <b> adding the scatterplot, regression plane, </b>and <b>95% confidence interval</b> around the predicted values will appear in the left navigation bar. </li>
						<li>	Which other plot feature options appear below these depend on whether there are regions of significance or crossover points within the observed range of data. </li>
							<ul>
								<li> If a region of significance falls within the range of data, the option to <b> shade the region of significance </b> on the regression plane will appear. Selecting the “Solid ROS” button will shade the region of significance on the plot in one color; selecting the “Gradient ROS” button will shade the region using a gradient that reflects the width of the 95% confidence band around the slope estimate.
								You can also view the confidence band estimates under the Johnson-Neyman tab. </li>
								<li> If a crossover point falls within the range of data, the “Crossover” checkbox will appear, which <b> adds the crossover point </b>to the figure.</li>
							</ul>
						</ul>
			<h5> Additional Output </h5>
					<ul>
						<li> <b> Johnson-Neyman Tab: </b> Marginal effects plots and Johnson-Neyman bounds regardless of whether they fall within the range of observed data; tables with slope estimates, 95% confidence bands and <i>p</i>-values at all values of the moderator.</li> 
						<li> <b> Descriptives Tab: </b> Table of descriptives for variables in model with standardization options applied; table for all variables in the dataset in raw form </li>
						<li> <b> Raw Data Tab: </b> Sortable data frame of raw data. </li>
					</ul>
			<h5> Saving and Sharing Plots </h5>
					There are three ways to save and share your plots:
					<ol>
						<li> Take screen/snapshots of the figure in various rotations (i.e., PrtScn on Windows; Shift + Command + 4 on Macs) or click the Plotly camera button to save a .png file. </li>
						<li> Record the plot while manually rotating the plot. To keep the mouse outside the plot while recording, click and drag above it. </li>
						<li> To create a “live” rotatable plot with a permanent link:
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
		DT::datatable(data = df(),
			options = list(pageLength = 10, lengthMenu = c(10, 25, 40)),
			rownames = FALSE)
	})
	
	# determine whether variables are numeric ------------
	non.numeric <- reactive({
		req(input$x1 != "--", input$x2 != "--", input$y != "--")
		mat <- as.matrix(df()[, c(input$y, input$x1, input$x2, input$cov)])
		if(is.numeric(mat) == TRUE) {
			non.numeric <- FALSE
		} else {
			non.numeric <- TRUE
		}
	})
	
	
	# determine whether same var input ----------------------------------------
	sameVar <- reactive({
		req(input$x1 != "--", input$x2 != "--", input$y != "--")
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
	
	# standardize variables and subset data --------
	df.std <- eventReactive(c(input$x1, input$x2, input$y, input$std.x1, input$std.x2, input$std.y, input$cov, input$uploadFile), {
		req(input$x1 != '--', input$x2 != '--', input$y != '--', sameVar() == FALSE, non.numeric() == FALSE)
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
		descdf <- describe(df())[c("n","mean","sd","median","min","max","skew","kurtosis")]
		descdf <- as.data.frame(descdf)
		Variable <- colnames(df())
		cbind(Variable, descdf)
	})
	
	output$desc.std <- renderTable({ # descriptives on variables in model with standardization option applied
		descdf <- describe(df.std())[c("n","mean","sd","median","min","max","skew","kurtosis")]
		descdf <- as.data.frame(descdf)
		Variable <- colnames(df.std())
		cbind(Variable, descdf)
	})
	
	# render descriptives output
	output$descTitle <- renderUI({
		req(input$uploadFile)
		h6("Descriptive Statistics on Raw Data")
	})
	
	output$desc.stdTitle <- renderUI({
		req(input$uploadFile, input$x1 != '--', input$x2 != '--', input$y != '--', sameVar() == FALSE, non.numeric() == FALSE)
		h6("Descriptive Statistics on Variables in Model with Standardization Options Applied")
	})
	
	
	# estimate linear model ---------------------------------------------------
	mod <- reactive({
		lm(paste(input$y, " ~ . + ", input$x1, " * ", input$x2), data = df.std())
	})
	
	modsum <- reactive({
		summary(mod())
	})
	
	# model output ------------------------------------------------------------
	output$modInfo <- renderUI({
		req(input$x1 != '--', input$x2 != '--', input$y != '--', sameVar() == FALSE, non.numeric() == FALSE)
		HTML(paste("<div style='font-weight: 700' >Model Summary</div>"))
	})
	
	output$modFit <- renderUI({
		req(input$x1 != '--', input$x2 != '--', input$y != '--', sameVar() == FALSE, non.numeric() == FALSE)
		dv <- paste0("Dependent variable: ", input$y)
		
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
	})
	
	modTable <- renderTable({
		req(input$x1 != '--', input$x2 != '--', input$y != '--', sameVar() == FALSE, non.numeric() == FALSE)
		modtable <- unlist(modsum()[4])
		modtable <- data.frame(matrix(modtable, ncol = 4))
		modtable <- format(round(modtable, digits = 3), nsmall = 3)
		modtable <- cbind(names(coef(mod())), modtable)
		colnames(modtable) <- c(" ", "Est.","SE", "t-value","p-value")
		return(modtable)
	})
	
	# test whether interaction term is significant ----------------------------
	intsig <- reactive({
		modsum <- summary(mod())
		modsum$coef[nrow(modsum$coef), 4] < .05 # access significance column for interaction term
	})
	
	
	# identify crossover points -----------------------------------------------
	ncoef <- reactive({
		length(coef(mod()))
	})
	
	# x1 as mod: -Bx1/Bx1*x2, plot on x2
	coonx2 <- reactive({
		mod <- mod()
		coonx2 <- -coef(mod)[2]/coef(mod)[ncoef()]
		if(coonx2 >= minx2() & coonx2 <= maxx2()) {
			coonx2 <- as.numeric(coonx2)
		} else {
			coonx2 <- NA
		}
	})
	
	# predicted y at crossover
	coonx2y <- reactive({
		data <- df.std()
		if (!is.na(coonx2())) {
			if (is.null(input$cov)) {
				newdatacoonx2y <- data.frame(vecx1(), coonx2())
				names(newdatacoonx2y) <- c(input$x1, input$x2)
			} else {
				if (length(input$cov) == 1) {
					x <- mean(data[,4])
					pmcov.l <- list(x) 
				} else {
					pmcov.l <- sapply(c(data[,4:ncol(data)]), mean, simplify = FALSE)
				}
				newdatacoonx2y <- data.frame(vecx1(), coonx2(), pmcov.l)
				names(newdatacoonx2y) <- c(input$x1, input$x2, input$cov)
			}
			predict.lm(mod(), newdatacoonx2y)
		} else {
			coonx2y <- NA
		}
	})
	
	# x2 as mod: -Bx2/Bx1*x2, plot on x1
	coonx1 <- reactive({
		mod <- mod()
		coonx1 <- -coef(mod)[3]/coef(mod)[ncoef()]
		if(coonx1 >= minx1() & coonx1 <= maxx1()) {
			coonx1 <- as.numeric(coonx1)
		} else {
			coonx1 <- NA
		}
	})
	
	# predicted y at crossover
	coonx1y <- reactive({
		data <- df.std()
		if (!is.na(coonx1())) {
			if (is.null(input$cov)) {
				newdatacoonx1y <- data.frame(coonx1(), vecx2())
				names(newdatacoonx1y) <- c(input$x1, input$x2)
			} else {
				if (length(input$cov) == 1) {
					x <- mean(data[,4])
					pmcov.l <- list(x) 
				} else {
					pmcov.l <- sapply(c(data[,4:ncol(data)]), mean, simplify = FALSE)
				}
				newdatacoonx1y <- data.frame(coonx1(), vecx2(), pmcov.l)
				names(newdatacoonx1y) <- c(input$x1, input$x2, input$cov)
			}
			predict.lm(mod(), newdatacoonx1y)
		} else {
			coonx1y <- NA
		}
	})
	
	
	# min/max of x1/x2 --------------------------------------------------------
	minx1 <- reactive({
		data <- df.std()
		min(data[,input$x1], na.rm = TRUE)
	})
	
	maxx1 <- reactive({
		data <- df.std()
		max(data[,input$x1], na.rm = TRUE)
	})
	
	minx2 <- reactive({
		data <- df.std()
		min(data[,input$x2], na.rm = TRUE)
	})
	
	maxx2 <- reactive({
		data <- df.std()
		max(data[,input$x2], na.rm = TRUE)
	})
	
	# johnson neyman analysis -------------------------------------------------
	# when x1 is mod
	jnmx1 <- reactive({
		johnson_neyman(mod(), pred = input$x2, modx = input$x1,
			sig.color = "#749702", insig.color = "#b3b2b2",
			mod.range = c(minx1(), maxx1()), alpha = 0.05, plot = TRUE)
	})
	
	# save jn values if within mod range
	jnmx1.1 <- reactive({
		if(jnmx1()$bounds[1] >= minx1() & jnmx1()$bounds[1] <= maxx1()) {
			jnmx1.1 <- as.numeric(jnmx1()$bounds[1])
		} else {
			jnmx1.1 <- NA
		}
	})
	
	jnmx1.2 <- reactive({
		if(jnmx1()$bounds[2] >= minx1() & jnmx1()$bounds[2] <= maxx1()) {
			jnmx1.2 <- as.numeric(jnmx1()$bounds[2])
		} else {
			jnmx1.2 <- NA
		}
	})
	
	# when x2 is mod
	jnmx2 <- reactive({
		johnson_neyman(mod(), pred = input$x1, modx = input$x2,
			sig.color = "#3366ff", insig.color = "#b3b2b2",
			mod.range = c(minx2(), maxx2()), alpha = 0.05, plot = TRUE)
	})
	
	# save jn values if within mod range
	jnmx2.1 <- reactive({
		if(jnmx2()$bounds[1] >= minx2() & jnmx2()$bounds[1] <= maxx2()) {
			jnmx2.1 <- as.numeric(jnmx2()$bounds[1])
		} else {
			jnmx2.1 <- NA
		}
	})
	
	jnmx2.2 <- reactive({
		if(jnmx2()$bounds[2] >= minx2() & jnmx2()$bounds[2] <= maxx2()) {
			jnmx2.2 <- as.numeric(jnmx2()$bounds[2])
		} else {
			jnmx2.2 <- NA
		}
	})
	
	# jn confidence band tables -----------------------------------------------
	jnmx1cb <- reactive({
		johnson_neyman_veccb2(mod(), pred = input$x2, modx = input$x1,  alpha = 0.05, 
			mod.range = c(minx1(), maxx1()), plot = FALSE, vec = vecx1())$cbands
	})
	
	jnmx2cb <- reactive({
		johnson_neyman_veccb2(mod(), pred = input$x1, modx = input$x2, alpha = 0.05, 
			mod.range = c(minx2(), maxx2()), plot = FALSE, vec = vecx2())$cbands
	})
	
	# johnson neyman output ---------------------------------------------------
	# marginal effects plots
	output$jnmx1Title <- renderUI({
		req(non.numeric() == FALSE, sameVar() == FALSE)
		h5(paste0(input$x1, " as Moderator"))
	})
	
	output$jnmx1Plot <- renderPlot({
		req(non.numeric() == FALSE, sameVar() == FALSE)
		if(intsig() == TRUE & mx1.ROS.none() == FALSE) {
			jnmx1()$plot
		} 
	})
	
	output$jnmx1Info <- renderPrint({
		req(non.numeric() == FALSE, sameVar() == FALSE)
		if(intsig() == TRUE & mx1.ROS.none() == FALSE) {
			jnmx1()
		} else {
			HTML(paste0("The interaction term is not significant and/or there is no region of significance within the range of ", input$x1,"."))
		}
	})
	
	output$jnmx2Title <- renderUI({
		req(non.numeric() == FALSE, sameVar() == FALSE)
		h5(paste0(input$x2, " as Moderator"))
	})
	
	output$jnmx2Plot <- renderPlot({
		req(non.numeric() == FALSE, sameVar() == FALSE)
		if(intsig() == TRUE & mx2.ROS.none() == FALSE) {
			jnmx2()$plot
		}
	})
	
	output$jnmx2Info <- renderPrint({
		req(non.numeric() == FALSE, sameVar() == FALSE)
		if(intsig() == TRUE & mx2.ROS.none() == FALSE) {
			jnmx2()
		} else {
			HTML(paste0("The interaction term is not significant and/or there is no region of significance within the range of ", input$x2,"."))
		}
	})
	
	
	# confidence band tables
	output$mx1cbTitle <- renderUI({
		req(non.numeric() == FALSE, sameVar() == FALSE)
		h5(paste0(input$x1, " as Moderator"))
	})
	
	output$jnmx1cbInfo <- DT::renderDataTable({
		req(non.numeric() == FALSE, sameVar() == FALSE)
		if(intsig() == TRUE & mx1.ROS.none() == FALSE) {
			jnmx1cb <- jnmx1cb()
			jnmx1cb[,1:4] <- format(round(jnmx1cb[,1:4], digits = 3), nsmall = 3)
			DT::datatable(jnmx1cb, options = list(lengthMenu = c(10,20,50, length(vecx1()))))
		}
	})
	output$jnmx1Text <- renderPrint({
		req(non.numeric() == FALSE, sameVar() == FALSE)
		if(intsig() == FALSE | mx1.ROS.none() == TRUE) {
			HTML(paste0("The interaction term is not significant and/or there is no region of significance within the range of ", input$x1,"."))
		}
	})
	
	output$mx2cbTitle <- renderUI({
		req(non.numeric() == FALSE, sameVar() == FALSE)
		h5(paste0(input$x2, " as Moderator"))
	})
	
	output$jnmx2cbInfo <- DT::renderDataTable({
		req(non.numeric() == FALSE, sameVar() == FALSE)
		if(intsig() == TRUE & mx2.ROS.none() == FALSE) {
			jnmx2cb <- jnmx2cb()
			jnmx2cb[,1:4] <- format(round(jnmx2cb[,1:4], digits = 3), nsmall = 3)
			DT::datatable(jnmx2cb, options = list(lengthMenu = c(10,20,50, length(vecx2()))))
		}
	})
	
	output$jnmx2Text <- renderPrint({
		req(non.numeric() == FALSE, sameVar() == FALSE)
		if(intsig() == FALSE | mx2.ROS.none() == TRUE) {
			HTML(paste0("The interaction term is not significant and/or there is no region of significance within the range of ", input$x2,"."))
		}
	})
	
	# johnson neyman and crossover text output for main tab ---------------------------------------------------------
	output$jnmx1.sum <- renderUI({
		req(input$x1 != '--', input$x2 != '--', input$y != '--', sameVar() == FALSE, non.numeric() == FALSE)
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
			obsorjnmx1 <- c(" (minimum)", " (maximum)", " (J-N bound)")
			
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
						jnmx1.prop <- prop.table(table(df.std()[,input$x1] >= jnmx1.2()))
						jnmx1.prop <- jnmx1.prop[dim(jnmx1.prop)]
						jnmx1.sum <- list(jnmx1.sl[2], jnmx1.slest[2], jnmx1.slestminmax[2], jnmx1.sign[1], jnmx1.val[2],
							jnmx1.sign[2], minmaxx1[2], obsorjnmx1[3], obsorjnmx1[2])
					} else {
						if(!is.na(jnmx1.1()) & is.na(jnmx1.2())) { # low band / outside
							jnmx1.prop <- prop.table(table(df.std()[,input$x1] <= jnmx1.1()))
							jnmx1.prop <- jnmx1.prop[dim(jnmx1.prop)]
							jnmx1.sum <- list(jnmx1.sl[1], jnmx1.slestminmax[1], jnmx1.slest[1], jnmx1.sign[1], minmaxx1[1],
								jnmx1.sign[2], jnmx1.val[1], obsorjnmx1[1], obsorjnmx1[3])
						} else { # high and low band / outside
							if(!is.na(jnmx1.1()) & !is.na(jnmx1.2())) {
								jnmx1.prop <- prop.table(table(df.std()[,input$x1] >= jnmx1.1() & 
										df.std()[,input$x1] <= jnmx1.2()))
								jnmx1.prop <- jnmx1.prop[dim(jnmx1.prop)]
								jnmx1.sum <- list(jnmx1.sl[1], jnmx1.sign[1], minmaxx1[1], jnmx1.sign[2], jnmx1.val[1], 
									jnmx1.slestminmax[1], jnmx1.slest[1], 
									jnmx1.sl[2], jnmx1.sign[1], jnmx1.val[2], jnmx1.sign[2], minmaxx1[2], 
									jnmx1.slest[2], jnmx1.slestminmax[2], obsorjnmx1[1], obsorjnmx1[3], obsorjnmx1[3], obsorjnmx1[2])
							}
						}
					}
				} else { # ROS falls inside JN values
					if(!is.na(jnmx1.1()) & !is.na(jnmx1.2())) { # (no squared terms allowed, center band not run)
						jnmx1.prop <- prop.table(table(df.std()[,input$x1] >= jnmx1.1() & 
								df.std()[,input$x1] <= jnmx1.2()))
						jnmx1.prop <- jnmx1.prop[dim(jnmx1.prop)]
						jnmx1.sum <- list(jnmx1.sl, jnmx1.slest[1], jnmx1.slest[2], 
							jnmx1.sign[1], jnmx1.val[1], jnmx1.sign[2], minmaxx1[2], obsorjnmx1[3], obsorjnmx1[3])
					} else {
						if(is.na(jnmx1.1()) & !is.na(jnmx1.2())) { # low band / inside
							jnmx1.prop <- prop.table(table(df.std()[,input$x1] <= jnmx1.2()))
							jnmx1.prop <- jnmx1.prop[dim(jnmx1.prop)]
							jnmx1.sum <- list(jnmx1.sl[2], jnmx1.slestminmax[1], jnmx1.slest[2], jnmx1.sign[1], minmaxx1[1],
								jnmx1.sign[2], jnmx1.val[2], obsorjnmx1[1], obsorjnmx1[3])
						} else {
							if (!is.na(jnmx1.1()) & is.na(jnmx1.2())) { # high band / inside
								jnmx1.prop <- prop.table(table(df.std()[,input$x1] >= jnmx1.1()))
								jnmx1.prop <- jnmx1.prop[dim(jnmx1.prop)]
								jnmx1.sum <- list(jnmx1.sl[1], jnmx1.slest[1], jnmx1.slestminmax[2], jnmx1.sign[1], jnmx1.val[1], 
									jnmx1.sign[2], minmaxx1[2], obsorjnmx1[3], obsorjnmx1[2])
							} 
						}
					}
				}
			}
			
				if(mx1.ROS.all() == TRUE) {
					jnmx1.text <- paste0("The relationship between ", input$x2, " and ", input$y,
						" is significantly ", jnmx1.sum[1], " across the entire range of ", input$x1, " (minimum: ",
						jnmx1.sum[4], "; maximum: ", jnmx1.sum[5], "; slope estimate range: ", jnmx1.sum[2], " to " , jnmx1.sum[3],
						"). This region of significance includes 100.00% of the sample.")
				} else {
					jnmx1.perc <- jnmx1.prop*100
					jnmx1.perc <- format(round(jnmx1.perc, digits = 2), nsmall = 2)
					if (jnmx1.inside() == FALSE & !is.na(jnmx1.1()) & !is.na(jnmx1.2())) { # two ROS (high and low) 
						jnmx1.text <- paste0("The relationship between ", input$x2, " and ", input$y,
							" is significantly ", jnmx1.sum[1], " when ", input$x1,
							jnmx1.sum[2], jnmx1.sum[3], jnmx1.sum[15], " and", jnmx1.sum[4], jnmx1.sum[5], jnmx1.sum[16],
							" (slope estimate range: = ", jnmx1.sum[6], " to ", jnmx1.sum[7], ") and significantly ", jnmx1.sum[8], " when ",
							input$x1, jnmx1.sum[9], jnmx1.sum[10],jnmx1.sum[17], " and ", jnmx1.sum[11], jnmx1.sum[12], jnmx1.sum[18],
							" (slope estimate range: ", jnmx1.sum[13], " to ", jnmx1.sum[14], 
							"). These regions of significance include ", jnmx1.perc, "% of the sample. ")
					} else { # one ROS
						jnmx1.text <- paste0("When ", input$x1,	jnmx1.sum[4], jnmx1.sum[5], jnmx1.sum[8], " and", jnmx1.sum[6], jnmx1.sum[7], jnmx1.sum[9],
							", the relationship between ", input$x2, " and ", input$y,
							" is significantly ", jnmx1.sum[1], " (slope estimate range: ", jnmx1.sum[2], 
							" to ", jnmx1.sum[3], "). This region of significance includes ", jnmx1.perc, "% of the sample. ")
					}
				}
		} else {
			jnmx1.text <- paste0("The relationship between ", input$x2, " and ", input$y, " is not significantly different from zero ",
				"at any value of ", input$x1, ".")
		}
		
		if (!is.na(coonx1())) {
			coonx1 <- format(round(coonx1(), digits = 2), nsmall = 2)
			coonx1.text <- paste0("Crossover point: ", coonx1)
		} else {
			coonx1.text <- paste0("The crossover point falls outside the range of ", input$x1, ".")
		}
		
		HTML(paste0("<div class='well' >", "<label>", input$x1, " as Moderator</label> <br/>",
			jnmx1.text, "<br/> <br/>", coonx1.text, "</div>"))
	})
	
	output$jnmx2.sum <- renderUI({
		req(input$x1 != '--', input$x2 != '--', input$y != '--', sameVar() == FALSE, non.numeric() == FALSE)
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
			obsorjnmx2 <- c(" (minimum)", " (maximum)", " (J-N bound)")
			
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
						jnmx2.prop <- prop.table(table(df.std()[,input$x2] >= jnmx2.2()))
						jnmx2.prop <- jnmx2.prop[dim(jnmx2.prop)]
						jnmx2.sum <- list(jnmx2.sl[2], jnmx2.slest[2], jnmx2.slestminmax[2], jnmx2.sign[1], jnmx2.val[2],
							jnmx2.sign[2], minmaxx2[2], obsorjnmx2[3], obsorjnmx2[2])
					} else {
						if(!is.na(jnmx2.1()) & is.na(jnmx2.2())) { # low band / outside
							jnmx2.prop <- prop.table(table(df.std()[,input$x2] <= jnmx2.1()))
							jnmx2.prop <- jnmx2.prop[dim(jnmx2.prop)]
							jnmx2.sum <- list(jnmx2.sl[1], jnmx2.slestminmax[1], jnmx2.slest[1], jnmx2.sign[1], minmaxx2[1],
								jnmx2.sign[2], jnmx2.val[1], obsorjnmx2[1], obsorjnmx2[3])
						} else { # high and low band / outside
							if(!is.na(jnmx2.1()) & !is.na(jnmx2.2())) {
								jnmx2.prop <- prop.table(table(df.std()[,input$x2] >= jnmx2.1() & 
										df.std()[,input$x2] <= jnmx2.2()))
								jnmx2.prop <- jnmx2.prop[dim(jnmx2.prop)]
								jnmx2.sum <- list(jnmx2.sl[1], jnmx2.sign[1], minmaxx2[1], jnmx2.sign[2], jnmx2.val[1], 
									jnmx2.slestminmax[1], jnmx2.slest[1], 
									jnmx2.sl[2], jnmx2.sign[1], jnmx2.val[2], jnmx2.sign[2], minmaxx2[2], 
									jnmx2.slest[2], jnmx2.slestminmax[2], obsorjnmx2[1], obsorjnmx2[3], obsorjnmx2[3], obsorjnmx2[2])
							}
						}
					}
				} else { # ROS falls inside JN values
					if(!is.na(jnmx2.1()) & !is.na(jnmx2.2())) { # (no squared terms allowed, center band not run)
						jnmx2.prop <- prop.table(table(df.std()[,input$x2] >= jnmx2.1() & 
								df.std()[,input$x2] <= jnmx2.2()))
						jnmx2.prop <- jnmx2.prop[dim(jnmx2.prop)]
						jnmx2.sum <- list(jnmx2.sl[1], jnmx2.slest[1], jnmx2.slest[2], 
							jnmx2.sign[1], jnmx2.val[1], jnmx2.sign[2], minmaxx2[2], obsorjnmx2[3], obsorjnmx2[3])
					} else {
						if(is.na(jnmx2.1()) & !is.na(jnmx2.2())) { # low band / inside
							jnmx2.prop <- prop.table(table(df.std()[,input$x2] <= jnmx2.2()))
							jnmx2.prop <- jnmx2.prop[dim(jnmx2.prop)]
							jnmx2.sum <- list(jnmx2.sl[2], jnmx2.slestminmax[1], jnmx2.slest[2], jnmx2.sign[1], minmaxx2[1],
								jnmx2.sign[2], jnmx2.val[2], obsorjnmx2[1], obsorjnmx2[3])
						} else {
							if (!is.na(jnmx2.1()) & is.na(jnmx2.2())) { # high band / inside
								jnmx2.prop <- prop.table(table(df.std()[,input$x2] >= jnmx2.1()))
								jnmx2.prop <- jnmx2.prop[dim(jnmx2.prop)]
								jnmx2.sum <- list(jnmx2.sl[1], jnmx2.slest[1], jnmx2.slestminmax[2], jnmx2.sign[1], jnmx2.val[1], 
									jnmx2.sign[2], minmaxx2[2], obsorjnmx2[3], obsorjnmx2[2])
							}
						}
					}
				}
			}
			
			if(mx2.ROS.all() == TRUE) {
				jnmx2.text <- paste0("The relationship between ", input$x1, " and ", input$y,
					" is significantly ", jnmx2.sum[1], " across the entire range of ", input$x2, " (minimum: ",
					jnmx2.sum[4], "; maximum: ", jnmx2.sum[5], "; slope estimate range: ", jnmx2.sum[2], " to " , jnmx2.sum[3],
					"). This region of significance includes 100.00% of the sample.")
			} else {
				jnmx2.perc <- jnmx2.prop*100
				jnmx2.perc <- format(round(jnmx2.perc, digits = 2), nsmall = 2)
				if (jnmx2.inside() == FALSE & !is.na(jnmx2.1()) & !is.na(jnmx2.2())) { # two ROS (high and low) # NEED TO FIX
					jnmx2.text <- paste0("The relationship between ", input$x1, " and ", input$y,
						" is significantly ", jnmx2.sum[1], " when ", input$x2,
						jnmx2.sum[2], jnmx2.sum[3], jnmx2.sum[15], " and", jnmx2.sum[4], jnmx2.sum[5], jnmx2.sum[16],
						" (slope estimate range: = ", jnmx2.sum[6], " to ", jnmx2.sum[7], ") and significantly ", jnmx2.sum[8], " when ",
						input$x2, jnmx2.sum[9], jnmx2.sum[10],jnmx2.sum[17], " and ", jnmx2.sum[11], jnmx2.sum[12], jnmx2.sum[18],
						" (slope estimate range: ", jnmx2.sum[13], " to ", jnmx2.sum[14], 
						"). These regions of significance include ", jnmx2.perc, "% of the sample. ")
				} else {
					jnmx2.text <- paste0("When ", input$x2,	jnmx2.sum[4], jnmx2.sum[5], jnmx2.sum[8], " and", jnmx2.sum[6], jnmx2.sum[7], jnmx2.sum[9],
						", the relationship between ", input$x1, " and ", input$y,
						" is significantly ", jnmx2.sum[1], " (slope estimate range: ", jnmx2.sum[2], 
						" to ", jnmx2.sum[3], "). This region of significance includes ", jnmx2.perc, "% of the sample. ")
				}
			}
		} else {
			jnmx2.text <- paste0("The relationship between ", input$x1, " and ", input$y, " is not significantly different from zero ",
				"at any value of ", input$x2, ".")
		}
		
		if (!is.na(coonx2())) {
			coonx2 <- format(round(coonx2(), digits = 2), nsmall = 2)
			coonx2.text <- paste0("Crossover point: ", coonx2)
		} else {
			coonx2.text <- paste0("The crossover point falls outside the range of ", input$x2, ".")
		}
		
		HTML(paste0("<div class='well' >", "<label>", input$x2, " as Moderator</label> <br/>",
			jnmx2.text, "<br/> <br/>", coonx2.text, "</div>"))
	})
	
	# vectors of x1 and x2 -----------------------------------
	vecx1 <- reactive({
		vecx1 <- seq(minx1(), maxx1(), length.out = 100)
		vecx1 <- c(vecx1, jnmx1.1(), jnmx1.2()) # add in j-n values 
		vecx1 <- sort(vecx1)
	})
	
	vecx2 <- reactive({
		vecx2 <- seq(minx2(), maxx2(), length.out = 100)
		vecx2 <- c(vecx2, jnmx2.1(), jnmx2.2())	# add in j-n values 
		vecx2 <- sort(vecx2)
	})
	
	# predicted values matrix -------------------------------------------------
	pm <- reactive({
		pm <- outer(vecx1(), vecx2(), function(pmx1, pmx2) {
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
			predict.lm(mod(), newdatapm)
		})
		rownames(pm) <- vecx1()
		colnames(pm) <- vecx2()
		pm
	})
	
	# lower ci of predicted values matrix -------------------------------------------------
	pmcil <- reactive({
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
	pmciu <- reactive({
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
	# identify whether ROS is inside j-n values
	jnmx1.inside <- reactive({
		johnson_neyman_inside(mod(), pred = input$x2, modx = input$x1, 
			mod.range = c(minx1(), maxx1()), alpha = 0.05)
	})
	
	# identify whether ROS covers whole range of x1
	mx1.ROS.all <- reactive({
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
	mx1.ROS.none <- reactive({
		if((is.na(jnmx1.1()) & is.na(jnmx1.2())) & mx1.ROS.all() == FALSE) {
			mx1.ROS.none <- TRUE
			mx1.ROS.none
		} else {
			mx1.ROS.none <- FALSE
			mx1.ROS.none
		}
	})
	
	# identify row number of jn values
	jnmx1.1rn <- reactive({
		jnmx1.1rn <- match(jnmx1.1(), vecx1())
	})
	
	jnmx1.2rn <- reactive({
		jnmx1.2rn <- match(jnmx1.2(), vecx1())
	})
	
	# make x1 ROS planes
	pm.mx1.wROS <- reactive({
		if (mx1.ROS.none() == FALSE & intsig() == TRUE) {
			
			pm <- pm()
			jnmx1.1 <- jnmx1.1()
			jnmx1.2 <- jnmx1.2()
			vecx1 <- vecx1()
			vecx2 <- vecx2()
			
			# ROS covers whole range
			if(mx1.ROS.all()==TRUE) {
				pm.mx1.wROS <- pm() 
			}
			
			# ROS inside JN bounds
			if(jnmx1.inside() == TRUE) {
				# inside both j-n values (center band) (not run, program does not accept squared terms)
				if(!is.na(jnmx1.1) & !is.na(jnmx1.2)) {
					pm.mx1.ROS1i <- pm[vecx1 >= jnmx1.1 & vecx1 <= jnmx1.2,]
					mx1.pad1 <- matrix(nrow = jnmx1.1rn(), ncol = ncol(pm)) # lower pad made of NAs
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
	# identify whether ROS is inside j-n values
	jnmx2.inside <- reactive({
		johnson_neyman_inside(mod(), pred = input$x1, modx = input$x2, 
			mod.range = c(minx2(), maxx2()), alpha = 0.05)
	})
	
	# identify whether ROS covers whole range of x2
	mx2.ROS.all <-reactive({
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
	mx2.ROS.none <- reactive({
		if((is.na(jnmx2.1()) & is.na(jnmx2.2())) & mx2.ROS.all() == FALSE) {
			mx2.ROS.none <- TRUE
			mx2.ROS.none
		} else {
			mx2.ROS.none <- FALSE
			mx2.ROS.none
		}
	})
	
	# identify column number of jn bounds
	jnmx2.1cn <- reactive({
		jnmx2.1cn <- match(jnmx2.1(), vecx2())
	})
	
	jnmx2.2cn <- reactive({
		jnmx2.2cn <- match(jnmx2.2(), vecx2())
	})
	
	
	# make x2 ROS planes
	pm.mx2.wROS <- reactive({
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
				# inside both j-n values (center band) (not run, program does not accept squared terms)
				if(!is.na(jnmx2.1) & !is.na(jnmx2.2)) {
					pm.mx2.ROS1i <- pm[,vecx2 >= jnmx2.1 & vecx2 <= jnmx2.2]
					mx2.pad1 <- matrix(nrow = nrow(pm), ncol = jnmx2.1cn()) # lower pad made of NAs
					mx2.pad2 <- matrix(nrow = nrow(pm), ncol = ncol(pm)- jnmx2.2cn()) # upper pad  made of NAs
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
		coonx1CB = FALSE, coonx2 = FALSE)
	
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
	
	
	# create 3D plot ---------------------------------------------------------
	plot3d <- reactive({
		req(input$x1 != "--", input$x2 != "--",input$y != "--", sameVar() == FALSE, non.numeric() == FALSE)
		
		data <- df.std()
		pmcolor <- array(1L, dim(pm())) # matrix for regression plane color
		
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
				x = ~data[,input$x1],
				y = ~data[,input$x2],
				z = ~data[,input$y],
				type = "scatter3d",
				mode = "markers",
				marker = list(size = scattersize,
					color = "grey",
					opacity = scatterop)) %>%
			# plot formatting
			layout(
				scene = list(
					camera = list(eye = list(x = 1.6, y = -1.6, z = 1.25)),
					xaxis = list(title = input$x1,
						titlefont = list(family = "Arial, sans-serif", color = "#3366ff", size = 10),
						tickfont = list(family = "Arial, sans-serif", size = 10)),
					yaxis = list(title = input$x2,
						titlefont = list(family = "Arial, sans-serif", color = "#567001", size = 10),
						tickfont = list(family = "Arial, sans-serif", size = 10)),
					zaxis = list(title = input$y,
						titlefont = list(family = "Arial, sans-serif", color = "#484747", size = 10),
						tickfont = list(family = "Arial, sans-serif",	size = 10)))) %>%
			hide_legend()
		
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
					colorscale = list(c(0, 1), c("#3366ff", 'white')),
					opacity = 0.7,
					showscale = FALSE)
			}
			if(input$pm.mx2.wROSCheckbox == "solid") {
				pm.mx2.wROSg <- pm.mx2.wROS()
				pm.mx2.wROSg[!is.na(pm.mx2.wROSg)] <- 1
				p <- add_surface(p, x = ~vecx1(), y = ~vecx2(), z = ~t(pm.mx2.wROS()),
					type = "surface",
					surfacecolor=t(pm.mx2.wROSg),
					colorscale = list(c(0, 1), c("#3366ff", 'white')),
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
					colorscale = list(c(0, 1), c("#567001", 'white')),
					opacity = 0.7,
					showscale = FALSE)
			}
			if(input$pm.mx1.wROSCheckbox == "solid") {
				pm.mx1.wROSg <- pm.mx1.wROS()
				pm.mx1.wROSg[!is.na(pm.mx1.wROSg)] <- 1
				p <- add_surface(p, x = ~vecx1(), y = ~vecx2(), z = ~t(pm.mx1.wROS()),
					type = "surface",
					surfacecolor=t(pm.mx1.wROSg),
					colorscale = list(c(0, 1), c("#567001", 'white')),
					opacity = 0.7,
					showscale = FALSE)
			}
		} else {
			NULL
		}
		
		# add plot features if checked --------------------------------------------
		# crossover line for x2 as mod
		if (!is.null(input$coonx1Checkbox)) {
			if(input$coonx1Checkbox == TRUE) {
				p <- add_trace(p, x = coonx1(), y = vecx2(), z = ~coonx1y(),
					type = "scatter3d",
					mode = "lines",
					line = list(width = 3, dash = "solid", color = "#567001"),
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
					line = list(width = 3, dash = "solid", color = "#3366ff"),
					opacity = 0.7)
			}
		} else {
			NULL
		}
		p
	})
	
	
	# scale for cb gradient ---------------------------------------------------
	# for x1 as mod	
	jnmx1cbd <- reactive({	
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
			}
		} else {
			NULL
		}
	})
	
	jnmx1cbd.print <- reactive({	
		jnmx1cbd <- seq(jnmx1cbd()[1], jnmx1cbd()[length(jnmx1cbd())], length.out = 10)
		jnmx1cbd <- sort(jnmx1cbd, decreasing = FALSE)
		jnmx1cbd <- format(round(jnmx1cbd, digits = 3), nsmall = 3)
		jnmx1cbd
	})
	
	output$mx1cbscale <- renderUI({
		req(non.numeric() == FALSE, sameVar() == FALSE)
		if(intsig() == TRUE & !is.null(input$pm.mx1.wROSCheckbox) & mx1.ROS.none() == FALSE) {
			if(input$pm.mx1.wROSCheckbox == "gradient") {
				HTML("<div style='background-image: linear-gradient(to right, #567001, white);  
							border-width: 1px; border-style: solid; border-color:gray;
							height: 20px;   text-align: justify' </div>")
			}
		} else {
			NULL
		}
	})
	
	
	output$mx1cbscale.title <- renderUI({
		req(non.numeric() == FALSE, sameVar() == FALSE)
		if(intsig() == TRUE & !is.null(input$pm.mx1.wROSCheckbox) & mx1.ROS.none() == FALSE) {
			if(input$pm.mx1.wROSCheckbox == "gradient") {
				HTML("<div style='height: 20px;   text-align: center' >", paste0("Width of 95% CB around ", input$x2, " Slope"), "</div>")
			}
		} else {
			NULL
		}
	})
	
	output$mx1cbscale.text1 <- renderUI({
		req(non.numeric() == FALSE, sameVar() == FALSE)
		if(intsig() == TRUE & !is.null(input$pm.mx1.wROSCheckbox) & mx1.ROS.none() == FALSE) {
			if(input$pm.mx1.wROSCheckbox == "gradient") {
				HTML("<div style='height: 20px;   text-align: left' >",  paste0("Min: ", min(jnmx1cbd.print())), "</div>")
			}
		} else {
			NULL
		}
	})
	
	output$mx1cbscale.text2 <- renderUI({
		req(non.numeric() == FALSE, sameVar() == FALSE)
		if(intsig() == TRUE & !is.null(input$pm.mx1.wROSCheckbox) & mx1.ROS.none() == FALSE) {
			if(input$pm.mx1.wROSCheckbox == "gradient") {
				HTML("<div style='height: 20px;   text-align: right' >",  paste0("Max: ", max(jnmx1cbd.print())), "</div>") 
			}
		} else {
			NULL
		}
	})
	
	# for x2 as mod
	jnmx2cbd <- reactive({	
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
			}
		} else {
			NULL
		}
	})
	
	jnmx2cbd.print <- reactive({	
		jnmx2cbd <- seq(jnmx2cbd()[1], jnmx2cbd()[length(jnmx2cbd())], length.out = 10)
		jnmx2cbd <- sort(jnmx2cbd, decreasing = FALSE)
		jnmx2cbd <- format(round(jnmx2cbd, digits = 3), nsmall = 3)
		jnmx2cbd
	})
	
	output$mx2cbscale <- renderUI({
		req(non.numeric() == FALSE, sameVar() == FALSE)
		if(intsig() == TRUE & !is.null(input$pm.mx2.wROSCheckbox) & mx2.ROS.none() == FALSE) {
			if(input$pm.mx2.wROSCheckbox == "gradient") {
				HTML("<div style='background-image: linear-gradient(to right, #3366ff, white);  
							border-width: 1px; border-style: solid; border-color:gray;
							height: 20px;   text-align: justify' </div>")
			}
		} else {
			NULL
		}
	})
	
	output$mx2cbscale.title <- renderUI({
		req(non.numeric() == FALSE, sameVar() == FALSE)
		if(intsig() == TRUE & !is.null(input$pm.mx2.wROSCheckbox) & mx2.ROS.none() == FALSE) {
			if(input$pm.mx2.wROSCheckbox == "gradient") {
				HTML("<div style='height: 20px;   text-align: center' >", paste0("Width of 95% CB around ", input$x1, " Slope"), "</div>")
			}
		} else {
			NULL
		}
	})
	
	output$mx2cbscale.text1 <- renderUI({
		req(non.numeric() == FALSE, sameVar() == FALSE)
		if(intsig() == TRUE & !is.null(input$pm.mx2.wROSCheckbox) & mx2.ROS.none() == FALSE) {
			if(input$pm.mx2.wROSCheckbox == "gradient") {
				HTML("<div style='height: 20px;   text-align: left' >",  paste0("Min: ", min(jnmx2cbd.print())), "</div>")
			}
		} else {
			NULL
		}
	})
	output$mx2cbscale.text2 <- renderUI({
		req(non.numeric() == FALSE, sameVar() == FALSE)
		if(intsig() == TRUE & !is.null(input$pm.mx2.wROSCheckbox) & mx2.ROS.none() == FALSE) {
			if(input$pm.mx2.wROSCheckbox == "gradient") {
				HTML("<div style='height: 20px;   text-align: right' >",  paste0("Max: ", max(jnmx2cbd.print())), "</div>") 
			}
		} else {
			NULL
		}
	})
	
	
	# render plot and UI for plot features -------------------------------------
	output$plot3dPlot <- renderPlotly({
		req(non.numeric() == FALSE, sameVar() == FALSE)
		plot3d()
	})
	
	output$plotfeatTitle <- renderText({
		req(input$x1 != '--', input$x2 != '--', input$y != '--', sameVar() == FALSE, non.numeric() == FALSE)
		HTML(paste0("<div style='font-weight: 700' >Plot Features</div>"))
	})
	
	output$scatterUI <- renderUI({
		req(input$x1 != '--', input$x2 != '--', input$y != '--', sameVar() == FALSE, non.numeric() == FALSE)
		checkboxInput(inputId = "scatterCheckbox",
			label = "Scatterplot",
			value = checkboxes$scatterCB)
	})
	output$regplaneUI <- renderUI({
		req(input$x1 != '--', input$x2 != '--', input$y != '--', sameVar() == FALSE, non.numeric() == FALSE)
		checkboxInput(inputId = "regplaneCheckbox",
			label = "Regression plane",
			value = checkboxes$regplaneCB)
	}) 
	output$predCIUI <- renderUI({
		req(input$x1 != '--', input$x2 != '--', input$y != '--', sameVar() == FALSE, non.numeric() == FALSE)
		checkboxInput(inputId = "predCICheckbox", 
			label = "95% CI",
			value = checkboxes$predCICB)
	})
	output$pm.mx1.wROSUI <- renderUI({
		if ((mx1.ROS.none() == FALSE) & intsig() == TRUE & non.numeric() == FALSE) {
			radioButtons(inputId = "pm.mx1.wROSCheckbox",
				label = NULL,
				choices = c("No ROS" = "none", "Solid ROS" = "solid", "Gradient ROS (Slope 95% CB)" = "gradient"),
				selected = checkboxes$pm.mx1.wROSCB, inline = FALSE)
		}
	})
	output$pm.mx2.wROSUI <- renderUI({
		if ((mx2.ROS.none() == FALSE) & intsig() == TRUE & non.numeric() == FALSE) {
			radioButtons(inputId = "pm.mx2.wROSCheckbox",
				label = NULL,
				choices = c("No ROS" = "none", "Solid ROS" = "solid", "Gradient ROS (Slope 95% CB)" = "gradient"),
				selected = checkboxes$pm.mx2.wROSCB, inline = FALSE)
		}
	})
	
	output$coonx1UI <- renderUI({
		if (!is.na(coonx1()) & non.numeric() == FALSE & sameVar() == FALSE) {
			checkboxInput(inputId = "coonx1Checkbox", 
				label = "Crossover",
				value = checkboxes$coonx1CB)
		} 
	})
	
	output$coonx2UI <- renderUI({
		if (!is.na(coonx2()) & non.numeric() == FALSE & sameVar() == FALSE) {
			checkboxInput(inputId = "coonx2Checkbox", 
				label = "Crossover",
				value = checkboxes$coonx2CB)
		} 
	})
	
	output$mx1CheckboxTitle <- renderText({
		req(input$x1 != '--', input$x2 != '--', input$y != '--', sameVar() == FALSE, non.numeric() == FALSE)
		if((intsig() == TRUE & mx1.ROS.none() == FALSE) | !is.na(coonx1())) {
			HTML(paste0("<div style='font-weight: 700' >", input$x1, " as Moderator </div>")) 
		} else {
			NULL
		}
	})
	
	output$mx2CheckboxTitle <- renderText({
		req(input$x1 != '--', input$x2 != '--', input$y != '--', sameVar() == FALSE, non.numeric() == FALSE)
		if((intsig() == TRUE & mx2.ROS.none() == FALSE) | !is.na(coonx2())) {
			HTML(paste0("<div style='font-weight: 700' >", input$x2, " as Moderator </div>")) 
		} else {
			NULL
		}
	})
	
	# render warnings/messages for sidebar ------------
	output$messageSameVar <- renderText({
		req(input$x1 != "--", input$x2 != "--", input$y != "--")
		if(sameVar() == TRUE) {
			HTML(paste0("<div style='background-color: #ff9999; border-radius: 5px;
				text-align: center; border-width: 1px; border-style: solid; border-color: #808080; padding:5px;
				font-size: 10px; font-color: #808080'> Check your model. The same variable was input twice. </div>"))
		} 
	})
	
	output$messageNonNum <- renderText({
		req(input$x1 != "--", input$x2 != "--", input$y != "--")
		if(non.numeric() == TRUE) {
			HTML(paste0("<div style='background-color: #ff9999; border-radius: 5px;
			text-align: center; border-width: 1px; border-style: solid; border-color: #808080; padding:5px;
			font-size: 12px; font-color: #808080'> Check your variables. At least one seems to be non-numeric. If
		  they are all numeric, it may be that one variable has only a few integer values or (if using an SPSS file) has value labels and is being read as a factor. Try 
			using a different file type and/or removing the value labels. </div>"))
		} 
	})
	
	output$messageInt <- renderText({
		req(input$x1 != "--", input$x2 != "--", input$y != "--", non.numeric() == FALSE, sameVar() == FALSE)
		if(intsig() == FALSE) {
			HTML(paste0("<div style='background-color: #ff9999; border-radius: 5px;
					text-align: center; border-width: 1px; border-style: solid; border-color:#808080;
					font-size: 10px; font-color: #808080; padding:10px;' >", 
				paste0("The interaction between ", input$x1, " and ", input$x2, " does not significantly predict ",
					input$y, ". To view all 3D plot features and the Johnson-Neyman test results, enter a model with a significant ",
					"interaction term. You can still see the regression plane and scatterplot on the <b> Main Output </b> tab."), "</div>"))
		} else {
			HTML(paste0("<div style='background-color: #d9d9d9; border-radius: 5px;
					text-align: center; border-width: 1px; border-style: solid; border-color: #808080; padding:8px;
					font-size: 10px; font-color: #808080'> Click the <b> Main Output </b> tab to see the 3D plot and output for your model. Scroll down to see all plot features. </div>"))
		}
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
	
}

shinyApp(ui, server)

