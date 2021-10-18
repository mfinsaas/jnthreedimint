# this code makes the output of johnson_neyman whether the ROS fall between jn values 
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

