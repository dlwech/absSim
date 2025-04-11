####################################################################################################################################
# Compute absolute similarity (reversed absolute difference) score from any two continuous measures
####################################################################################################################################
getZ <- function(x) (x - mean(x, na.rm=T)) / sd(x, na.rm=T)

roundUp <- function(x, to=10) to * (x %/% to + as.logical(x %% to))

absSim <- function(df, varX, varY, varAux=NULL, plots=F, histograms=F, report=T, debug=F, completePairs=T) {
  if (!is.data.frame(df)) stop("Argument df must contain dataframe")
  if (!is.character(varX) | !is.character(varY)) stop("Arguments varX and varY must contain variable names as strings")
  if (!all(c(varX, varY) %in% names(df))) stop("Dataframe must contain varX and varY")
  if (all(is.na(df[, varX])) | all(is.na(df[, varY]))) stop("One or more variables has no data")
  if (!is.numeric(df[, varX]) | !is.numeric(df[, varY])) stop("Variables varX and varY must be numeric")
  if (sum(!is.na(unique(df[, varX]))) < 3 | sum(!is.na(unique(df[, varY]))) < 3) stop("Variables varX and varY must not be binary")
  if (!is.null(varAux)) {
    if (any(!is.character(varAux))) stop("Argument varAux must contain auxiliary variable names as strings")
    if (!all(varAux %in% names(df))) stop("Dataframe must contain all auxiliary variables")}
  
  bar <- "--------------------------------------------------------------------------------------\n"
  if (debug) {
    cat(paste0(bar, "Calculating similarity scores with debugging output. See readme for options\n"))
    report <- T; plots <- T; histograms <- T
    }
  else cat(paste0(bar, "Calculating similarity scores. See readme for options and debugging\n"))
  
  if (is.null(varAux)) {
    cat(bar); warning("
  If any portion of the sample is missing data on intended outcome/other 
  analytic variables, entering these as auxiliary variables is advisable to 
  ensure that standardised scores used to calculate absolute similarity are 
  based only on cases to be included in analyses. This should help avoid 
  ceiling/floor effects and reduce heteroscedasticity",
  "\nOne or more auxiliary variable names can be entered in the varAux argument\n")
    }
  if (any(duplicated(df[!is.na(df[, 1]), 1]))) {
    cat(bar); warning("
  Values in first column are expected to be unique IDs, but contain duplicates. This 
  function assumes a wide format dataset, with each row containing scores for unique 
  pairs of individuals. Similarity scores may be erroneous if data are structured in 
  long format or have duplicated rows\n")
    }
  
  # Warn if data contain integers with narrow range, or floats with few unique values
  intVars <- c()
  warningString <- "\n"
  for (i in c(varX, varY)) {
    v <- df[!is.na(df[, i]), i]
    if (all(v == as.integer(v))) {
      r <- range(df[, i], finite=T)
      ar <- abs(r[1] - r[2])
      if (ar < 10) warningString <- paste0(warningString, "  ", i, " is discrete and has a range < 10\n")
      }
    else if (length(unique(v))) warningString <- paste0(warningString, "  ", i, " is continuous but has < 10 unique values\n")
    }
  if (warningString != "\n") {
    warningString <- paste0(warningString, 
  "This function is designed for continuous questionnaire data with a reasonable level of 
  variablility. It should not be used for categorical or ordinal data. If scores are 
  continuous measures with low variability, exercise caution and consult diagnostics 
  during analysis to check for issues with outliers, heteroscedasticity, and 
  ceiling/floor effects\n")
    cat(bar); warning(warningString)
    }
  
  allXY <- sum(rowSums(is.na(df[, c(varX, varY)])) == 0)
  anyXY <- sum(rowSums(is.na(df[, c(varX, varY)])) < 2)
  cat(paste0(bar, "Dataframe input with ", nrow(df)," rows\n  ", anyXY, " cases with any varX or varY data; ", allXY, " complete pairs\n"))
  
  # Remove cases missing on auxiliary variables
  if (!is.null(varAux)) {
    df <- df[, c(varX, varY, varAux)]
    cat(paste0(bar, "Checking auxiliary data for ", paste(varAux, collapse=", "), "\n"))
    exc <- 0
    for (a in varAux) {
      if (all(is.na(df[, a]))) stop(paste0("Auxiliary variable ", a, " has no data"))
      missingAux <- sum(is.na(df[, a]) & rowSums(!is.na(df[, c(varX, varY)])) > 0)
      if (missingAux == 0) cat(paste0("  All cases have ", a, " data\n"))
      else {
        df[is.na(df[, a]), c(varX, varY)] <- NA
        exc <- exc + missingAux
        cat(paste0("  ", missingAux, " cases excluded due to missing ", a, " data\n"))
      }
    }
    if (exc == 0) cat(paste0("No cases excluded\n"))
    else {
      cat(paste0(exc, " total cases excluded\n"))
      allXY <- sum(rowSums(is.na(df[, c(varX, varY)])) == 0)
      anyXY <- sum(rowSums(is.na(df[, c(varX, varY)])) < 2)
      cat(paste0("  ", anyXY, " cases remaining; ", allXY, " complete pairs\n"))
      }
    }
  
  # Remove unpaired scores (or warn if completePairs == FALSE)
  df <- df[, c(varX, varY)]
  if (allXY != anyXY) {
    if (completePairs) {
      df[rowSums(!is.na(df[, c(varX, varY)])) != 2, ] <- NA
      cat(paste0(bar, "Raw scores standardised against data from ", allXY," complete pairs\n"))
      }
    else {
      onlyX <- sum(!is.na(df[, varX]) & is.na(df[, varY]))
      onlyY <- sum(!is.na(df[, varY]) & is.na(df[, varX]))
      warningString <- "\n"
      if (onlyY > 0) warningString <- paste0(warningString, "  ", onlyY, " cases missing ", varX, " data\n")
      if (onlyX > 0) warningString <- paste0(warningString, "  ", onlyX, " cases missing ", varY, " data\n")
      cat(bar); warning(paste0(warningString, "completePairs set to FALSE
  Raw scores will be standardised against all available data for each variable,
  including ", onlyX + onlyY, " unpaired scores. This method should be used with careful consideration 
  and assessment of diagnostics during analysis. Depending on sample characteristics, 
  it could contribute to ceiling/floor effects, outliers, and heteroscedasticity\n"))
      }
    }
  
  # Calculate absolute similarity scores
  x <- df[, varX]
  y <- df[, varY]
  Zx <- getZ(x)
  Zy <- getZ(y)
  DiffZ <- Zx - Zy
  AbsDiffZ <- abs(DiffZ)
  AbsSimZ <- max(AbsDiffZ, na.rm=T) - AbsDiffZ
  
  # For debugging - unstandardised equivalents in case Z-scoring is causing issues
  if (debug) {
    Diff <- x - y
    AbsDiff <- abs(Diff)
    AbsSim <- max(AbsDiff, na.rm=T) - AbsDiff
    }
  
  if (histograms) {
    if (!nzchar(system.file(package = "ragg"))) cat(paste0(bar, "For histograms, run again with ragg package installed\n"))
    else {
      cat(paste0(bar, "Histograms\n", bar))
      graphics::par(mfrow = c(1, 1))
      if (debug) graphics::hist(AbsSim, col=rgb(0,0,0,1/4), breaks=20, xlim=c(0, roundUp(max(AbsSim, na.rm=T), 10)))
      graphics::hist(AbsSimZ, col=rgb(0,0,0,1/4), breaks=20, xlim=c(0, roundUp(max(AbsSimZ, na.rm=T), 1)))
      }
    }
  
  # Plots, reports, debugging outputs if needed for data checks
  if (plots) {
    if (!nzchar(system.file(package = "ggplot2")) | !nzchar(system.file(package = "ragg"))) cat(paste0(bar, "For plots, run again with ggplot2 and ragg packages installed\n"))
    else {
      cat(paste0(bar, "Plots\n", bar))
      ZxG <- Zx[!is.na(AbsSimZ)]
      ZyG <- Zy[!is.na(AbsSimZ)]
      AbsSimZG <- AbsSimZ[!is.na(AbsSimZ)]
      wrapper <- function(x, ...) paste(strwrap(x, ...), collapse = "\n")
      st <- "Lines on the y-axis represent absolute differences (i.e., distances) between standardised scores on varX (blue dots) and varY (red dots). Each pair's position on the x-axis represents its absolute similarity (i.e., reversed absolute difference) score. More similar pairs therefore aggregate on the right, regardless of who scores higher, or whether they score overall higher/lower than other pairs."
      diffplot <- ggplot2::ggplot() + 
        ggplot2::geom_linerange(ggplot2::aes(x=AbsSimZG, ymin=ZyG, ymax=ZxG)) + 
        ggplot2::ylab("Std. scores") + ggplot2::xlab("AbsSimZ") +
        ggplot2::geom_point(ggplot2::aes(x=AbsSimZG, y=ZxG), colour="#0077FF", size=0.9) +
        ggplot2::geom_point(ggplot2::aes(x=AbsSimZG, y=ZyG), colour="#DD1111", size=0.9) +
        ggplot2::labs(title="Distance/similarity plot", subtitle=wrapper(st, width=66)) + 
        ggplot2::theme(plot.title=ggplot2::element_text(face="bold", size=12, hjust=0.5)) + 
        ggplot2::theme(plot.subtitle=ggplot2::element_text(size=11))
      print(diffplot)
      
      if (debug) {
        st <- "Here, each pair's position on the x-axis represents their score on a conventional interaction term. Any pairs scoring closer to the mean aggregate in the centre (i.e., 'average' similarity) even if their standardised scores are nearly identical. Only those with similar deviations from the mean aggregate on the right, and only those with opposing deviations from the mean aggregate on the left."
        intplot <- ggplot2::ggplot() + 
          ggplot2::geom_linerange(ggplot2::aes(x=ZxG*ZyG, ymin=ZyG, ymax=ZxG)) + 
          ggplot2::ylab("Std. scores") + ggplot2::xlab("Zx * Zy") +
          ggplot2::geom_point(ggplot2::aes(x=ZxG*ZyG, y=ZxG), colour="#0077FF", size=0.9) +
          ggplot2::geom_point(ggplot2::aes(x=ZxG*ZyG, y=ZyG), colour="#DD1111", size=0.9) +
        ggplot2::labs(title="Distance/interaction term plot", subtitle=wrapper(st, width=66)) + 
        ggplot2::theme(plot.title=ggplot2::element_text(face="bold", size=12, hjust=0.5)) +
        ggplot2::theme(plot.subtitle=ggplot2::element_text(size=11))
        print(intplot)
        }
      }
    }
  
  if (debug) {
    cat(paste0(bar, "[Debug] Full calculation steps (including unstandardised alternatives):\n"))
    cat(paste0(bar, varX, " raw scores (x):\n")); print(x)
    cat(paste0(bar, varY, " raw scores (y):\n")); print(y)
    cat(paste0(bar, "Difference (Diff):\n")); print(Diff)
    cat(paste0(bar, "Abs. difference (AbsDiff):\n")); print(AbsDiff)
    cat(paste0(bar, "Abs. similarity (AbsSim):\n")); print(AbsSim)
    cat(paste0(bar, varX, " Z-scores (Zx):\n")); print(round(Zx, 2))
    cat(paste0(bar, varY, " Z-scores (Zy):\n")); print(round(Zy, 2))
    cat(paste0(bar, "Z-score difference (DiffZ):\n")); print(round(DiffZ, 2))
    cat(paste0(bar, "Abs. Z-score difference (AbsDiffZ):\n")); print(round(AbsDiffZ, 2))
    cat(paste0(bar, "Abs. Z-score similarity (AbsSimZ):\n")); print(round(AbsSimZ, 2))
    }
  
  # Report of descriptives of absolute similarity measure, and correlations with other measures
  if (report) {
    cors <- list()
    vars <- list(standardised = AbsSimZ)
    if (debug) vars <- c(vars, list(unstandardised = AbsSim))
    for (var in seq_along(vars)) {
      v <- vars[[var]]
      if (var == 2 & debug) Ss <- "[Debug]\n" else Ss <- ""
      Ds <- paste0(bar, Ss, "Summary of ", names(vars)[var], " absolute similarity scores (", varX," - ", varY, ")",
        "\n  N = ", sum(!is.na(v)), "; Range: ", round(min(v[!is.na(v)]), 2), " - ", round(max(v[!is.na(v)]), 2), "; ",
        "Mean (SD) = ", round(mean(v[!is.na(v)]), 2), " (", round(sd(v[!is.na(v)]), 2), ")")
      if (nzchar(system.file(package = "moments"))) {
        Ds <- paste0(Ds, "\n  Skewness, Kurtosis = ", round(moments::skewness(v[!is.na(v)]), 2), 
          ", ", round(moments::kurtosis(v[!is.na(v)]), 2))}
      else Ds <- paste0(Ds, "\n  For normality statistics, run again with moments package installed")
      cat(Ds)
      Rx <- cor.test(v, x, method="pearson")
      Ry <- cor.test(v, y, method="pearson")
      RxE <- round(as.numeric(Rx$estimate), 2); RxP <- signif(Rx$p.value, 3)
      RyE <- round(as.numeric(Ry$estimate), 2); RyP <- signif(Ry$p.value, 3)
      cat(paste0("\n    Pearson's R with ", varX, " = ", RxE, " (p = ", RxP, ")"))
      cat(paste0("\n    Pearson's R with ", varY, " = ", RyE, " (p = ", RyP, ")"))
      cors <- c(cors, list(list(RxE = RxE, RyE = RyE)))
      cat("\n")
      }
    
    if (debug) {
      R <- cor.test(vars[['unstandardised']], vars[['standardised']], method="pearson")
      RE <- round(as.numeric(R$estimate), 3)
      RP <- signif(as.numeric(R$p.value), 3)
      cat(paste0(bar, "[Debug]\nPearson's R between AbsSim and AbsSimZ: ", RE, " (p = ", RP, ")\n"))
      
      # Return AbsSimZ variable only, unless debugging, in which case return all calculation steps as one variable
      cat(paste0(bar, "[Debug]\nReturning all calculation steps for data checking purposes. 
  NOTE: Using any intermediate score may lead to erroneous inference, particularly 
  when scores are calculated from undstandardised variables not on the same scale\n", bar))
      return(cbind(df, cbind(Diff, AbsDiff, AbsSim, Zx, Zy, DiffZ, AbsDiffZ, AbsSimZ)))
      }
    }
  cat(bar)
  AbsSimZ
  }
