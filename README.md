# absSim
Function for computing an absolute similarity (i.e., reversed absolute difference) score from any two continuous measures within an R dataframe. 

Tested on R versions 4.3.2 and 4.4.3 from within RStudio 2024.09.01 (Build 394). Optional plots and histograms implemented for RMarkdown with block settings `fig.width=5, fig.height=5` and default HTML output. 

## Usage
**Warning:** This tool is not guaranteed to be appropriate for all types of data and research question. Before drawing any inferences from findings, it is crucial to thoroughly consider what any similarity score may represent in your specific sample, and to rigorously assess diagnostics of any subsequent analyses, ensuring all statistical assumptions are met. It is always advisable to consider whether there is a potentially more effective or appropriate way to answer a given research question. 

---

The function assumes a wide format dataset wherein each row contains both individuals' scores (and any auxiliary variable data, i.e., outcome measures, covariates, and any other analytic variables; see `varAux` argument). A warning is given if the first column of the dataframe contains duplicated values, as these should be unique IDs. 

The first argument must be a dataframe containing all variables, and the second and third arguments must be the names of the two variables from which to compute similarity scores.

Optional arguments are available for providing auxiliary variables (recommended), producing plots and histograms, and more verbose output including all calculation steps (for debugging/diagnostics).
```
source('absSim.R')
myData$absSimilarity = absSim(df=myData, varX="var1", varY="var2")
```
## Optional arguments

### varAux
If any portion of the sample is missing data on outcome or other variables to be included in subsequent analyses (i.e., there are more data on varX and/or varY than there are on outcomes, covariates, or other analytic variables), entering these as auxiliary variables may be advisable. This will enforce listwise deletion of partially missing cases before calculating Z-scores which form the basis for absolute similarity scores. In many cases (though not always), this will help avoid ceiling/floor effects and reduce heteroscedasticity. Thoroughly assessing diagnostics of any subsequent analyses is necessary to ensure assumptions are met before drawing any inference from findings.

The name of one or more auxiliary variables can be passed as a string or a vector of strings, as in the below examples. This will enforce listwise deletion of partially missing cases, ensureing that Z-scores rely only on the group means/SDs of individuals who are not missing auxiliary data. 

```
varAux = "outcome"
```
```
varAux = c("outcome", "covariate1", "covariate2")
```

### plots
* Requires `ggplot2` and `ragg` packages. 

If set to `TRUE`, produces an annotated distance/similarity plot showing the standardised scores of each pair, with distances (i.e., absolute differences) represented by connecting lines on the y-axis. These distances are inversely proportional to absolute similarity scores of each pair, which is represented by their position on the x-axis. 

If `debug=TRUE`, this setting is forced to `TRUE`, and additionally produces an equivalent distance/interaction term plot where each pair's position on the x-axis represents the value of a conventional interaction term computed from their scores. 

### histograms
* Requires `ragg` package. 

If set to `TRUE`, produces a histogram of the computed absolute similarity score.

If `debug=TRUE`, this setting is forced to `TRUE`, and additionally produces a histogram of an absolute similarity score based on unstandardised scores. 

### report (default is TRUE)
* Normality statistics require `moments` package.

If set to `FALSE`, the function will not print a descriptive and correlational summary of the computed absolute similarity score.

If `debug=TRUE`, this setting is forced to `TRUE`, and several additions are introduced. A descriptive and correlational summary of an unstandardised similarity score equivalent is added, and a correlation is estimated between the unstandardised and standardised scores. 

### debug
Verbose output for diagnostic purposes only, and/or where data issues are suspected.

If set to `TRUE`:
* Coerces `plots`, `histograms`, and `report` to `TRUE` and extends their functionality (see respective descriptions for summaries)
* Returns a dataframe containing all intermediate variables (i.e., standardised scores on varX and varY, differences, absolute differences, and absolute similarities) along with unstandardised equivalents
  * Also prints these intermediate variables during computation 

### completePairs (default is TRUE)
Enforces listwise deletion of partially missing cases before calculating Z-scores which form the basis for absolute similarity scores. This ensures they rely only on the group means/SDs of individuals who are in a pair. 

In some cases, it may be appropriate to base Z-scores on a wider sample of individuals who are not in pairs. This should be done with careful consideration and thorough assessment of diagnostics during any subsequent analysis, as it could cause ceiling/floor effects and contribute to heteroscedasticity. 

If set to `FALSE`, Z-scores will be based on all available scores for varX and varY (after any exclusions due to missing auxiliary data), even those not in a pair. 
