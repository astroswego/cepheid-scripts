#!/usr/bin/env Rscript

library(optparse)

get_options <- function() {
    option_list <- list(
        make_option(c("-i", "--input"), metavar = "INFILE",
                    action = "store", type = "character",
                    help = paste(
                        "Input table with header containing",
                        "mag, ext, logP, PC1, and PC2."
                        )
                    ),
        make_option(c("-o", "--output"), metavar = "OUTDIR",
                    action = "store", type = "character",
                    help = "Output directory."),
        make_option("--period-split", metavar = "P", dest = "period_split",
                    action = "store", type = "double", default = 1.0,
                    help = "Period to split dataset at."),
        make_option("--dist-mod", metavar = "M_0", dest = "M_0",
                    action = "store", type = "double", default = 0.0,
                    help = paste(
                        "Distance modulus of target, to convert magnitudes",
                        "from apparent to absolute."
                        )
                    ),
        make_option("--conversion", metavar = "CONV", dest = "conversion",
                    action = "stor", type = "double", default = 3.24,
                    help = "Conversion from E(V-I)/1.4 to A_V"),
        make_option("--avg-ext", metavar = "AVG", dest = "avg_ext",
                    action = "store", type = "double", default = 0.0,
                    help = "Average extinction for dataset."),
        make_option("--use-zeroes", dest = "zeroes",
                    action = "store_true", default = FALSE,
                    help = "Use stars with no known extinction values.")
        )
    parse_args(OptionParser(prog = "pl",
                            option_list = option_list))
}

split_period <- function(data, split = 1.0) {
    s <- data$logP < split
    low_data <- data[s,]
    high_data <- data[!s,]

    list(low_data, high_data)
}

correct_extinction <- function(data, use_zeroes = TRUE, average = 0.0,
                               conversion = 3.24) {
    if (use_zeroes)
        data$A_V[data$A_V == 0.0] <- average
    else
        data <- data[data$A_V != 0.0,]

    IV_ratio <- 0.604937
    data$A0_I <- data$A0_I - IV_ratio*data$A_V
    data$A0_V <- data$A0_V - data$A_V
    data$color <- data$color - (1.4/conversion)*data$A_V

    data
}

formula_string <- function(ystr, ...) {
    # string representation of right hand side of equation
    rhs <- paste(..., sep = " + ")
    # string representation of entire equation
    paste(ystr, rhs, sep = " ~ ")
}

do_fit <- function(outdir, fname, formula, data, labels) {
    fit <- lm(formula)

    # instead of printing the summary using its print method, extract the
    # data from the summary and output it in a table
    print(summary(fit))

    # 2D plot
    fpath <- file.path(outdir, fname)
    png(fpath)
    
    plot(fit$model[[2]], fit$model[[1]],
         xlab = labels[[2]], ylab = labels[[1]])
    points(fit$model[[2]], fit$fitted.values, col="red")
    dev.off()
}

name_string <- function(prefix, tag) {
    paste(prefix, "-", tag, ".png", sep = "")
}

do_fits <- function(outdir, data, tag) {
    do_fit(outdir, name_string("PL", tag),
           data$A0 ~ data$logP, data,
           list("A0", "logP"))
    do_fit(outdir, name_string("PLC", tag),
           data$A0 ~ data$logP + data$color, data,
           list("A0", "logP", "V-I"))
    do_fit(outdir, name_string("PLR21", tag),
           data$A0 ~ data$logP + data$R21, data,
           list("A0", "logP", "R21"))
    do_fit(outdir, name_string("PLR31", tag),
           data$A0 ~ data$logP + data$R31, data,
           list("A0", "logP", "R31"))
    do_fit(outdir, name_string("PLPC1", tag),
           data$A0 ~ data$logP + data$PC1, data,
           list("A0", "logP", "PC1"))
    do_fit(outdir, name_string("PLPC2", tag),
           data$A0 ~ data$logP + data$PC2, data,
           list("A0", "logP", "PC2"))
    do_fit(outdir, name_string("PLPC3", tag),
           data$A0 ~ data$logP + data$PC3, data,
           list("A0", "logP", "PC3"))
}

opts <- get_options()
data <- read.table(opts$input, header = TRUE)

data$color <- data$A0_V - data$A0_I
data$A0_I <- data$A0_I - opts$M_0
data$A0_V <- data$A0_V - opts$M_0
data <- correct_extinction(data,
                           use_zeroes = opts$zeroes,
                           average = opts$avg_ext,
                           conversion = opts$conversion)
#split_data <- split_period(data, split = opts$period_split)
#low_data <- split_data[[1]]
#high_data <- split_data[[2]]

data_I <- data.frame(
    logP  = data$logP,
    A0    = data$A0_I,
    color = data$color,
    R21   = data$R21_I,
    R31   = data$R31_I,
    PC1   = data$PC1_I,
    PC2   = data$PC2_I,
    PC3   = data$PC3_I
)
data_V <- data.frame(
    logP  = data$logP,
    A0    = data$A0_V,
    color = data$color,
    R21   = data$R21_V,
    R31   = data$R31_V,
    PC1   = data$PC1_V,
    PC2   = data$PC2_V,
    PC3   = data$PC3_V
)
## low_data_I <- data.frame(
##     logP  = low_data$logP,
##     A0    = low_data$A0_I,
##     color = low_data$color,
##     R21   = low_data$R21_I,
##     R31   = low_data$R31_I,
##     PC1   = low_data$PC1_I,
##     PC2   = low_data$PC2_I,
##     PC3   = low_data$PC3_I
## )
## low_data_V <- data.frame(
##     logP  = low_data$logP,
##     A0    = low_data$A0_V,
##     color = low_data$color,
##     R21   = low_data$R21_V,
##     R31   = low_data$R31_V,
##     PC1   = low_data$PC1_V,
##     PC2   = low_data$PC2_V,
##     PC3   = low_data$PC3_V
## )
## high_data_I <- data.frame(
##     logP  = high_data$logP,
##     A0    = high_data$A0_I,
##     color = high_data$color,
##     R21   = high_data$R21_I,
##     R31   = high_data$R31_I,
##     PC1   = high_data$PC1_I,
##     PC2   = high_data$PC2_I,
##     PC3   = high_data$PC3_I
## )
## high_data_V <- data.frame(
##     logP  = high_data$logP,
##     A0    = high_data$A0_V,
##     color = high_data$color,
##     R21   = high_data$R21_V,
##     R31   = high_data$R31_V,
##     PC1   = high_data$PC1_V,
##     PC2   = high_data$PC2_V,
##     PC3   = high_data$PC3_V
## )
    
do_fits(opts$output, data_I, "I")
do_fits(opts$output, data_V, "V")
