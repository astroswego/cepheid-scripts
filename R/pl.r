#!/usr/bin/env Rscript

library(optparse)
library(scatterplot3d)

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

correct_extinction <- function(data, use_zeroes = TRUE, average = 0.0) {
    if (use_zeroes) {
        data$ext[data$ext == 0.0] <- average
    }
    else {
        data <- data[data$ext != 0.0,]
    }
    
    data$mag <- data$mag - data$ext

    data
}

formula_string <- function(ystr, ...) {
    # string representation of right hand side of equation
    rhs <- paste(..., sep = " + ")
    # string representation of entire equation
    paste(ystr, rhs, sep = " ~ ")
}

do_fit <- function(outdir, fname, data, params, ystr, ...) {
    # formula from string representation of equation
    eqn <- as.formula(formula_string(ystr, ...))

    model <- lm(eqn)
    print(summary(model))

    # 2D plot
    fpath <- file.path(outdir, fname)
    png(fpath)
    plot(data$logP, data$mag)
    points(data$logP, model$fitted.values, col="red")
    dev.off()
    if (length(params) == 2) {
        # 3D plot
        fpath3d <- file.path(outdir, paste("3D", fname, sep = "-"))
        png(fpath3d)
        s3d <- scatterplot3d(data$logP, params[[2]], data$mag)
        s3d$plane3d(model)
        dev.off()
    }
}

do_fits <- function(outdir, datastr, data, tag) {
    magstr <- paste(datastr, "mag", sep = "$")
    logPstr <- paste(datastr, "logP", sep = "$")
    PC1str <- paste(datastr, "PC1", sep = "$")
    PC2str <- paste(datastr, "PC2", sep = "$")
    
    do_fit(outdir,
           paste("PL-", tag, ".png", sep = ""),
           data, list(data$logP), magstr, logPstr)
    do_fit(outdir,
           paste("PLPC1-", tag, ".png", sep = ""),
           data, list(data$logP, data$PC1), magstr, logPstr, PC1str)
    do_fit(outdir,
           paste("PLPC2-", tag, ".png", sep = ""),
           data, list(data$logP, data$PC2), magstr, logPstr, PC2str)
    do_fit(outdir,
           paste("PLPC12-", tag, ".png", sep = ""),
           data, list(data$logP, data$PC1, data$PC2),
           magstr, logPstr, PC1str, PC2str)
}


opts <- get_options()
data <- read.table(opts$input, header = TRUE)
data$mag <- data$mag - opts$M_0
data <- correct_extinction(data,
                           use_zeroes = opts$zeroes,
                           average = opts$avg_ext)
split_data <- split_period(data, split = opts$period_split)
low_data <- split_data[[1]]
high_data <- split_data[[2]]

cat("All logP")
do_fits(opts$output, "data", data, "all")

cat(paste("logP < ", opts$period_split, "\n", sep = ""))
do_fits(opts$output, "low_data", low_data, "low")

cat(paste("logP >= ", opts$period_split), "\n", sep = "")
do_fits(opts$output, "high_data", high_data, "high")
