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
                    help = "Use stars with no known extinction values."),
        make_option("--find-distance-modulus", dest = "find_distance_modulus",
                    action = "store_true", default = FALSE,
                    help = "Output distance modulus for each star.")
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
#    print(summary(fit))
    print_model(fit, fname)

    # 2D plot
    fpath <- file.path(outdir, paste(fname, ".png", sep = ""))
    png(fpath)
    
    plot(fit$model[[2]], fit$model[[1]],
         xlab = labels[[2]], ylab = labels[[1]])
    points(fit$model[[2]], fit$fitted.values, col="red")
    graphics.off()
}

print_model <- function(fit, fname) {
    sumry <- summary(fit)
    ncoefs <- nrow(sumry$coefficients)

    cat(fname)
    cat("\t")

    cat(sumry$coefficients[1,1])
    cat("\t")
    cat(sumry$coefficients[1,2])
    cat("\t")
    cat(sumry$coefficients[1,3])
    cat("\t")
    cat(sumry$coefficients[1,4])
    cat("\t")

    cat(sumry$coefficients[2,1])
    cat("\t")
    cat(sumry$coefficients[2,2])
    cat("\t")
    cat(sumry$coefficients[2,3])
    cat("\t")
    cat(sumry$coefficients[2,4])
    cat("\t")

    if (ncoefs > 2) {
        cat(sumry$coefficients[3,1])
        cat("\t")
        cat(sumry$coefficients[3,2])
        cat("\t")
        cat(sumry$coefficients[3,3])
        cat("\t")
        cat(sumry$coefficients[3,4])
        cat("\t")
    }
    else
        cat("NA\tNA\tNA\tNA\t")
    if (ncoefs > 3) {
        cat(sumry$coefficients[4,1])
        cat("\t")
        cat(sumry$coefficients[4,2])
        cat("\t")
        cat(sumry$coefficients[4,3])
        cat("\t")
        cat(sumry$coefficients[4,4])
        cat("\t")
    }
    else
        cat("NA\tNA\tNA\tNA\t")

    cat(sumry$sigma)
    cat("\t")

    cat(sumry$fstatistic[1])
    cat("\t")

    cat(sumry$r.squared)
    cat("\n")
}

name_string <- function(prefix, tag) {
    paste(prefix, "-", tag, sep = "")
}

do_fits <- function(outdir, data, tag) {
    # print heading
    cat("model\t")
    cat("a\ta_err\ta_t\ta_Pr\t")
    cat("b\tb_err\tb_t\tb_Pr\t")
    cat("c\tc_err\tc_t\tc_Pr\t")
    cat("d\td_err\td_t\td_Pr\t")
    cat("sigma\tF\tR2\n")

    do_fit(outdir, name_string("PL", tag),
           data$M0 ~ data$logP, data,
           list("M", "logP"))
    do_fit(outdir, name_string("PLC", tag),
           data$M0 ~ data$logP + data$color, data,
           list("M", "logP", "V-I"))
    do_fit(outdir, name_string("PLR21", tag),
           data$M0 ~ data$logP + data$R21, data,
           list("M", "logP", "R21"))
    do_fit(outdir, name_string("PLR31", tag),
           data$M0 ~ data$logP + data$R31, data,
           list("M", "logP", "R31"))
    do_fit(outdir, name_string("PLPC1", tag),
           data$M0 ~ data$logP + data$PC1, data,
           list("M", "logP", "PC1"))
    do_fit(outdir, name_string("PLPC2", tag),
           data$M0 ~ data$logP + data$PC2, data,
           list("M", "logP", "PC2"))
    ## do_fit(outdir, name_string("PLPC3", tag),
    ##        data$M0 ~ data$logP + data$PC3, data,
    ##        list("M", "logP", "PC3"))
    ## do_fit(outdir, name_string("PLA1", tag),
    ##        data$M0 ~ data$logP + data$A1, data,
    ##        list("M", "logP", "A1"))
    ## do_fit(outdir, name_string("PLA2", tag),
    ##        data$M0 ~ data$logP + data$A2, data,
    ##        list("M", "logP", "A2"))
    ## do_fit(outdir, name_string("PLA3", tag),
    ##        data$M0 ~ data$logP + data$A3, data,
    ##        list("M", "logP", "A3"))
    ## do_fit(outdir, name_string("PLA4", tag),
    ##        data$M0 ~ data$logP + data$A4, data,
    ##        list("M", "logP", "A4"))
    ## do_fit(outdir, name_string("PLA13", tag),
    ##        data$M0 ~ data$logP + data$A1 + data$A3, data,
    ##        list("M", "logP", "A1", "A3"))
}

distance_modulus <- function(data) {
    fit <- lm(data$M0 ~ data$logP + data$color)
    data$A0 - fit$fitted.values
}

distance_kPc <- function(distance_modulus) {
    Pc <- 10 * 10^(distance_modulus/5)
    Pc / 1000
}

opts <- get_options()
data <- read.table(opts$input, header = TRUE)

data$color <- data$A0_V - data$A0_I
data <- correct_extinction(data,
                           use_zeroes = opts$zeroes,
                           average = opts$avg_ext,
                           conversion = opts$conversion)
data$M_I <- data$A0_I - opts$M_0
data$M_V <- data$A0_V - opts$M_0



if (opts$find_distance_modulus) {
    data_I <- data.frame(
        logP  = data$logP,
        A0    = data$A0_I,
        M0    = data$M_I,
        color = data$color
    )
    data_V <- data.frame(
        logP  = data$logP,
        A0    = data$A0_V,
        M0    = data$M_V,
        color = data$color
    )
    
    I_modulus <- distance_modulus(data_I)
    V_modulus <- distance_modulus(data_V)

    I_kPc <- distance_kPc(I_modulus)
    V_kPc <- distance_kPc(V_modulus)

    write.table(
        data.frame(
            ID  = data$ID,
            d_V = V_kPc,
            d_I = I_kPc
        ),
        quote = FALSE,
        row.names = FALSE
    )
    ## print(summary(lm(V_kPc ~ I_kPc)))
    ## fpath <- file.path(opts$output, "I_vs_V.png")
    ## png(fpath)
    ## plot(I_kPc, V_kPc,
    ##      xlab = "I", ylab = "V")
    ## graphics.off()
} else {
    data_I <- data.frame(
        logP  = data$logP,
        A0    = data$A0_I,
        M0    = data$M_I,
#        A1    = data$A1_I,
#        A2    = data$A2_I,
#        A3    = data$A3_I,
#        A4    = data$A4_I,
        color = data$color,
        R21   = data$R21_I,
        R31   = data$R31_I,
        PC1   = data$PC1_I,
        PC2   = data$PC2_I
#        PC3   = data$PC3_I
    )
    data_V <- data.frame(
        logP  = data$logP,
        A0    = data$A0_V,
        M0    = data$M_V,
#        A1    = data$A1_V,
#        A2    = data$A2_V,
#        A3    = data$A3_V,
#        A4    = data$A4_V,
        color = data$color,
        R21   = data$R21_V,
        R31   = data$R31_V,
        PC1   = data$PC1_V,
        PC2   = data$PC2_V
#        PC3   = data$PC3_V
    )
    do_fits(opts$output, data_I, "I")
    do_fits(opts$output, data_V, "V")
}
