#!/usr/bin/env Rscript
path <- commandArgs(trailingOnly=TRUE)
coverage <- scan(path, double(), quote = "")
print(median(coverage, na.rm=TRUE))