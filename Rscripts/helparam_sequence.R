#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library(gsubfn)
require(stringr)
require(Hmisc)
if (length(args)==0) {
        stop("Usage: ./helparam_sequence.R <helparam.ser> <sequence>", call.=FALSE)
}
srcFile <- args[1]
seq <- args[2]

bpstep_param <- c("shift","slide","rise","tilt","roll","twist")
bp_param <- c("shear","stagger","stretch","buckle","propeller","opening")
param <- strapplyc(srcFile, "_(\\w+)")[[1]]
if (param %in% bpstep_param) {
        steps <- c(substring(seq, seq(1, nchar(seq)-1, 1), seq(2, nchar(seq), 1)),"")
} else {
        steps <- c(substring(seq, seq(1, nchar(seq), 1), seq(1, nchar(seq), 1)),"")
}
data <- read.table(srcFile,header=F)
data <- as.data.frame(data)
cdata <- data.frame(x <- seq(1:(ncol(data)-1)),
                    y <- sapply(data[,2:ncol(data)],mean),
                    sd <- sapply(data[,2:ncol(data)],sd))
png(file=paste(param,".png",sep = "", collapse = NULL),width=1000,height=500,res=72)
par(cex=2)
par(mar=c(5,5,2,3))
plt <- function(x,y,xmax,ymin,ymax,clr,lty) {
        plot(x,y,
             type="l",
             xlim=c(0,xmax),ylim=c(ymin,ymax),
             xlab="sequence",
             ylab=param,
             col=clr,axes=F,
             xaxs='i',yaxs='i',
             lwd=3,cex.lab=2,lty=lty,
             pch=20
        )
        with (
                data = cdata
                , expr = errbar(x, y, y+sd, y-sd, add=T, pch=1, cap=.02)
        )
}
ymin <- min(data[,2:ncol(data)])
ymax <- max(data[,2:ncol(data)])
plt(cdata$x,cdata$y,length(steps),ymin,ymax,1,1)
axis(2,cex.axis=1)
if (param %in% bpstep_param) {
        axis(1, 1:nchar(seq),
             labels=steps,
             cex.axis=1)
} else {
        axis(1, 1:(nchar(seq)+1),
             labels=steps,
             cex.axis=1)
}       
box(lwd=4)
dev.off()