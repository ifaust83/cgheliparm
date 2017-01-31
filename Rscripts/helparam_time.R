#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library(gsubfn)
if (length(args)==0) {
        stop("Usage: ./helparam_time.R <helparam.ser> <base pair (step)>", call.=FALSE)
}
srcFile <- args[1]
i <- as.integer(args[2])

param <- strapplyc(srcFile, "_(\\w+)")[[1]]
data <- read.table(srcFile,header=F)
png(file=paste(param,"_",i,".png",sep = "", collapse = NULL),width=1000,height=500,res=72)
par(cex=2)
par(mar=c(5,5,2,3))
plt <- function(x,y,xmax,ymin,ymax,clr,lty) {
        plot(x,y,
             type="l",
             xlim=c(0,xmax),ylim=c(ymin,ymax),
             xlab="time",
             ylab=param,
             col=clr,axes=F,
             xaxs='i',yaxs='i',
             lwd=1,cex.lab=2,lty=lty
        )
}
xmax <- max(data[,1])
ymin <- min(data[,i+1])
ymax <- max(data[,i+1])
plt(data[,1],data[,i+1],xmax,ymin,ymax,1,2)
axis(2,cex.axis=1.5)
axis(1,tck=0.05,lwd=4,cex.axis=1.5)
box(lwd=4)
dev.off()
