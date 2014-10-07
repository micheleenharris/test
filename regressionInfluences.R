# Load libraries

library(car)
library(lattice) #in case I need any trellis functions like xyplot and bwplot
library(gplots) # for smartlegend

# Read in files with user defined parameters below:
working_dir = "V:/mharris/tissueSpecific/compareInteract/pro/"
#working_dir <- "V:/mharris/mouseEcoli/"
#working_dir <- "V:/mharris/mouseSalm/"
#working_dir <- "V:/mharris/fluNS1/"

#working_dir = "/serum/analysis/mharris/tissueSpecific/compareInteract/pro/"

#infile <- "interactParser_corra_target.tsv"
infile <- "interactParser_pro_orbi_ltq.tsv"
#infile <- "interactParser_WT_Mutant.tsv"
RT1_lab <- 'pro_orbi'
RT2_lab <- 'pro_ltq'


setwd(working_dir)
#pandata <- read.table(infile, na.strings = "NA", header=T, sep="\t")
pandata <- read.delim(infile, header=T, sep="\t")



dimnames(pandata)[[2]] <- c('index','peptide',RT1_lab,RT2_lab,'protein_descr1', 'protein_descr2')
attach(pandata)
# Draw distribution graphs


hist(pandata[,RT2_lab],xlab="Retention Times (sec)", main=paste("Distribution of Retention Times of ",RT2_lab), xlim=c(0,6000))

# Create output file and print png of histogram to file
histPlotFilePNG <- paste( "histPlot_",RT1_lab,".png", sep="" ) 
GDD(file = histPlotFilePNG, type = "png", width = 600, height = 400, ps = 12, bg = "white")
hist(pandata[,RT1_lab],xlab="Retention Times (sec)", main=paste("Distribution of Retention Times of ",RT1_lab), xlim=c(0,6000))
dev.off()

# Regression Deletion Analysis

xh <- pandata[ ,RT1_lab] #first set of RTs
yh <- pandata[ ,RT2_lab] #second set of RTs
summary(lmH <- lm(yh ~ xh)) #lmH is the model
(im <- influence.measures(lmH))
which(apply(im$is.inf, 1, any)) # which observations 'are' influential
summary(im) # only these

# Cook's D plot

plot(cookd(lmH))
abline(h=4/length(pandata[,1]),lty=2,col=2)
#identify(1:length(pandata[,1]), cookd(lmH), row.names(pandata))

# Remove the outliers

outlist <- summary(im)

#So far there are three ways to calculate outliers (here are two):
#remove <- which.names(row.names(outlist), pandata) # all influential points
#remove <- influencePlot(lmH,col.identify=4,cex.identify=.7)  # this plot is interactive!

# Manually remove outliers based upon Cook's dist which is contained in outlist

cutoff = 4/(length(pandata[,1])) # cook's dist general cutoff
len <- length(outlist[,1])
tmp <- outlist

j <- 1

for (i in 1:len) {

	# if the cook.d (col 5) is less than the cutoff it is an outliers	
	if (outlist[i,5] < cutoff) { tmp <- tmp[-j,]; j = j-1 }

	j <- j+1
}
outlist <- tmp

remove <- which.names(row.names(outlist), pandata)


# remove now stores the outliers which will be removed before recalc linear
# regression

# Recalculate linear regression model with outliers excluded

remove
lmH.2 <- update(lmH, subset=-remove)
summary(lmH.2)

# To get the peptide list for outliers
pandata[remove,2]

# Replot the graph indicating the outliers

#nf <- layout(matrix(c(1,2), 2, 1, byrow=TRUE), respect=TRUE)
#layout.show(nf)
lowlim <- if( min(xh) < min(yh) ) min(xh)-500 else min(yh)-500
uplim <- if ( max(xh) > max(yh) ) max(xh)+500 else max(yh)+500

# First regression plot, print to file

xyPlotFilePNG <- paste( "regressionPlot_before",".png", sep="" ) 
GDD(file = xyPlotFilePNG, type = "png", width = 600, height = 600, ps = 12, bg = "white")

plot(xh,yh, xlim=c(lowlim,uplim), ylim=c(lowlim,uplim), xlab=paste("Retention Times of ", RT1_lab, " (sec)"), ylab=paste("Retention Times of ", RT2_lab, " (sec)"), main="Regression Analysis of Retention Times", col=1)

xout <-  pandata[remove,3]
yout <- pandata[remove,4]
abline(lmH, lty=1, col=1, lwd=2) # original regression line
abline(lmH.2,lty=2, col=2, lwd=2) # adjusted regression line w/ outliers removed
points(xout, yout, pch=20, col=2)

# Add a legend
#legend(locator(1), lty=1:2, col=1:2, lwd=2, legend=c("All Cases", "Outliers Excluded"))
smartlegend(x = c("left", "center", "right"), 
	y = c("top", "center", "bottom"),
	inset = 0, lty=1:2, col=1:2, lwd=2, 
	legend=c("All Cases", "Outliers Excluded"))

# Get the slope and intercept of new regression

df <- data.frame(t(coef(lmH.2)))
#dimnames(df)[[2]] <- c('b','m') # change column names to b and m (intercept and slope)
#attach(df) # link the column names to the data
b <- df[1,1] # y-intercept
m <- df[1,2] # slope


# Remove outliers from pandata and store in new variable
pandata.2 <- pandata[-remove,]
dimnames(pandata.2)[[2]] <- c('index','peptide',RT1_lab,RT2_lab,'protein_descr1','protein_descr2')
attach(pandata.2)
xh.2 <- pandata.2[,RT1_lab]
yh.2 <- pandata.2[,RT2_lab]

#  Make a file with the outliers info
write.table(pandata[remove,],file=paste("Outliers_",RT1_lab,"_",RT2_lab,".tsv",sep=""),sep="\t",eol="\n",row.names=FALSE)

# The correlation coefficient
r.2 <- cor(xh.2, yh.2)
r2.2 <- r.2^2
#legend(locator(1), paste(c("slope: ","y-intercept: ","R-squared: "), c( round(m,3),round(b,3), round(r2.2,3) )), lty=2, col=2, lwd=2, cex=0.8)

smartlegend(x = c("right"), y = c("bottom"),
	inset = 0, paste(c("slope: ","y-intercept: ","R-squared: "),
	c( round(m,3),round(b,3), round(r2.2,3) )), lty=2, col=4, lwd=2, cex=0.8)

dev.off()




# Make a stat file with slope, yintercept, r-squared, & number outliers
dfstat <- data.frame(cbind(m,b,r2.2,length(remove)))
dimnames(dfstat)[[2]] <- c("slope","yintercept","r_sq","num_outliers")
write.table(dfstat,file=paste("Stats_",RT1_lab,"_",RT2_lab,".tsv",sep=""),sep="\t",eol="\n",row.names=FALSE)


# Shift RTs for data to fit the regression curve equation
ys <- (yh.2-b)/m
lmS <- lm(ys ~ xh.2)
dfS <- data.frame(t(coef(lmS)))
bS <- dfS[1,1] # y-intercept
mS <- dfS[1,2] # slope

# The correlation coefficient
rS <- cor(xh.2, ys)
rS2 <- rS^2


plot(xh.2,ys, xlim=c(lowlim,uplim), ylim=c(lowlim,uplim), xlab=paste("Retention Times of ", RT1_lab, " (sec)"), ylab=paste("Retention Times of ", RT2_lab, " (sec)"), main="Correlation Between Retention Times\n with Correction from Regression Analysis", col=1)
abline(lmS,lty=4,col=4,lwd=2)
abline(lmH.2,lty=1, col=1, lwd=2) # adjusted regression line w/ outliers removed
#legend(locator(1), lty=c(4,1), col=c(4,1), lwd=2, legend=c("Shifted RTs","Original Regression"))

smartlegend(x = c("left", "center", "right"), y = c("top", "center", "bottom"),
	inset = 0, lty=c(4,1), col=c(4,1), lwd=2, 
	legend=c("Shifted RTs","Original Regression"))


#legend(locator(1), paste(c("slope: ","y-intercept: ","R-squared: "), 
#	c( round(mS,3),round(bS,3), round(rS2,3) )), lty=2, col=4, lwd=2, cex=0.8)


smartlegend(x = c("right"), y = c("bottom"),
	inset = 0, paste(c("slope: ","y-intercept: ","R-squared: "),
	c( round(mS,3),round(bS,3), round(rS2,3) )), lty=2, col=4, lwd=2, cex=0.8)


