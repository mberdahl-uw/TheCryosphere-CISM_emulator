##############################################################
## This code was originally written by Nathan, Jan 2020
## I've made some small changes and commenting it as best I can.
## It is the FIRST script in the following order of scripts:
##1) Validate emulator with 80% of the data. It saves the validation emulator as an R script, but we dont really use it again.
##  This script also plots the 1:1 plot of the true vs. predicted values for 2100 and 2200. You can cahnge the year you are interested in in the script.
##2) Build the emulator and save it as an Rdata file.  This uses the entire CISM ensemble to build the emulator, and its only output right now is the Rdata file.  Can choose 2100 or 2200.
##3) Use the emulator saved in step two above to run a set of priors to predict the SLR distribution as a PDF for 2100 or 2200. This script outputs the SLR predicted data as a csv, and also makes a plot of the PDF of the predicted SLRs for the given priors.
############################################################

#uncomment this the first time you run in R so that it installs the appropriate package
# install.packages("GPfit")
library(GPfit)

setwd("~/Documents/UW/CMIP_meltrates/Timmermann/PriorWork/")

# read the design metadata for the CISM ensemble, basin by basin.
amundsen = read.csv("design_metadata_amundsen.csv")
eais = read.csv("design_metadata_eais.csv")
peninsula = read.csv("design_metadata_peninsula.csv")
ross = read.csv("design_metadata_ross.csv")
weddell = read.csv("design_metadata_weddell.csv")

sigpars = c("m2200","t0","tau")

design = cbind(amundsen[sigpars], eais[sigpars], peninsula[sigpars], ross[sigpars], weddell[sigpars])
succeeded = amundsen$flag == 0
design = design[succeeded,]

normalize = function(x) (x-min(x))/(max(x)-min(x))
design = as.data.frame(lapply(design, normalize))

# read all the SLR data from CISM ensemble
slr.data = read.csv("SLR.csv", comment.char="#", header=FALSE)
slr.data = slr.data[succeeded,2:ncol(slr.data)]

slr = data.frame(matrix(NA, nrow=nrow(slr.data), ncol=ncol(slr.data)))
for(i in 1:nrow(slr))
	slr[i,] = -as.numeric(slr.data[i,] - slr.data[i,1])
yr = 2200 # choose the year to validate on
baseyr = 2000; iyr = yr - baseyr
slr.yr = slr[,iyr]
N = length(slr.yr)

idx = sample(1:N)
frac_train = 0.8  #train on only a fraction of the data, test with the rest.
Ntrain = floor(frac_train * N)
Ntest = N - Ntrain
itrain = idx[1:Ntrain]
itest = idx[(Ntrain+1):N]

design.train = design[itrain,]
design.test = design[itest,]
slr.yr.train = slr.yr[itrain]
slr.yr.test = slr.yr[itest]

print(dim(design.train))
print(dim(design.test))

t0 = Sys.time()
emu = GP_fit(design.train, slr.yr.train)
print(Sys.time()-t0)  #prints the time it took to run the code.

predict = function(emu, newx) {
	p = predict.GP(emu, xnew=newx)
	pred = t(as.matrix(p$Y_hat))
	pred.error = sqrt(p$MSE)
	return(list(pred=pred, pred.error=pred.error))
}

pred = predict(emu, design.test)
pred.test = pred$pred
pred.test.sd = pred$pred.error

save.image(sprintf("emu_slr%d_valid.RData", yr)) #save the validation emulator (built on only 80% of the ensemble SLR data) to an Rdata Structure

# plot the 1:1 line of emulator validation (true vs predicted) with the points and error bars associated.
plot(pred.test, slr.yr.test, pch=19, col="blue", main=sprintf("SLR %d emulator validation", yr), xlab="Predicted", ylab="True", xlim=c(0,500), ylim=c(0,500), frame=FALSE)
segments(pred.test-2*pred.test.sd, slr.yr.test, pred.test+2*pred.test.sd, slr.yr.test, col="blue")
lines(slr.yr.train, slr.yr.train, lty="dashed", lwd=2)
