##############################################################
## This code was originally written by Nathan, Jan 2020
## I've made some small changes and commenting it as best I can.
## It is the SECOND script in the following order of scripts:
##1) Validate emulator with 80% of the data. It saves the validation emulator as an R script, but we dont really use it again.
##  This script also plots the 1:1 plot of the true vs. predicted values for 2100 and 2200. You can cahnge the year you are interested in in the script.
##2) Build the emulator and save it as an Rdata file.  This uses the entire CISM ensemble to build the emulator, and its only output right now is the Rdata file.  Can choose 2100 or 2200.
##3) Use the emulator saved in step two above to run a set of priors to predict the SLR distribution as a PDF for 2100 or 2200. This script outputs the SLR predicted data as a csv, and also makes a plot of the PDF of the predicted SLRs for the given priors.
############################################################

#uncomment this the first time you run in R so that it installs the appropriate package
# install.packages("GPfit")

library(GPfit)

setwd("~/Documents/UW/CMIP_meltrates/Timmermann/PriorWork/")

yr = 2200  #choose the year we are interested in building an emulator for.

train.emulator = TRUE  #set this to true if you want to train the emulator now.  If you've already done this before and have it saved as an RData file, then you can say this is false and it'll jump straight to the else line near the end of the file.
if(train.emulator) {
    
## load the design metadata for each basin.
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

# read the SLR data from CISM ensemble
slr.data = read.csv("SLR.csv", comment.char="#", header=FALSE)
slr.data = slr.data[succeeded,2:ncol(slr.data)]

slr = data.frame(matrix(NA, nrow=nrow(slr.data), ncol=ncol(slr.data)))
for(i in 1:nrow(slr))
	slr[i,] = -as.numeric(slr.data[i,] - slr.data[i,1])
baseyr = 2000; iyr = yr - baseyr
slr.yr = slr[,iyr]
N = length(slr.yr)

design.train = design
slr.yr.train = slr.yr

print(dim(design.train))  # print dimensions.

t0 = Sys.time()
emu = GP_fit(design.train, slr.yr.train)
print(Sys.time()-t0)

save.image(sprintf("emu_slr%d.RData", yr))  #save the emulator as an Rdata structrue
} else {
    load("emulator_slr%d.RData", yr)  # if the if statement at the top was false, then just load the emulator here instead of building it all in the if statement
}

predict = function(emu, newx) {
	p = predict.GP(emu, xnew=newx)
	pred = t(as.matrix(p$Y_hat))
	pred.error = sqrt(p$MSE)
	return(list(pred=pred, pred.error=pred.error))
}
