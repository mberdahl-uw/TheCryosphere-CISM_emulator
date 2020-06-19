##############################################################
## This code was originally written by Nathan, Jan 2020
## I've made some small changes and commenting it as best I can.
## It is the THIRD script in the following order of scripts:
##1) Validate emulator with 80% of the data. It saves the validation emulator as an R script, but we dont really use it again.
##  This script also plots the 1:1 plot of the true vs. predicted values for 2100 and 2200. You can cahnge the year you are interested in in the script.
##2) Build the emulator and save it as an Rdata file.  This uses the entire CISM ensemble to build the emulator, and its only output right now is the Rdata file.  Can choose 2100 or 2200.
##3) Use the emulator saved in step two above to run a set of priors to predict the SLR distribution as a PDF for 2100 or 2200. This script outputs the SLR predicted data as a csv, and also makes a plot of the PDF of the predicted SLRs for the given priors.
############################################################


install.packages("GPfit")
#uncomment this the first time you run in R so that it installs the appropriate package.
library(GPfit)

setwd("~/Documents/UW/CMIP_meltrates/Timmermann/PriorWork/")


##### prior PDF and bounds on sigmoid parameters

### normal priors####################################
#####################################################
# these values are based on just the A1B scenarios from Cornford and Timmermann
# the excel spreadsheet named BestFitParameters_Cornford_Timmermann.xlsx has all of the actual data that leads to these
# mean and std values.  For ease, these are currently hardcoded values.
#####################################################
m2200.amundsen.min = 0.1; m2200.amundsen.max = 50.0
m2200.eais.min = 0.1; m2200.eais.max = 36.0
m2200.peninsula.min = 0.1; m2200.peninsula.max = 12.0
m2200.ross.min = 0.1; m2200.ross.max = 20.0
m2200.weddell.min = 0.1; m2200.weddell.max = 16.0


## old values, based on BEST fits.

#m2200.amundsen.mean = 14.33446548; m2200.amundsen.sd = 8.19604831
#m2200.eais.mean = 11.03972597; m2200.eais.sd = 4.03074875
#m2200.peninsula.mean = 1.59706821; m2200.peninsula.sd = 1.59706821
#m2200.ross.mean = 2.64408506; m2200.ross.sd = 2.67049618
#m2200.weddell.mean = 4.5640759; m2200.weddell.sd = 1.45733891

#t0.amundsen.mean = 153.6 ; t0.amundsen.sd = 57.68396658
#t0.eais.mean = 192.0; t0.eais.sd = 53.14132102
#t0.peninsula.mean = 100.5; t0.peninsula.sd = 0.5
#t0.ross.mean = 181.66666667; t0.ross.sd = 57.76004001
#t0.weddell.mean = 193.5; t0.weddell.sd = 52.82754963

#tau.amundsen.mean = 34.6 ; tau.amundsen.sd = 18.39130229
#tau.eais.mean = 34.25; tau.eais.sd = 7.85413904
#tau.peninsula.mean = 36.0; tau.peninsula.sd = 26.0
#tau.ross.mean = 23.33333333; tau.ross.sd = 9.46337971
#tau.weddell.mean = 29.0; tau.weddell.sd = 13.05756486


## new values based on adjusted fits so taht we're not edging.
m2200.amundsen.mean = 12.58326; m2200.amundsen.sd = 6.663675
m2200.eais.mean = 12.39287; m2200.eais.sd = 5.499974
m2200.peninsula.mean = 1.70747; m2200.peninsula.sd = 1.70747
m2200.ross.mean = 2.848176; m2200.ross.sd = 2.426784
m2200.weddell.mean = 4.98193; m2200.weddell.sd = 1.059854

t0.amundsen.mean = 132.3363 ; t0.amundsen.sd = 27.36718
t0.eais.mean = 163.3008; t0.eais.sd = 24.82787
t0.peninsula.mean = 162.0785; t0.peninsula.sd = 12.9215
t0.ross.mean = 163.4154; t0.ross.sd = 20.2603
t0.weddell.mean = 162.4607; t0.weddell.sd = 35.32445

tau.amundsen.mean = 26.37856; tau.amundsen.sd = 20.00555
tau.eais.mean =38.016; tau.eais.sd = 19.32499
tau.peninsula.mean = 47.5055; tau.peninsula.sd = 22.5055
tau.ross.mean = 32.63466; tau.ross.sd = 18.66036
tau.weddell.mean = 24.28617; tau.weddell.sd = 6.269313


### put these all together in the following order: Amundsen, EAIS, Peninsula, Ross, Weddell
par.mean = c(m2200.amundsen.mean, t0.amundsen.mean, tau.amundsen.mean, m2200.eais.mean, t0.eais.mean, tau.eais.mean, m2200.peninsula.mean, t0.peninsula.mean, tau.peninsula.mean, m2200.ross.mean, t0.ross.mean, tau.ross.mean, m2200.weddell.mean, t0.weddell.mean, tau.weddell.mean)
par.sd = c(m2200.amundsen.sd, t0.amundsen.sd, tau.amundsen.sd, m2200.eais.sd, t0.eais.sd, tau.eais.sd, m2200.peninsula.sd, t0.peninsula.sd, tau.peninsula.sd, m2200.ross.sd, t0.ross.sd, tau.ross.sd, m2200.weddell.sd, t0.weddell.sd, tau.weddell.sd)

########################
### parameter bounds####
########################

t0.min = 100; t0.max = 225
tau.min = 10; tau.max = 75
idx.basin = seq(0, 14, by=3)

par.min = rep(NA, 15)
par.min[1+idx.basin] = c(m2200.amundsen.min, m2200.eais.min, m2200.peninsula.min, m2200.ross.min, m2200.weddell.min)
par.min[2+idx.basin] = t0.min
par.min[3+idx.basin] = tau.min

par.max = rep(NA, 15)
par.max[1+idx.basin] = c(m2200.amundsen.max, m2200.eais.max, m2200.peninsula.max, m2200.ross.max, m2200.weddell.max)
par.max[2+idx.basin] = t0.max
par.max[3+idx.basin] = tau.max

#######################################
##### sample from truncated normal prior
########################################

sample.prior = function(par.mean, par.sd) {
	p = rnorm(length(par.mean), par.mean, par.sd)
	p[p<par.min] = par.min[p<par.min]
	p[p>par.max] = par.max[p>par.max]
	return(p)
}

#########################
##### load SLR emulator
#########################
yr = 2100  #choose the year we are interested in.
load(sprintf("emu_slr%d.RData", yr))  # this was generated in the script cism_slr_emulator.R
design.norm = design

# these give the basin specific metadata########
amundsen = read.csv("design_metadata_amundsen.csv")
eais = read.csv("design_metadata_eais.csv")
peninsula = read.csv("design_metadata_peninsula.csv")
ross = read.csv("design_metadata_ross.csv")
weddell = read.csv("design_metadata_weddell.csv")

### the sigmoid parameters
sigpars = c("m2200","t0","tau")

design = cbind(amundsen[sigpars], eais[sigpars], peninsula[sigpars], ross[sigpars], weddell[sigpars])
succeeded = amundsen$flag == 0
design = design[succeeded,]

dmin = apply(design, 2, min)
dmax = apply(design, 2, max)
normalize.input = function(x) {
	# first, scale the first parameter (M2200) in each basin to [0,1] by the basin maxima
	x[1+idx.basin] = x[1+idx.basin] / par.max[1+idx.basin]
	# then scale everything by the training ensemble range
	return((x-dmin)/(dmax-dmin))
}

predict = function(emu, newx) {
	if(is.vector(newx)) newx = t(as.matrix(newx))
	for(i in 1:nrow(newx))
		newx[i,] = normalize.input(newx[i,])
	p = predict.GP(emu, xnew=newx)
	pred = t(as.matrix(p$Y_hat))
	pred.error = sqrt(p$MSE)
	return(list(pred=pred, pred.error=pred.error))
}

##### sigmoid #################

sigmoid = function(t, m2200, t0, tau) {
	mmax = m2200 * (1 + exp(-((2200-2000)-t0)/tau))
	traj = mmax / (1 + exp(-(t-t0)/tau))
	return(traj)
}

######################################
##### propagate prior through emulator
########################################

parnames = c("M2200.amundsen", "t0.amundsen", "tau.amundsen", "M2200.eais", "t0.eais", "tau.eais", "M2200.peninsula", "t0.peninsula", "tau.peninsula", "M2200.ross", "t0.ross", "tau.ross", "M2200.weddell", "t0.weddell", "tau.weddell")

prior.predictive = function(Nsamp, p.mean, p.sd) {
	t = 2000:2200
	tstd = t - t[1]
	
	par = matrix(nrow=Nsamp, ncol=15)
	colnames(par) = parnames
	
	melt.amundsen = matrix(nrow=Nsamp, ncol=length(t))
	melt.eais = matrix(nrow=Nsamp, ncol=length(t))
	melt.peninsula = matrix(nrow=Nsamp, ncol=length(t))
	melt.ross = matrix(nrow=Nsamp, ncol=length(t))
	melt.weddell = matrix(nrow=Nsamp, ncol=length(t))

	slr = rep(NA, Nsamp)
	slr.err = rep(NA, Nsamp)	
	for(i in 1:Nsamp) {
		p = sample.prior(p.mean, p.sd)
		par[i,] = p
		
		melt.amundsen[i,] = sigmoid(tstd, p[1], p[2], p[3])
		melt.eais[i,] = sigmoid(tstd, p[4], p[5], p[6])
		melt.peninsula[i,] = sigmoid(tstd, p[7], p[8], p[9])
		melt.ross[i,] = sigmoid(tstd, p[10], p[11], p[12])
		melt.weddell[i,] = sigmoid(tstd, p[13], p[14], p[15])
	}
			
	emu.pred = predict(emu, par)
	slr = emu.pred$pred
	slr.err = emu.pred$pred.error
	    
	return(list(par=par, melt.amundsen=melt.amundsen, melt.eais=melt.eais, melt.peninsula=melt.peninsula, melt.ross=melt.ross, melt.weddell=melt.weddell, slr=slr, slr.err=slr.err))
}

######################################

Nsamp = 1000
pred = prior.predictive(Nsamp, par.mean, par.sd)

slr = as.vector(pred$slr)

par(mar=c(5,5,2.5,1))
# make a plot of the PDF of SLR predicted by the emulator.
plot(density(slr, from=0, to=1.1*max(slr)), lwd=5, col="red", frame=FALSE, xlab="Sea level rise (mm)", ylab="Probability density", main=sprintf("Antarctic sea level rise projection in %d", yr), cex.main=1.7, cex.lab=1.6, cex.axis=1.4)



# save the output######################
write.csv(slr, 'slr_emu_2100.csv')
