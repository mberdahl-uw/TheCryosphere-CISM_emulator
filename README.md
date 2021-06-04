# EmulatorPaper
Data and Code for Emulator Paper

There are three R scripts, and some CSV files needed to run these.  Together these scripts validate & build the emulator, and then finally run example priors through the emulator. The RData files are the actual emulators (year 100 and 200) that are built and used by the final (third) script.

There are three R scripts that should be run in the following order:
First, cism_slr_emulator_valid.R:
1) Validate emulator with 80% of the data. It saves the validation emulator as an R script, but we dont really use it again.
This script also plots the 1:1 plot of the true vs. predicted values for 2100 and 2200. You can change the year you are interested in the script.

Next, CISM_slr_emulator.R:
2) Build the emulator and save it as an Rdata file.  This uses the entire CISM ensemble to build the emulator, and its only output right now is the Rdata file.  Can choose year 2100 or 2200.

And finally, prior_predictive.R
3) Use the emulator saved in step two above to run a set of priors to predict the SLR distribution as a PDF for 2100 or 2200. This script outputs the SLR predicted data as a csv, and also makes a plot of the PDF of the predicted SLRs for the given priors.
