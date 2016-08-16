## VSLAM.R
# Description: Tool to evaluate visual tracking performance
# Author: Trevor Stanhope

# http://stackoverflow.com/questions/30548217/modifying-a-tukey-hsd-95-family-wise-cl-plot-in-rstudio

## Snippets
#summary(all_surfaces);
#boxplot(data_set$Error ~ data_set$Surface)
#fit = aov(Error ~ Surface, data=data_set)
#summary(fit)
#model.tables(fit, "means") # For type = "means", returns the mean for each term
#hsd <- TukeyHSD(fit, conf.level=0.90) # the confidence leve
#hsd # View table for significant differences within groups
#plot(hsd) # Plot to see results
#par(mfrow=c(3 ,2)) 
#data_set = by_surface[grep("Asphalt", by_surface$Surface), ];
#hist(data_set$Error, border="blue", col="green", xlim=c(100,700), breaks=5);

## Cleanup
rm(list = ls())

## Libraries
library("multcomp")
library("hydroGOF")
library(hexbin)

## Functions
printf <- function(...) { cat(sprintf(...)) } 
loss1 <- function(pars, xObs, yObs) {
  b0 <- pars[1]
  b1 <- pars[2]
  sigma <- pars[3]
  yFit <- b0 + b1 * xObs
  e <- yObs - yFit
  nll <- -sum(dnorm(e, 0, sigma, log=T))
  return(nll)
}

## Constants
NUM_REPLICATES = 5;
DIR = "C:\\Users\\Bioresource\\Documents\\Trevor Docs\\Repos\\project-R";
SURFACES = c("asphault", "gravel", "residue", "grass", "corn", "hay");
SURFACES_FIXED = c("Asphalt", "Gravel", "Seedlings", "Grass", "Corn", "Pasture");
HEIGHT = c("Low", "Low", "Medium", "Medium", "High", "High");
THRESHOLDS = c(125, 250, 500, 750, 1000);
FILE_EXTENSION = ".csv";
CONFIDENCE = 0.95;
PERCENTILE = 0.95;

## Workspace Environment
setwd(DIR);

## Iterate All
data_pct = data.frame(Surface=c(), Freq=c(), Error=c(), Height=c());
data_rmse = data.frame(Surface=c(), Freq=c(), Error=c(), Height=c());
for (k in 1:length(THRESHOLDS)) {
  threshold = THRESHOLDS[k];
  for (j in 1:length(SURFACES)) {
    surface = SURFACES[j];
    height = HEIGHT[j];
    surface_fixed = SURFACES_FIXED[j];
    for (i in 1:NUM_REPLICATES) {
      rtk_path = paste(DIR, "\\data\\RTK\\", sep="");
      cv_path = paste(DIR, "\\data\\2NN\\", threshold, "\\", sep="");
      trial_name = paste(surface, i, sep="-");
      file_name = paste(trial_name, FILE_EXTENSION, sep="");
      rtk_data = read.csv(paste(rtk_path, file_name, sep=""));
      cv_data = read.csv(paste(cv_path, file_name, sep=""));
      speed = rtk_data['speed'];
      error_v = speed - cv_data['v']; # find the error between the two sensors
      hz = cv_data['hz'];
      n = lengths(speed);
      error_pct = quantile(abs(error_v), PERCENTILE, na.rm = TRUE);
      error_rmse = rmse(speed, cv_data['v']);
      freq_rmse = quantile(abs(hz), PERCENTILE, na.rm = TRUE);
      freq_pct = rmse(speed, hz);
      data_pct = rbind(data_pct, data.frame(Error=error_pct, Surface=surface_fixed, Threshold=threshold, Freq=freq_rmse, Height=height));
      data_rmse = rbind(data_rmse, data.frame(Error=error_rmse, Surface=surface_fixed, Threshold=threshold, Freq=freq_rmse, Height=height));
    }
  }
}

## Iterate By Surface (Fixed Threshold)
threshold = 500;
by_surface = data.frame(Surface=c(), Freq=c(), Error=c(), Speed=c());
for (j in 1:length(SURFACES)) {
  surface = SURFACES[j];
  surface_fixed = SURFACES_FIXED[j];
  for (i in 1:NUM_REPLICATES) {
    rtk_path = paste(DIR, "\\data\\RTK\\", sep="");
    cv_path = paste(DIR, "\\data\\1NN\\", threshold, "\\", sep="");
    trial_name = paste(surface, i, sep="-");
    file_name = paste(trial_name, FILE_EXTENSION, sep="");
    rtk_data = read.csv(paste(rtk_path, file_name, sep=""));
    cv_data = read.csv(paste(cv_path, file_name, sep=""));
    speed = rtk_data['speed'];
    error_v = speed - cv_data['v']; # find the error between the two sensors
    hz = cv_data['hz'];
    n = lengths(speed);
    trial = data.frame(Error=unname(error_v), Speed=unname(speed), Freq=unname(hz), Surface=rep(surface_fixed, n))
    #plot(trial$Speed, trial$Error, xlab="True Speed (km/h)", ylab="Error (km/h)", ylim=c(-1,1));
    by_surface = rbind(by_surface, trial);
  }
}

## Iterate By Threshold (Fixed Surface)
surface = "residue"
by_threshold = data.frame(Threshold=c(), Error=c(), Freq=c(), Speed=c());
for (j in 1:length(THRESHOLDS)) {
  threshold = THRESHOLDS[j];
  for (i in 1:NUM_REPLICATES) {
    cv_path = paste(DIR, "\\data\\1NN\\", threshold, "\\", sep="");
    rtk_path = paste(DIR, "\\data\\RTK\\", sep="");
    trial_name = paste(surface, i, sep="-");
    file_name = paste(trial_name, FILE_EXTENSION, sep="");
    rtk_data = read.csv(paste(rtk_path, file_name, sep=""));
    cv_data = read.csv(paste(cv_path, file_name, sep=""));
    speed = rtk_data['speed'];
    error = speed - cv_data['v']; # find the error between the two sensors
    hz = cv_data['hz'];
    n = lengths(hz);
    trial = data.frame(Error=unname(error), Speed=unname(speed), Freq=unname(hz), Threshold=rep(toString(threshold), n));
    trial.lm = lm(Error ~ Speed, data=trial);
    trial.res = resid(trial.lm);
    plot(trial.res);
    by_threshold = rbind(by_threshold, trial); # append trial to data set
  }
}

## Analyze Error ~ Surface
boxplot(by_surface$Error ~ by_surface$Surface, show.names=TRUE, ylim=c(-10, 10), xlab="Surface", ylab="Error (km/h)");
fit_surface = aov(Error ~ Surface, data=by_surface);
summary(fit_surface);
model.tables(fit_surface, "means") # For type = "means", returns the mean for each term
hsd_surface <- TukeyHSD(fit_surface, conf.level=CONFIDENCE); # the confidence level 
hsd_surface # View table for significant differences within groups
plot(hsd_surface); # Plot to see results
lm_error_surface = lm(Error ~ Surface, data=by_surface);
summary(lm_error_surface);

## Analyze Error ~ Threshold
boxplot(by_threshold$Error ~ by_threshold$Threshold, show.names=TRUE, ylim=c(-10, 10), xlab="Threshold (Keypoints per Image)", ylab="Error (km/h)");
fit_error_threshold = aov(Error ~ Threshold, data=by_threshold);
summary(fit_error_threshold);
model.tables(fit_error_threshold, "means"); # For type = "means", returns the mean for each term
hsd_error_threshold <- TukeyHSD(fit_error_threshold, conf.level=CONFIDENCE); # the confidence level 
hsd_error_threshold # View table for significant differences within groups
plot(hsd_error_threshold); # Plot to see results
lm_error_threshold = lm(Error ~ Threshold, data=by_threshold);
summary(lm_error_threshold);

## Analyze Frequency ~ Surface 
boxplot(by_surface$Freq ~ by_surface$Surface, show.names=TRUE, xlab="Surface", ylab="Frequency (Hz)", ylim=c(0,40));
fit_freq_surface = aov(Freq ~ Surface, data=by_surface);
summary(fit_freq_surface);
model.tables(fit_freq_surface, "means"); # For type = "means", returns the mean for each term
hsd_freq_surface <- TukeyHSD(fit_freq_surface, conf.level=CONFIDENCE); # the confidence level 
hsd_freq_surface # View table for significant differences within groups
plot(hsd_freq_surface); # Plot to see results
lm_freq_surface = lm(Freq ~ Surface, data=by_surface);
summary(lm_freq_surface);

## Analyze Frequency ~ Threshold
boxplot(data_rmse$Freq ~ data_rmse$Threshold, show.names=TRUE, xlab="Threshold (Keypoints per Image)", ylab="Frequency (Hz)");
fit_freq_threshold = aov(Freq ~ Threshold, data=data_rmse);
summary(fit_freq_threshold);
model.tables(fit_freq_threshold, "means"); # For type = "means", returns the mean for each term
hsd_freq_threshold <- TukeyHSD(fit_freq_threshold, conf.level=CONFIDENCE); # the confidence level 
hsd_freq_threshold # View table for significant differences within groups
plot(hsd_freq_threshold); # Plot to see results
lm_freq_threshold = lm(Freq ~ Threshold, data=data_rmse);
summary(lm_freq_threshold);

## RMSE ~ Surface
boxplot(data_rmse$Error ~ data_rmse$Surface, xlab="Surface", ylab="95th Percentile (km/h)", ylim=c(0,15))
fit_rmse_threshold = aov(Error ~ Surface + Threshold, data=data_rmse);
summary(fit_rmse_threshold);
model.tables(fit_rmse_threshold, "means"); # For type = "means", returns the mean for each term
hsd_rmse_threshold <- TukeyHSD(fit_rmse_threshold, conf.level=CONFIDENCE); # the confidence level 
hsd_rmse_threshold # View table for significant differences within groups
plot(hsd_rmse_threshold); # Plot to see results
lm_rmse_threshold = lm(Error ~ Threshold, data=by_threshold);
summary(lm_rmse_threshold);

## RMSE ~ Height
boxplot(data_rmse$Error ~ data_rmse$Height, xlab="Height", ylab="RMSE (km/h)", ylim=c(0,5))
fit_rmse_height = aov(Error ~ Height * Threshold, data=data_rmse);
summary(fit_rmse_height);
model.tables(fit_rmse_height, "means"); # For type = "means", returns the mean for each term
hsd_rmse_height <- TukeyHSD(fit_rmse_height, conf.level=CONFIDENCE); # the confidence level 
hsd_rmse_height # View table for significant differences within groups
plot(hsd_rmse_height); # Plot to see results
lm_rmse_height = lm(Error ~ Height, data=data_rmse);
summary(lm_rmse_height);

## Two-Way ANOVAs
int_mult = aov(Error ~ Surface * Threshold, data=data_rmse);
summary(int_mult);
int_add = aov(Error ~ Surface + Threshold, data=data_rmse);
summary(int_add);
no_int = aov(Error ~ Surface, data=data_rmse);
summary(no_int);

hydroGOhydroGO

## Examples
data_subset = by_surface[grep("Asphalt", by_surface$Surface), ];
plot(data_subset$Speed, data_subset$Error, xlab="True Speed (km/h)", ylab="Error (km/h)", ylim=c(-1,1))
data_subset = by_surface[grep("Gravel", by_surface$Surface), ];
plot(data_subset$Speed, data_subset$Error, xlab="True Speed (km/h)", ylab="Error (km/h)", ylim=c(-1,1))