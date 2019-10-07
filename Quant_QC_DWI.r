#! /opt/quarantine/R/3.4.3/build2/bin/Rscript

#This script conducts the quantitative QC of diffusion weighted images

args <- commandArgs(TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (csv)", call.=FALSE)
} 
# else if (length(args)==1) { #Ccan use something like this so the SD can be changeable
#   # default SD size is 2
#   args[2] = 2
# }

metric_file = args[1]
# SD_size = args[2]


#loading necessary libraries
# library("zoo", lib.loc="/scratch/nforde/homotopic/bin/R_lib")
# library("lmtest", lib.loc="/scratch/nforde/homotopic/bin/R_lib")
# library("forecast", lib.loc="/scratch/nforde/homotopic/bin/R_lib")

library(tidyverse)
library(cowplot)
library(ggpubr)
library(gridExtra)
library(grid)
library(rmarkdown)


library(glue)
library(plotly)
library(table1)
library(wesanderson)

library(gridGraphics)

## set paths & load csv
outdir <- getwd()
df <- read.csv(metric_file, header=TRUE) 


#bvals can vary across studies and there are multiple for multishell data. This will make a list of the cnr variables you need 
avg_cnr <- df %>% select(starts_with("avg_cnr")) %>% names()

plot_cnr <- function(roi, df) {
  df2 <- df %>% filter(get(roi) < 9999)
  plot <- ggplot(df2, aes_string(x = roi)) +
    geom_density(size = 2, aes(color = "#440154FF"), fill="lightblue") +
    geom_vline(data = df2, aes(xintercept = mean(get(roi)), color = "lightblue"), linetype = "dashed", size = 1) +
    labs(title = paste0("Distribution of ", roi," CNR Values"), x = paste0(roi, " CNR Values"), y = "Density") +
    theme_bw() + theme(legend.position="none")
  return(plot)
}

density_plot_cnr <- lapply(avg_cnr, function(f) plot_cnr(f, df)) # this should create a list of 1 or more plots depending on the number of bvals

density_plot_snr <- function(df) {
  df2 <- df %>% filter(avg_snr_0 < 9999)
  plot <- ggplot(df2, aes(x = avg_snr_0)) +
    geom_density(size = 2, aes(color = "#440154FF"), fill = "lightblue") +
    geom_vline(data = df2, aes(xintercept = mean(avg_snr_0), color = "lightblue"), linetype = "dashed", size = 1) +
    labs(title = "Distribution of SNR Values", x = "SNR Values", y = "Density") +
    theme_bw() + theme(legend.position="none")
  return(plot)
}


density_plot_rel_mot <- function(df) {
  plot <- ggplot(df, aes(x = qc_mot_rel)) +
    geom_density(size = 2, aes(color = "#440154FF"), fill = "lightblue") +
    geom_vline(data = df, aes(xintercept = mean(qc_mot_rel), color = "lightblue"), linetype = "dashed", size = 1) +
    labs(title = "Distribution of Relative Motion Values", x = "Relative Motion Values", y = "Density") +
    theme_bw() + theme(legend.position="none")
  return(plot)
}

density_plot_outliers <- function(df) {
  plot <- ggplot(df, aes(x = qc_outliers_tot)) +
    geom_density(size = 2, aes(color = "#440154FF"), fill = "lightblue") +
    geom_vline(data = df, aes(xintercept = mean(qc_outliers_tot), color = "lightblue"), linetype = "dashed", size = 1) +
    labs(title = "Distribution of Outlier Values", x = "Outlier Values", y = "Density") +
    theme_bw() + theme(legend.position="none")
  return(plot)
}

plot_grid(density_plot_outliers(df), density_plot_rel_mot(df), density_plot_snr(df) + rremove("x.text"), 
          labels = c("A", "B", "C"),
          ncol = 1, nrow = 3)

plot_grid(plotlist=density_plot_cnr, ncol = 1)


#Excluding based on bad CNR

bad_cnr_1000 <- df %>% filter(avg_cnr_1000 < 9999) %>%
  mutate(cnr_std = scale(avg_cnr_1000)) %>% #This provides us with their STD so that we can determine if the distribution is Guassian
  filter(cnr_std < -2 | cnr_std > 2) #this will inform you of the subjects who had CNR that was greater and less than 2 standard deviations away from the mean. Anyone less than 2 STD should be excluded.
  #Those with 2 STD greater than the mean should be checked to make sure their other metrics are in the normal range and their visual QC is alright, however, they shouldn't be excluded unless there are other problems with their acquisitions.

#Note: if its multishell then you will have different average CNRs for every B values, and in that case, be sure to check all of them
#In the case of multishell acquisitions, use this if loop and feel free to change the B values to suit your study specific acquisitions

bad_cnr_1600 <- df %>%
  mutate(cnr_std_1600 = scale(avg_cnr_1600)) %>% 
  filter(cnr_std_1600 < -2 | cnr_std_1600 > 2)

bad_cnr_2600 <- df %>%
  mutate(cnr_std_2600 = scale(avg_cnr_2600)) %>% 
  filter(cnr_std_2600 < -2 | cnr_std_2600 > 2)


#Excluding based on bad SNR
bad_snr <- df %>%
  mutate(snr_std = scale(avg_snr_0)) %>%
  filter(snr_std < -2 | snr_std > 2)


#Excluding based on too much relative motion
bad_rel_mot <- df %>%
  mutate(rel_mot_std = scale(qc_mot_rel)) %>%
  filter(rel_mot_std > 2) #for motion, we're only concerned about participants with too much motion


#Excluding based on too many outliers
bad_outliers <- df %>%
  mutate(outliers_std = scale(qc_outliers_tot)) %>%
  filter(outliers_std > 2)  #for outliers, we're only concerned about participants with too many outliers


#Creating a dataframe all the subjects that failed the Quantitative QC and their problematic values

library(plyr)
Failed_Quant_QC <- rbind.fill(bad_cnr, bad_cnr_2600, bad_cnr_1600, bad_snr, bad_rel_mot, bad_outliers)

Failed_Quant_QC_DF <- Failed_Quant_QC[,c("subject_id", "avg_cnr_1000", "avg_cnr_1600", "avg_cnr_2600", "avg_snr_0", "qc_mot_rel", "qc_outliers_tot")]

library(kableExtra)
DWI.Metrics.Table <- as.data.table(Failed_Quant_QC_DF, keep.rownames = TRUE)
DWI.Metrics.Table %>%
  kable() %>%
  kable_styling()




### this line graph will plot the relationship between age
