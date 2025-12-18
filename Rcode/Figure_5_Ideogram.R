#########################################################################
####          Figure 5 : Ideogram                                    ####
#### Goal: Understand the distribution of the SNPs in the genome     ####
#########################################################################


#Ideogram
# setwd("path//to//working//dir")

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("GenomicRanges")
# 
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("regioneR", force = TRUE)

library(karyoploteR)
library(GenomicRanges)
library(regioneR)
library(readxl)
library(openxlsx2)
library(dplyr)

#Import the data
data <- load("Rdata/Figure_5_Ideogram_input.RData") # All locus information required for plotting

head(Loc_all)

#create seqinfo object
chr_names <- as.character(chromosome$CHR)
chr_length <- chromosome$chromEnd

seqinfo_obj <- Seqinfo(
  seqnames = chr_names,
  seqlengths = chr_length,
  isCircular = rep(FALSE, length(chr_names))
)


# Define the chromosome names and ranges
chr <- as.character(chromosome$CHR)
start <- c(chromosome$chromStart)
end <- c(chr_length)

# Create the GRanges object
gr_obj <- GRanges(seqnames = Rle(chr), ranges = IRanges(start = start, end = end))


# Create the GRanges object for locus info
data.points <- GRanges(seqnames = Rle(Loc_all$CHR), ranges = IRanges(start = Loc_all$Start, end = Loc_all$End))
data.points$y <- 0.3
data.points$CWD_captive <- Loc_all$CWD_marker_y

# Create the karyotype plot
kp <- plotKaryotype(genome = gr_obj, ideogram.plotter=kpAddCytobands)
# Add background color to data panel 1
kpDataBackground(kp, data.panel = 1, r0=0, r1=0.80, col="#fefeea")

# Plot points with specified aesthetics
kpPoints(kp, data=data.points, pch=25, cex=0.8, bg="grey75", col ="grey12", r0=0, r1=0.8, alpha = 0.6)
kpPoints(kp, data=data.points, y= data.points$CWD_captive, pch=23, cex=1.5, col="white", bg ="blue4",
         alpha = 0.5, r0=0, r1=0.8)

# Details added using an editor


# 
# #Save high res figure
# library(ggplot2)
# 
# #save the plot
# # Define the file path where the plot should be saved
# file_path <- "./Ideogram.png"
# 
# # Open a PNG device with specified dimensions and resolution
# png(filename = file_path, width = 2700, height = 1950, res = 300)
# 
# # Create the karyotype plot
# kp <- plotKaryotype(genome = gr_obj, ideogram.plotter = kpAddCytobands)
# 
# # Add background color to data panel 1
# kpDataBackground(kp, data.panel = 1, r0 = 0, r1 = 0.80, col = "#fefeea")
# 
# # Plot points with specified aesthetics
# kpPoints(kp, data = data.points, pch = 25, cex = 0.8, bg = "#00254d", col = "#00254d", r0 = 0, r1 = 0.8)
# kpPoints(kp, data = data.points, pch = 25, cex = 0.8, bg = "#ADD8E6", col = "#00254d", r0 = 0, r1 = 0.8, alpha = 0.6)
# kpPoints(kp, data = data.points, y = data.points$CWD_captive, pch = 25, cex = 0.8, col = "#FF4433", bg = "#F88379", r0 = 0, r1 = 0.8)
# 
# # Close the graphics device to save the plot
# dev.off()
















