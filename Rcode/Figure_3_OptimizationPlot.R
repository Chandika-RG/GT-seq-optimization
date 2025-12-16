#####################################################################################
####          Script used to generate optimization plots (Figure 3)              ####
#### Goals: Identify over-amplifiers &  Check uniformity of amplification        ####
#####################################################################################
#
#Use the GTscore perl script AmpliconReadCounter.pl to extract the data
# Download link - https://github.com/gjmckinney/GTscore 
#Reference:Garrett J. McKinney, Carita E. Pascal, William D. Templin, Sara E. Gilk-Baumer, 
# Tyler H. Dann, Lisa W. Seeb, and James E. Seeb. 2020. Dense SNP panels resolve closely 
# related Chinook salmon populations. Canadian Journal of Fisheries and Aquatic Sciences. 
# 77(3): 451-461. https://doi.org/10.1139/cjfas-2019-0067

#Packages required
library(dplyr)
library(ggrepel)
library(ggplot2)

#For the optimization plot we require the GTscore_locusSummary.txt from the perl script output

#Load the data associated with the figure in the manuscript use - 
load("./Rdata/all_rounds_locus_summary.RData") #six data.frame object will load, each indicating locusSummary for an optimization round.
# GTscore_locusSummary_Opt1
# GTscore_locusSummary_Opt2a
# GTscore_locusSummary_Opt2b
# GTscore_locusSummary_Opt3
# GTscore_locusSummary_Opt4
# GTscore_locusSummary_Opt5

#Set working directory
# setwd("working//directory")

# 1. Pick an optimization round to plot
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

opt_round <- "Opt #" #Change the opt round; will remain in the optimization plot

# Choose an optimization round to plot -
# GTscore_locusSummary_Opt1
# GTscore_locusSummary_Opt2a
# GTscore_locusSummary_Opt2b
# GTscore_locusSummary_Opt3
# GTscore_locusSummary_Opt4
# GTscore_locusSummary_Opt5

GTScore_locussummary <- GTscore_locusSummary_Opt1 #To streamline the plotting and reusing the code, only change the Opt# on this line

#If importing from text file of your extracted data
# GTScore_locussummary <- read.delim("GTScore_locussummary.txt",header=TRUE) 

head(GTScore_locussummary)
#Primer Read is the barcode. Primer probe read is our target sequences. 

#    2. Append and sort 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Total primer read count
Tot_primer_read_count <- sum(GTScore_locussummary$Primer.Reads)

#Total primer probe read count
Tot_primer_probe_read_count <- sum(GTScore_locussummary$Primer.Probe.Reads)

#Add a column with percent total primer read count in GTScore_locussummary_sorted
GTScore_locussummary$Percent.totalprimer.reads <- (GTScore_locussummary$Primer.Reads*100)/Tot_primer_read_count
GTScore_locussummary$Percent.totalprimerprobe.reads <- (GTScore_locussummary$Primer.Probe.Reads*100)/Tot_primer_probe_read_count
head(GTScore_locussummary)

#Sort the data with decreasing Percent.totalprimer.reads

GTScore_locussummary_sorted <- GTScore_locussummary %>% 
  arrange(desc(Percent.totalprimer.reads))
head(GTScore_locussummary_sorted)
tail(GTScore_locussummary_sorted)

#     3. Calculate relative and log relative read count
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

GTScore_locussummary_sorted$Relative_Primer_Read_Count = (GTScore_locussummary_sorted$Primer.Reads)/max((GTScore_locussummary_sorted$Primer.Reads))
GTScore_locussummary_sorted$Relative_Primer_Probe_Read_Count = (GTScore_locussummary_sorted$Primer.Probe.Reads)/max((GTScore_locussummary_sorted$Primer.Probe.Reads))

GTScore_locussummary_sorted$log_Relative_Primer_Read_Count = log10(GTScore_locussummary_sorted$Primer.Reads)/max(log10(GTScore_locussummary_sorted$Primer.Reads))
GTScore_locussummary_sorted$log_Relative_Primer_Probe_Read_Count = log10(GTScore_locussummary_sorted$Primer.Probe.Reads)/max(log10(GTScore_locussummary_sorted$Primer.Probe.Reads))

head(GTScore_locussummary_sorted)

sum(GTScore_locussummary_sorted$Primer.Reads)
sum(GTScore_locussummary_sorted$Primer.Probe.Reads)

# 4. Get the values
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#i. number of loci that account for top 50% of the cummulative sum of the total primer read count

GTScore_locussummary_sorted$cumSum_primerReadPercent <-  cumsum(GTScore_locussummary_sorted$Percent.totalprimer.reads)
#Identify the number of values to reach less than or equal to 50
top_loci_cs50percent_primers<- length(which(GTScore_locussummary_sorted$cumSum_primerReadPercent  <= 50))+1
top_loci_cs50percent_primers

#ii. number of loci that account for top 75% of the cummulative sum of the total primer read count

#Identify the number of values to reach less than or equal to 75
top_loci_cs_75percent_primers<- length(which(GTScore_locussummary_sorted$cumSum_primerReadPercent  <= 75))+1
top_loci_cs_75percent_primers

top_loci_cs_95percent_primers<- length(which(GTScore_locussummary_sorted$cumSum_primerReadPercent  <= 95))+1
top_loci_cs_95percent_primers

#iii. number of loci that account for top 50% of the cummulative sum of the total primer probe read count
GTScore_locussummary_sorted$cumSum_primerProbeReadPercent <-  cumsum(GTScore_locussummary_sorted$Percent.totalprimerprobe.reads)
#Identify the number of values to reach less than or equal to 50
top_loci_cs50percent_primersProbe<- length(which(GTScore_locussummary_sorted$cumSum_primerProbeReadPercent  <= 50))+1
top_loci_cs50percent_primersProbe

#iv.number of loci that account for top 75% of the cummulative sum of the total primer probe read count
#Identify the number of values to reach less than or equal to 75
top_loci_cs_75percent_primersProbe<- length(which(GTScore_locussummary_sorted$cumSum_primerProbeReadPercent  <= 75))+1
top_loci_cs_75percent_primersProbe

top_loci_cs_95percent_primersProbe<- length(which(GTScore_locussummary_sorted$cumSum_primerProbeReadPercent  <= 95))+1
top_loci_cs_95percent_primersProbe


Pseudo_96Variant <- which(GTScore_locussummary_sorted$Locus %in%  ("Pseudo_96Variant"))
PRNP_PRNPy <- which(GTScore_locussummary_sorted$Locus %in%  ("PRNP-PRNPy "))
rows_of_interest_indices <- c(Pseudo_96Variant, PRNP_PRNPy)

# Create a data frame for rows of interest
rows_of_interest <- GTScore_locussummary_sorted[rows_of_interest_indices, ]
# plot relative primer probe propotion to look at read distribution

#5. Plot
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Set figure dimensions (in inches) and resolution (DPI)
options(repr.plot.width = 10, repr.plot.height = 8, repr.plot.res = 300)

plot <- ggplot(GTScore_locussummary_sorted)+
  #primer counts
  #probe proportion data layers
  geom_point(size = 2, color = "#33cc00", alpha = 0.4,
             aes(x = 1:nrow(GTScore_locussummary_sorted), y = log_Relative_Primer_Probe_Read_Count))+
  #I add this layer again, with a smaller size to get the concentric points with different colors
  geom_point(size = 1, color = "#666666", alpha = 0.6,
             aes(x = 1:nrow(GTScore_locussummary_sorted), y = log_Relative_Primer_Probe_Read_Count))+
  
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "right",
        text = element_text(size = 25),
        axis.title = element_text(size = 15),  # Adjust axis label size
        axis.text.y = element_text(size = 10, face = "plain", color= "grey12"),
        axis.text.x = element_text(size = 10, face = "plain", color= "grey1"),
        panel.background = element_rect(fill = 'white', colour = 'grey13'))+
  scale_x_continuous(name = "Locus",
                     breaks = c(0,100,200,300,400,500,600),
                     labels = c(0,100,200,300,400,500,600),
                     limits = c(0,651),
                     expand = c(0,0))+
  scale_y_continuous(name = "Log Relative Read Count",
                     breaks = c(0,0.5,1), 
                     labels = c(0,0.5,1),
                     limits = c(0,1),
                     expand = c(0,0))+
  # labs(title = "Opt1")+
  #theme(plot.title = element_text(hjust = 0.5))+
  geom_smooth(se = F, color = "#050201", alpha = 0.2, linetype = "longdash", size=0.8,
              aes(x = 1:nrow(GTScore_locussummary_sorted), y = log_Relative_Primer_Probe_Read_Count))+
  geom_point(size = 0.5,
             aes(x = 1:nrow(GTScore_locussummary_sorted), y = log_Relative_Primer_Read_Count))+
  
  geom_hline(yintercept = 0)+          
  geom_hline(yintercept=0.5, linetype="dashed", color = "red")+
  
  
  # Highlight rows of interest
  geom_point(data = rows_of_interest, 
             aes(x = as.numeric(row.names(rows_of_interest)), 
                 y = log_Relative_Primer_Probe_Read_Count), 
             color = "blue", shape = 18, size = 8) +
  # Add labels for the locus names
  geom_text(data = rows_of_interest, 
            aes(x = as.numeric(row.names(rows_of_interest)), 
                y = log_Relative_Primer_Probe_Read_Count, 
                label = Locus), 
            color = "blue4", vjust = -3, size = 5, fontface = "bold")+
  #geom_label_repel is throwing error. FIX IT
  
  geom_vline(xintercept = top_loci_cs50percent_primers, linetype="dotdash",
             color = "red", size=0.7)+
  geom_label_repel(data = GTScore_locussummary_sorted[top_loci_cs50percent_primers,],
                   aes(x = c(top_loci_cs50percent_primers),
                       y = log_Relative_Primer_Read_Count,
                       label = top_loci_cs50percent_primers),
                   size =10,
                   nudge_y = 0,
                   direction = "x",
                   nudge_x = 50,
                   box.padding = unit(0.5, "lines"),
                   point.padding = unit(0.5, "lines"),
                   label.padding = unit(0.2, "lines"),
                   label.r = unit(0.1, "lines"),  # Radius of rounded corners
                   fill = "white",# Background color
                   alpha = 0.8,
                   color = "red4",  # Border color
                   family = "Arial",
                   fontface="bold"
  )+
  geom_vline(xintercept = top_loci_cs_75percent_primers, linetype="dotdash",
             color = "orange", size=0.7)+
  geom_label_repel(data = GTScore_locussummary_sorted[top_loci_cs_75percent_primers,],
                   aes(x = c(top_loci_cs_75percent_primers),
                       y = log_Relative_Primer_Read_Count,
                       label = top_loci_cs_75percent_primers),
                   size =10,
                   nudge_y = 0,
                   direction = "x",
                   nudge_x = 50,
                   box.padding = unit(0.5, "lines"),
                   point.padding = unit(0.5, "lines"),
                   label.padding = unit(0.2, "lines"),
                   label.r = unit(0.1, "lines"),  # Radius of rounded corners
                   fill = "white",# Background color
                   alpha = 0.8,
                   color = "orange3",  # Border color
                   family = "Arial",
                   fontface="bold"
  )+
  geom_vline(xintercept = top_loci_cs_95percent_primers, linetype="dotdash",
             color = "green", size=0.7)+
  geom_label_repel(data = GTScore_locussummary_sorted[top_loci_cs_95percent_primers,],
                   aes(x = c(top_loci_cs_95percent_primers),
                       y = log_Relative_Primer_Read_Count,
                       label = top_loci_cs_95percent_primers),
                   size =10,
                   nudge_y = 0,
                   direction = "x",
                   nudge_x = 50,
                   box.padding = unit(0.5, "lines"),
                   point.padding = unit(0.5, "lines"),
                   label.padding = unit(0.2, "lines"),
                   label.r = unit(0.1, "lines"),  # Radius of rounded corners
                   fill = "white",# Background color
                   alpha = 0.8,
                   color = "green4",  # Border color
                   family = "Arial",
                   fontface="bold")+  
     annotate(
    "label",
    x = Inf, y = Inf,
    label = opt_round,
    hjust = 1.1, vjust = 1.2,   # nudges inward from the top-right corner
    size = 6,
    fontface = "bold",
    color = "grey32"
  ) +
  theme(axis.title = element_blank(),  # Remove both x and y axis titles
        axis.text.x = element_text(size = 12, face = "bold", color = "grey12"),
        axis.text.y = element_text(size = 12, face = "bold", color = "grey12"),
        legend.position = "none")

plot

# To save the plot
# # Define the file path where the plot should be saved
# file_path <- file.path("figures", paste0(opt_round_file, "_OptimizationPlot.png"))
# 
# # Save the plot with the specified dimensions and 300 DPI
# ggsave(
#   filename = file_path,
#   plot = plot,       # The plot object to save
#   width = 2700/ 300, # Width in inches (900*3 pixels / 300 dpi) 
#   height = 1950/ 300, # Height in inches (650*3 pixels / 300 dpi) 
#   dpi = 300          # Resolution (dots per inch)
# )
# 
# file_path

#Repeat with other optimization rounds. 
