##########################################################################################
####          Figure 4 : Distribution of Locus and Sample genotyping summary          ####
#### Question: Do the summary statistics improve after each round of optimization?    ####
##########################################################################################


#Import data and required packages
load("./Rdata/Figure_4_Genotyping_summary.RData") #Locus and sample genotyping summary from all rounds 

library(readxl)
library(ggplot2)
library(reshape2)
library(scales)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Figure 4a
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
head(All_locus_summary)

# Convert PrimerProbeProportion to numeric if necessary
All_locus_summary$PrimerProbeProportion <- as.numeric(All_locus_summary$PrimerProbeProportion)

# Reshape the data to long format for easier plotting
long_summary <- melt(All_locus_summary, 
                     id.vars = c("OptRound", "Locus"), 
                     measure.vars = c("PrimerReads", "PrimerProbeReads", "PrimerProbeProportion"),
                     variable.name = "Metric", 
                     value.name = "Value")

# Create the plot
plot_locus <- ggplot(All_locus_summary, aes(x = OptRound, y = PrimerProbeProportion)) +
  geom_violin(trim = TRUE, alpha = 0.8, color = "black", fill = "grey55") +
  
  geom_jitter(width = 0.2, size = 2, alpha = 0.5, color = "grey49") +
  stat_summary(fun = mean, geom = "crossbar", width = 0.72, color = "grey41", size = 0.5) +
  labs(title = "Violin Plot of Primer Probe Proportions by OptRound",
       x = "Optimization Round",
       y = "Primer Probe Proportion") +
  theme_minimal(base_size = 15) +
  theme(axis.text.x = element_text(size = 12, face = "bold", color = "grey12"),
        axis.text.y = element_text(size = 12, face = "bold", color = "grey12"),
        legend.position = "none") +
  geom_violin(trim = TRUE, alpha = 0.4, color = "black", fill = "grey7") +
  scale_y_continuous(labels = comma)
plot_locus

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Figure 4b
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
head(All_individual_summary)
All_individual_summary <- as.data.frame(All_individual_summary)

# Convert Off_target_Reads to numeric if it's not already
All_individual_summary$PrimerProbeProportion <- as.numeric(All_individual_summary$PrimerProbeProportion)

# Create the boxplot with jittered points and formatted y-axis
plot_indiv <- ggplot(All_individual_summary, aes(x = OptRound, y = PrimerProbeProportion, fill = OptRound)) +
  geom_violin(trim = TRUE, alpha = 0.8, color = "black", fill = "grey55") +
  geom_jitter(width = 0.2, size = 2, alpha = 0.5, color = "grey6") +
  labs(title = "Comparison of PrimerProbeProportion Among Samples by OptRound",
       x = "Optimization Round",
       y = "Primer Probe Proportion") +
  theme_minimal(base_size = 15) +
  theme(axis.text.x = element_text(size = 12, face = "bold", color = "grey12"),
        axis.text.y = element_text(size = 12, face = "bold", color = "grey12"),
        legend.position = "none") +
  geom_violin(alpha = 0.4, color = "black", fill = "grey73") +
  stat_summary(fun = mean, geom = "crossbar", width = 0.72, color = "grey3", size = 0.5) +
  scale_y_continuous(labels = comma)

plot_indiv


