####################################################################################
####                             Figure S2                                      ####
####  Question: Do modifications help in making PRNP and Pseudogene work?       ####
####################################################################################
#Import data with modifications per well.
load("./Rdata/Figure_S2.RData")

library(ggplot2)
library(reshape2)
library(scales)

#Figure S2a
#Optimization 4 individual summary  

head(All_individual_summary)

All_individual_summary <- as.data.frame(All_individual_summary)

# Fit the ANOVA model
anova_model <- aov(PrimerProbeProportion ~ Data_group, data = All_individual_summary)
summary(anova_model)

# Convert Off_target_Reads to numeric if it's not already
All_individual_summary$PrimerProbeProportion <- as.numeric(All_individual_summary$PrimerProbeProportion)

# Create the boxplot with jittered points and formatted y-axis
plot <- ggplot(All_individual_summary, aes(x = Data_group, y = PrimerProbeProportion, fill = Data_group)) +
  geom_violin(trim = TRUE, alpha = 0.8, color = "black", fill = "grey55") +
  
  geom_jitter(width = 0.2, size = 2, alpha = 0.5, color = "grey6") +
  labs(title = "Comparison of PrimerProbeProportion Among Samples by OptRound4_groups",
       x = "Group",
       y = "Primer Probe Proportion") +
  theme_minimal(base_size = 15) +
  
  theme(axis.text.x = element_text(size = 10, face = "bold", color = "grey12", angle = 90),
        axis.text.y = element_text(size = 10, face = "bold", color = "grey12"),
        legend.position = "none") +
  geom_violin(alpha = 0.4, color = "black", fill = "grey73") +
  stat_summary(fun = mean, geom = "crossbar", width = 0.72, color = "grey3", size = 0.5) +
  scale_y_continuous(labels = comma)

plot


###################################################################
###################################################################


#Figure S2b
head(All_locus_summary)
# Fit the ANOVA model
anova_model <- aov(PrimerProbeProportion ~ Data_Group, data = All_locus_summary)
summary(anova_model)

Pseudo_96Variant <- which(All_locus_summary$Locus %in%  ("Pseudo_96Variant"))
PRNP_PRNPy <- which(All_locus_summary$Locus %in%  ("PRNP-PRNPy "))
rows_of_interest_indices <- c(Pseudo_96Variant, PRNP_PRNPy)

# Create the plot
plot <- ggplot(All_locus_summary, aes(x = Augmentation, y = PrimerProbeProportion)) +
  geom_violin(trim = TRUE, alpha = 0.2, color = "black", fill = "grey75") +
  geom_jitter(width = 0.2, size = 2, alpha = 0.5, color = "grey49") +
  stat_summary(fun = mean, geom = "crossbar", width = 0.72, color = "grey13", size = 0.5) +
  labs(title = "Violin Plots of group-wise Locus Primer probe summary",
       x = "Augmentation",
       y = "Primer Probe Proportion") +
  theme_minimal(base_size = 15) +
  theme(axis.text.x = element_text(size = 12, face = "bold", color = "grey12"),
        axis.text.y = element_text(size = 12, face = "bold", color = "grey12"),
        legend.position = "none") +
  geom_violin(trim = TRUE, alpha = 0.4, color = "black", fill = "grey7") +
  scale_y_continuous(labels = scales::comma) +
  
  # Add blue diamonds for the rows of interest
  geom_point(data = All_locus_summary[rows_of_interest_indices, ], 
             aes(x = Augmentation, y = PrimerProbeProportion), 
             color = "white", fill = "white",  shape = 23, size = 4.5)+
  
  # Add blue diamonds for the rows of interest
  geom_point(data = All_locus_summary[rows_of_interest_indices, ], 
             aes(x = Augmentation, y = PrimerProbeProportion), 
             color = "blue", shape = 18, size = 4)

# Display the plot
print(plot)
