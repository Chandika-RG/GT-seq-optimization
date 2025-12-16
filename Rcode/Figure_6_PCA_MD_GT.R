############################################################################
####                     Figure 6: PCA                                  ####
####     Question: Do the OVSNP60 genotypes match GT-seq genotypes?     ####
############################################################################
library(readxl)
library(adegenet)
library(dplyr) 
library(ggplot2)
library(ggrepel)
library(dplyr)

#Load data
load("./Rdata/Figure_6_PCA_geno.RData")
head(genepop_data)

########################################################
########################################################
#PCA
sum(is.na(genepop_data$tab))#number of NAs

#replace NAs
Data <- scaleGen(genepop_data, NA.method = "mean")
dim(Data)

Data[1:5, 1:5]

pca <- dudi.pca(Data, cent=FALSE,scale=FALSE,scannf=FALSE,nf=15)

#ScreePlot
barplot(pca$eig[1:100],main="PCA eigenvalues", col=heat.colors(100))

pca_data <- as.data.frame(pca$li)
pca_data <- cbind(SampleID = rownames(pca_data), pca_data)
head(pca_data)

pca_data$PanelType <- sapply(strsplit(pca_data$SampleID, "_"), `[`, 2)
pca_data$Sample_num <- sapply(strsplit(pca_data$SampleID, "_"), `[`, 1)
head(pca_data)

# Contribution of each PC
contributions <- pca$eig / sum(pca$eig) * 100
sum(contributions)
# Display the contributions in a table
contrib_table <- data.frame(PC = 1:length(contributions), Contribution = contributions)
print(contrib_table)

col = c("#0B5345", "#0E6251", "#117864", "#117A65", "#145A32", "#154360",
        "#196F3D", "#1A5276", "#1B4F72", "#1F618D", "#239B56", "#283747",
        "#2874A6", "#2E4053", "#4A235A", "#512E5F", "#512E5F", "#512E5F",
        "#566573", "#5B2C6F", "#6C3483", "#6C3483", "#641E16", "#6E2C00",
        "#76448A", "#73281A", "#7B241C", "#7B441C", "#7D3C98", "#7D6608",
        "#7E5109", "#922B21", "#943126", "#9C640C", "#9C640C", "#A93226",
        "#B03A2E", "#B9770E", "#C0392B", "#C0392B", "#D35400", "#D35400",
        "#3D2C8D", "#FF5733", "#008080", "#8B0000") 



shape = c(21,24, 22)

plot <- ggplot(
  pca_data,
  aes(x = Axis1, y = Axis2, color = Sample_num, shape = PanelType, fill = Sample_num)
) +
  geom_point(size = 5, stroke = 0.4, color = "grey18", alpha = 0.8) +
  
  scale_shape_manual(values = c(21, 24, 22)) +
  scale_color_manual(values = col) +
  scale_fill_manual(values = col) +
    labs(
    x = paste("PC1 (", round(contrib_table[1, 2], 2), "%)"),
    y = paste("PC2 (", round(contrib_table[2, 2], 2), "%)")) +
    theme(
    axis.text = element_text(size = 8, face = "bold", color = "grey12"),
    axis.title = element_text(size = 8, face = "bold", color = "grey12"),
    legend.position = "none") +
    geom_label_repel(
    data = pca_data,
    aes(
      label = SampleID,
      fill = after_scale(alpha(fill, 0.3))),
    size = 4,
    nudge_y = 0.35,
    nudge_x = 0.35,
    direction = "both",
    box.padding = unit(0.25, "lines"),
    point.padding = unit(0.25, "lines"),
    label.padding = unit(0.15, "lines"),
    label.r = unit(0.1, "lines"),
    family = "Arial",
    fontface = "bold",
    max.overlaps = 50,
    segment.curvature = 0.3,
    segment.ncp = 3,
    segment.angle = 45,
    label.box = "none") +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "grey80", size = 0.5),
    panel.grid.minor = element_line(color = "grey90", size = 0.3),
    axis.line = element_line(color = "black")
  )

plot 




