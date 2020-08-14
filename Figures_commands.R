if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
devtools::install_github("jbisanz/qiime2R") # current version was 0.99.20

library(tidyverse)
library(qiime2R)
library(cowplot)
library(scales)

metadata=read_q2metadata("../qiime2_2020_2_analysis/metadata/DPSW_qiime2_R_metadata.txt")
SVs = read_qza("../qiime2_2020_2_analysis/feature_tables/table-gg-99-no-chlo-mito-minfreq20-no-b-d-S25.qza")
taxonomy = read_qza("../qiime2_2020_2_analysis/taxonomy/taxonomy-gg-99.qza")
tree = read_qza("../qiime2_2020_2_analysis/trees/fragment_insertion_out/tree.qza")

shannon = read_qza("../qiime2_2020_2_analysis/core-metrics-results-15161/alpha_diversity/shannon_vector.qza")
shannon = shannon$data %>% rownames_to_column("SampleID") # this moves the sample names to a new column that matches the metadata and allows them to be merged

faith = read_qza("../qiime2_2020_2_analysis/core-metrics-results-15161/alpha_diversity/faith_pd_vector.qza")
faith = faith$data %>% rownames_to_column("SampleID") # this moves the sample names to a new column that matches the metadata and allows them to be merged

evenness = read_qza("../qiime2_2020_2_analysis/core-metrics-results-15161/alpha_diversity/evenness_vector.qza")
evenness = evenness$data %>% rownames_to_column("SampleID") # this moves the sample names to a new column that matches the metadata and allows them to be merged

observed_otus = read_qza("../qiime2_2020_2_analysis/core-metrics-results-15161/alpha_diversity/observed_otus_vector.qza")
observed_otus = observed_otus$data %>% rownames_to_column("SampleID") # this moves the sample names to a new column that matches the metadata and allows them to be merged

unweighted = read_qza("../qiime2_2020_2_analysis/core-metrics-results-15161/beta_diversity/unweighted_unifrac_pcoa_results.qza")
weighted = read_qza("../qiime2_2020_2_analysis/core-metrics-results-15161/beta_diversity/weighted_unifrac_pcoa_results.qza")

# Shows you how many samples in each category
# 4 metadata samples do not have Shannon values
gplots::venn(list(metadata=metadata$SampleID, shannon=shannon$SampleID))
gplots::venn(list(metadata=metadata$SampleID, observed_otus=observed_otus$SampleID))

# Add alpha diversity values to metadata
metadata = metadata %>% left_join(shannon)
metadata = metadata %>% left_join(faith)
metadata = metadata %>% left_join(evenness)
metadata = metadata %>% left_join(observed_otus)
head(metadata)

# Alpha diveristy
# shannon
metadata %>%
  filter(!is.na(shannon)) %>%
  ggplot(aes(x=`pmi_months`, y=shannon, color=`sample_detail`)) +
  ylim(2,12) +
  stat_summary(geom="errorbar", fun.data=mean_se, width=0, orientation="x") +
  stat_summary(geom="line", fun.data=mean_se, orientation="x", size=1) +
  stat_summary(geom="point", fun.data=mean_se, orientation="x", size=1) +
  labs(x = "Time (Months)", y = "Shannon's Diversity Index", title = "Shannon's Diversity over Decomposition") +
  theme_q2r() + 
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=11),
        plot.title = element_text(size=12, hjust = 0.5),
        legend.text = element_text(size=7),
        legend.title = element_text(size=7),
        legend.position = c(0.7,0.04),
        legend.direction = "horizontal",
        legend.key.height = unit(0.01, "cm"),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill="white",color="black")) +
  scale_colour_manual(name='Location:', values=c('1 Meter North'='mediumorchid4', 'Underneath Carcass'='royalblue3', 'Soil Control'='mediumseagreen')) 
ggsave("shannon_plot.png",dpi=300, width=7, height=5, units="in")

# faith
metadata %>%
  filter(!is.na(faith_pd)) %>%
  ggplot(aes(x=`pmi_months`, y=faith_pd, color=`sample_detail`)) +
  stat_summary(geom="errorbar", fun.data=mean_se, width=0, orientation="x") +
  stat_summary(geom="line", fun.data=mean_se, orientation="x", size=1) +
  stat_summary(geom="point", fun.data=mean_se, orientation="x", size=1) +
  labs(x = "Time (Months)", y = "Faith's Phylogenetic Diversity", title = "Faith's Phylogenetic Diversity over Decomposition") +
  theme_q2r() + 
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=11),
        plot.title = element_text(size=12, hjust = 0.5),
        legend.text = element_text(size=7),
        legend.title = element_text(size=7),
        legend.position = c(0.7,0.04),
        legend.direction = "horizontal",
        legend.key.height = unit(0.01, "cm"),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill="white",color="black")) +
  scale_colour_manual(name='Location:', values=c('1 Meter North'='mediumorchid4', 'Underneath Carcass'='royalblue3', 'Soil Control'='mediumseagreen')) 
ggsave("faith_plot.png",dpi=300, width=7, height=5, units="in")

# evenness
metadata %>%
  filter(!is.na(pielou_e)) %>%
  ggplot(aes(x=`pmi_months`, y=pielou_e, color=`sample_detail`)) +
  stat_summary(geom="errorbar", fun.data=mean_se, width=0, orientation="x") +
  stat_summary(geom="line", fun.data=mean_se, orientation="x", size=1) +
  stat_summary(geom="point", fun.data=mean_se, orientation="x", size=1) +
  labs(x = "Time (Months)", y = "Pielou's Evenness Index", title = "Pielou's Evenness Index over Decomposition") +
  theme_q2r() + 
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=11),
        plot.title = element_text(size=12, hjust = 0.5),
        legend.text = element_text(size=7),
        legend.title = element_text(size=7),
        legend.position = c(0.7,0.04),
        legend.direction = "horizontal",
        legend.key.height = unit(0.01, "cm"),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill="white",color="black")) +
  scale_colour_manual(name='Location:', values=c('1 Meter North'='mediumorchid4', 'Underneath Carcass'='royalblue3', 'Soil Control'='mediumseagreen')) 
ggsave("evenness_plot.png",dpi=300, width=7, height=5, units="in")

# observed_otus
metadata %>%
  filter(!is.na(observed_otus)) %>%
  ggplot(aes(x=`pmi_months`, y=observed_otus, color=`sample_detail`)) +
  stat_summary(geom="errorbar", fun.data=mean_se, width=0, orientation="x") +
  stat_summary(geom="line", fun.data=mean_se, orientation="x", size=1) +
  stat_summary(geom="point", fun.data=mean_se, orientation="x", size=1) +
  labs(x = "Time (Months)", y = "ASV Richness", title = "ASV Richness over Decomposition") +
  theme_q2r() + 
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=11),
        plot.title = element_text(size=12, hjust = 0.5),
        legend.text = element_text(size=7),
        legend.title = element_text(size=7),
        legend.position = c(0.7,0.04),
        legend.direction = "horizontal",
        legend.key.height = unit(0.01, "cm"),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill="white",color="black")) +
  scale_colour_manual(name='Location:', values=c('1 Meter North'='mediumorchid4', 'Underneath Carcass'='royalblue3', 'Soil Control'='mediumseagreen')) 
ggsave("observed_otus_plot.png",dpi=300, width=7, height=5, units="in")

############### MULTIPLOT
# Alpha diveristy

p1=metadata %>%
  filter(!is.na(shannon)) %>%
  ggplot(aes(x=`pmi_months`, y=shannon, color=`sample_detail`)) +
  ylim(2,12) +
  stat_summary(geom="errorbar", fun.data=mean_se, width=0, orientation="x") +
  stat_summary(geom="line", fun.data=mean_se, orientation="x", size=1) +
  stat_summary(geom="point", fun.data=mean_se, orientation="x", size=1) +
  labs(x = "Time (Months)", y = "Shannon's Diversity Index", title = "Shannon's Diversity") +
  guides(color = guide_legend(override.aes = list(size=2))) +
  theme_q2r() + 
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10),
        plot.title = element_text(size=12, hjust = 0.5, face = "bold"),
        strip.text = element_text(size=6.5),
        legend.text = element_text(size=10),
        legend.title = element_text(size=11, hjust = 0.5, face = "bold"),
        legend.position = c(0.858, 0.11),
        legend.direction = "vertical",
        legend.key.height = unit(0.01, "cm"),
        legend.key.width = unit(0.01, "cm"),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill="white",color="black")) +
  scale_colour_manual(name='Location', values=c('1 Meter North'='mediumorchid4', 'Underneath Carcass'='royalblue3', 'Soil Control'='mediumseagreen')) +
  geom_hline(yintercept=mean(na.omit(metadata$shannon))) +
  geom_hline(yintercept=(mean(na.omit(metadata$shannon))+2*(sd(na.omit(metadata$shannon)))), linetype="dashed") +
  geom_hline(yintercept=(mean(na.omit(metadata$shannon))-2*(sd(na.omit(metadata$shannon)))), linetype="dashed")
p1

p2=metadata %>%
  filter(!is.na(faith_pd)) %>%
  ggplot(aes(x=`pmi_months`, y=faith_pd, color=`sample_detail`)) +
  stat_summary(geom="errorbar", fun.data=mean_se, width=0, orientation="x") +
  stat_summary(geom="line", fun.data=mean_se, orientation="x", size=1) +
  stat_summary(geom="point", fun.data=mean_se, orientation="x", size=1) +
  labs(x = "Time (Months)", y = "Faith's Phylogenetic Diversity", title = "Faith's Phylogenetic Diversity") +
  theme_q2r() + 
  theme(axis.text=element_text(size=8),
        axis.title.x = element_blank(),
        axis.title=element_text(size=10),
        plot.title = element_text(size=12, hjust = 0.5, face = "bold"),
        strip.text = element_text(size=6.5),
        legend.text = element_text(size=7),
        legend.title = element_text(size=7),
        legend.position = "none",
        legend.direction = "horizontal",
        legend.key.height = unit(0.01, "cm"),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill="white",color="black")) +
  scale_colour_manual(name='Location:', values=c('1 Meter North'='mediumorchid4', 'Underneath Carcass'='royalblue3', 'Soil Control'='mediumseagreen')) +
  geom_hline(yintercept=mean(na.omit(metadata$faith_pd))) +
  geom_hline(yintercept=(mean(na.omit(metadata$faith_pd))+2*(sd(na.omit(metadata$faith_pd)))), linetype="dashed") +
  geom_hline(yintercept=(mean(na.omit(metadata$faith_pd))-2*(sd(na.omit(metadata$faith_pd)))), linetype="dashed")
p2

p3=metadata %>%
  filter(!is.na(pielou_e)) %>%
  ggplot(aes(x=`pmi_months`, y=pielou_e, color=`sample_detail`)) +
  stat_summary(geom="errorbar", fun.data=mean_se, width=0, orientation="x") +
  stat_summary(geom="line", fun.data=mean_se, orientation="x", size=1) +
  stat_summary(geom="point", fun.data=mean_se, orientation="x", size=1) +
  labs(x = "Time (Months)", y = "Pielou's Evenness Index", title = "Pielou's Evenness") +
  theme_q2r() + 
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10),
        plot.title = element_text(size=12, hjust = 0.5, face = "bold"),
        strip.text = element_text(size=6.5),
        legend.text = element_text(size=7),
        legend.title = element_text(size=7),
        legend.position = "none",
        legend.direction = "horizontal",
        legend.key.height = unit(0.01, "cm"),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill="white",color="black")) +
  scale_colour_manual(name='Location:', values=c('1 Meter North'='mediumorchid4', 'Underneath Carcass'='royalblue3', 'Soil Control'='mediumseagreen')) +
  geom_hline(yintercept=mean(na.omit(metadata$pielou_e))) +
  geom_hline(yintercept=(mean(na.omit(metadata$pielou_e))+2*(sd(na.omit(metadata$pielou_e)))), linetype="dashed") +
  geom_hline(yintercept=(mean(na.omit(metadata$pielou_e))-2*(sd(na.omit(metadata$pielou_e)))), linetype="dashed")
p3

p4=metadata %>%
  filter(!is.na(observed_otus)) %>%
  ggplot(aes(x=`pmi_months`, y=observed_otus, color=`sample_detail`)) +
  stat_summary(geom="errorbar", fun.data=mean_se, width=0, orientation="x") +
  stat_summary(geom="line", fun.data=mean_se, orientation="x", size=1) +
  stat_summary(geom="point", fun.data=mean_se, orientation="x", size=1) +
  labs(x = "Time (Months)", y = "ASV Richness", title = "ASV Richness") +
  theme_q2r() + 
  theme(axis.text=element_text(size=8),
        axis.title.x = element_blank(),
        axis.title=element_text(size=10),
        plot.title = element_text(size=12, hjust = 0.5, face = "bold"),
        strip.text = element_text(size=6.5),
        legend.text = element_text(size=7),
        legend.title = element_text(size=7),
        legend.position = "none",
        legend.direction = "horizontal",
        legend.key.height = unit(0.01, "cm"),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill="white",color="black")) +
  scale_colour_manual(name='Location:', values=c('1 Meter North'='mediumorchid4', 'Underneath Carcass'='royalblue3', 'Soil Control'='mediumseagreen'))  +
  geom_hline(yintercept=mean(na.omit(metadata$observed_otus))) +
  geom_hline(yintercept=(mean(na.omit(metadata$observed_otus))+2*(sd(na.omit(metadata$observed_otus)))), linetype="dashed") +
  geom_hline(yintercept=(mean(na.omit(metadata$observed_otus))-2*(sd(na.omit(metadata$observed_otus)))), linetype="dashed")
p4

#Multiplot
p5 = plot_grid(p4, p2, p3, p1, labels=c("a.", "b.", "c.", "d."), ncol = 2, nrow = 2)
p5
save_plot("alpha_multiplot.png", p5, base_height = 8, base_width = 12,dpi=300)

# Beta diversity
shannon_data = read_qza("../qiime2_2020_2_analysis/core-metrics-results-15161/alpha_diversity/shannon_vector.qza")$data %>% rownames_to_column("SampleID") 
# Unweighted
unweighted$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  left_join(metadata) %>%
  left_join(shannon_data) %>%
  ggplot(aes(x=PC1, y=PC2, fill=sample_detail, size=pmi_months)) +
  geom_point(colour="black",pch=21) +
  stat_ellipse(type="t", aes(color = sample_detail)) +
  guides(size = guide_legend(override.aes = list(linetype = 0)),
         fill = guide_legend(override.aes = list(linetype = 0,size=3)),
         color = FALSE) +
  labs(fill = "Location",
       size = "Time (Months)",
       x = paste("PC1: ", round(100*unweighted$data$ProportionExplained[1]), "%"), 
       y = paste("PC2: ", round(100*unweighted$data$ProportionExplained[2]), "%"), 
       title = "Unweighted UniFrac PCoA") +
  theme_bw() +
  theme(axis.text=element_text(size=11),
        axis.title=element_text(size=12),
        plot.title = element_text(size=14, hjust = 0.5),
        legend.position = c(0.855,0.105),
        legend.box = "horizontal",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(color = "black"),
        legend.background = element_rect(fill="white",color="black")) +
  scale_fill_manual(name='Location:', values=c('1 Meter North'='mediumorchid4', 'Underneath Carcass'='royalblue3', 'Soil Control'='mediumseagreen')) +
  scale_colour_manual(name='Location:', values=c('1 Meter North'='mediumorchid4', 'Underneath Carcass'='royalblue3', 'Soil Control'='mediumseagreen'))
ggsave("unweighted_unifrac_plot.png",dpi=300, width=6, height=6, units="in")


# Beta diversity
# Weighted
weighted$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  left_join(metadata) %>%
  left_join(shannon_data) %>%
  ggplot(aes(x=PC1, y=PC2, fill=sample_detail, size=pmi_months)) +
  geom_point(colour="black",pch=21) +
  stat_ellipse(type="t", aes(color = sample_detail)) +
  guides(size = guide_legend(override.aes = list(linetype = 0)),
         fill = guide_legend(override.aes = list(linetype = 0,size=3)),
         color = FALSE) +
  labs(fill = "Location",
       size = "Time (Months)",
       x = paste("PC1: ", round(100*weighted$data$ProportionExplained[1]), "%"), 
       y = paste("PC2: ", round(100*weighted$data$ProportionExplained[2]), "%"), 
       title = "Weighted UniFrac PCoA") +
  theme_bw() +
  theme(axis.text=element_text(size=11),
        axis.title=element_text(size=12),
        plot.title = element_text(size=14, hjust = 0.5),
        legend.position = c(0.855,0.105),
        legend.box = "horizontal",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(color = "black"),
        legend.background = element_rect(fill="white",color="black")) +
  scale_fill_manual(name='Location:', values=c('1 Meter North'='mediumorchid4', 'Underneath Carcass'='royalblue3', 'Soil Control'='mediumseagreen')) +
  scale_colour_manual(name='Location:', values=c('1 Meter North'='mediumorchid4', 'Underneath Carcass'='royalblue3', 'Soil Control'='mediumseagreen'))
ggsave("weighted_unifrac_plot.png",dpi=300, width=6, height=6, units="in")

############### MULTIPLOT
# Beta diversity
# Unweighted
p6=unweighted$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  left_join(metadata) %>%
  left_join(shannon_data) %>%
  ggplot(aes(x=PC1, y=PC2, fill=sample_detail, size=pmi_months)) +
  geom_point(colour="black",pch=21) +
  stat_ellipse(type="t", aes(color = sample_detail)) +
  guides(size = FALSE,
         fill = FALSE,
         color = FALSE) +
  labs(fill = "Location",
       size = "Time (Months)",
       x = paste("PC1: ", round(100*unweighted$data$ProportionExplained[1]), "%"), 
       y = paste("PC2: ", round(100*unweighted$data$ProportionExplained[2]), "%"), 
       title = "Taxonomic: Unweighted UniFrac") +
  theme_bw() +
  theme(axis.text=element_text(size=6.5),
        axis.title=element_text(size=8),
        plot.title = element_text(size=9, hjust = 0.5, face = "bold"),
        legend.position = c(0.855,0.105),
        legend.box = "horizontal",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(color = "black"),
        legend.background = element_rect(fill="white",color="black"),
        legend.key.height = unit(0.01, "cm"),
        legend.key.width = unit(0.01, "cm")) +
  scale_fill_manual(name='Location', values=c('1 Meter North'='mediumorchid4', 'Underneath Carcass'='royalblue3', 'Soil Control'='mediumseagreen')) +
  scale_colour_manual(name='Location', values=c('1 Meter North'='mediumorchid4', 'Underneath Carcass'='royalblue3', 'Soil Control'='mediumseagreen'))
p6

# Beta diversity
# Weighted
p7=weighted$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  left_join(metadata) %>%
  left_join(shannon_data) %>%
  ggplot(aes(x=PC1, y=PC2, fill=sample_detail, size=pmi_months)) +
  geom_point(colour="black",pch=21) +
  stat_ellipse(type="t", aes(color = sample_detail)) +
  guides(size = guide_legend(override.aes = list(linetype = 0)),
         fill = guide_legend(override.aes = list(linetype = 0,size=2)),
         color = FALSE) +
  labs(fill = "Location",
       size = "Time (Months)",
       x = paste("PC1: ", round(100*weighted$data$ProportionExplained[1]), "%"), 
       y = paste("PC2: ", round(100*weighted$data$ProportionExplained[2]), "%"), 
       title = "Taxonomic: Weighted UniFrac") +
  scale_size_continuous(breaks=c(0,0.5,1,6,12,24,72,84,120),
                        limits=c(0,120)) +
  theme_bw() +
  theme(axis.text=element_text(size=6.5),
        axis.title=element_text(size=8),
        plot.title = element_text(size=9, hjust = 0.5, face = "bold"),
        legend.text = element_text(size=7),
        legend.title = element_text(size=8, face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(color = "black"),
        legend.background = element_rect(fill="white",color="white"),
        legend.key.height = unit(0.4, "cm"),
        legend.key.width = unit(0.01, "cm")) +
  scale_fill_manual(name='Location', values=c('1 Meter North'='mediumorchid4', 'Underneath Carcass'='royalblue3', 'Soil Control'='mediumseagreen')) +
  scale_colour_manual(name='Location', values=c('1 Meter North'='mediumorchid4', 'Underneath Carcass'='royalblue3', 'Soil Control'='mediumseagreen'))
p7



## Beta Diversity Longitudinal 
# Unweighted UniFrac
pcoa_vol_unwunifrac = read.table("../qiime2_2020_2_analysis/longitudinal/all/pcoa-vol-unwunifrac.tsv", header = TRUE, sep = "\t", row.names = 1)

p16 = ggplot(pcoa_vol_unwunifrac, aes(x=pmi_months, y=-(Axis1), color=sample_detail)) +
  scale_y_continuous(breaks=pretty_breaks()) +
  stat_summary(geom="errorbar", fun.data=mean_se, width=0, orientation="x") +
  stat_summary(geom="line", fun.data=mean_se, orientation="x", size=1) +
  stat_summary(geom="point", fun.data=mean_se, orientation="x", size=1) +
  labs(x = "Time (Months)", y = "Axis 1 Distance", title = "Unweighted UniFrac Distance over Decomposition") +
  theme_q2r() + 
  theme(axis.text=element_text(size=8),
        axis.title.x = element_blank(),
        axis.title=element_text(size=9),
        plot.title = element_text(size=12, hjust = 0.5, face = "bold"),
        strip.text = element_text(size=6.5),
        legend.text = element_text(size=7),
        legend.title = element_text(size=8, hjust = 0.5, face = "bold"),
        legend.position = c(0.874, 0.11),
        legend.direction = "vertical",
        legend.key.height = unit(0.01, "cm"),
        legend.key.width = unit(0.01, "cm"),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill="white",color="black")) +
  scale_colour_manual(name='Location', values=c('1 Meter North'='mediumorchid4', 'Underneath Carcass'='royalblue3', 'Soil Control'='mediumseagreen')) + 
  geom_hline(yintercept=mean(pcoa_vol_unwunifrac$Axis1)) +
  geom_hline(yintercept=mean(pcoa_vol_unwunifrac$Axis1)+2*(sd(pcoa_vol_unwunifrac$Axis1)), linetype="dashed") +
  geom_hline(yintercept=mean(pcoa_vol_unwunifrac$Axis1)-2*(sd(pcoa_vol_unwunifrac$Axis1)), linetype="dashed")
p16

# Weighted UniFrac
pcoa_vol_wunifrac = read.table("../qiime2_2020_2_analysis/longitudinal/all/pcoa-vol-wunifrac.tsv", header = TRUE, sep = "\t", row.names = 1)

p17 = ggplot(pcoa_vol_wunifrac, aes(x=pmi_months, y=Axis1, color=sample_detail)) +
  scale_y_continuous(breaks=pretty_breaks()) +
  stat_summary(geom="errorbar", fun.data=mean_se, width=0, orientation="x") +
  stat_summary(geom="line", fun.data=mean_se, orientation="x", size=1) +
  stat_summary(geom="point", fun.data=mean_se, orientation="x", size=1) +
  labs(x = "Time (Months)", y = "Axis 1 Distance", title = "Weighted UniFrac Distance over Decomposition") +
  theme_q2r() + 
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=9),
        plot.title = element_text(size=12, hjust = 0.5, face = "bold"),
        strip.text = element_text(size=6.5),
        legend.text = element_text(size=7),
        legend.title = element_text(size=8, hjust = 0.5, face = "bold"),
        legend.position = c(0.905, 0.16),
        legend.direction = "vertical",
        legend.key.height = unit(0.01, "cm"),
        legend.key.width = unit(0.01, "cm"),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill="white",color="black")) +
  scale_colour_manual(name='Location', values=c('1 Meter North'='mediumorchid4', 'Underneath Carcass'='royalblue3', 'Soil Control'='mediumseagreen')) + 
  geom_hline(yintercept=mean(pcoa_vol_wunifrac$Axis1)) +
  geom_hline(yintercept=mean(pcoa_vol_wunifrac$Axis1)+2*(sd(pcoa_vol_wunifrac$Axis1)), linetype="dashed") +
  geom_hline(yintercept=mean(pcoa_vol_wunifrac$Axis1)-2*(sd(pcoa_vol_wunifrac$Axis1)), linetype="dashed")
p17

#Multiplot
p18 = plot_grid(p16+theme(legend.position="none"),p17, labels=c("a.", "b."), ncol = 1, nrow = 2)
p18
save_plot("beta_long_multiplot.png", p18, base_height = 5.5, base_width = 7,dpi=300)

## Important Genera Longitudinal 
genera_vol = read.table("../qiime2_2020_2_analysis/longitudinal/all/genera_lme/volatility-data.tsv", header = TRUE, sep = "\t", row.names = 1)

# k__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhodospirillales;f__Rhodospirillaceae;g__
my_title <- expression(paste(italic("Rhodospirillaceae genus")))

p19 = ggplot(genera_vol, aes(x=pmi_months, y=Rhodospirillaceae_genus, color=sample_detail)) +
  scale_y_continuous(breaks=pretty_breaks()) +
  stat_summary(geom="errorbar", fun.data=mean_se, width=0, orientation="x") +
  stat_summary(geom="line", fun.data=mean_se, orientation="x", size=1) +
  stat_summary(geom="point", fun.data=mean_se, orientation="x", size=1) +
  labs(x = "Time (Months)", y = "Relative Abundance", title = my_title) +
  theme_q2r() + 
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=9),
        plot.title = element_text(size=11, hjust = 0.5, face = "bold"),
        strip.text = element_text(size=6.5),
        legend.text = element_text(size=7),
        legend.title = element_text(size=8, hjust = 0.5, face = "bold"),
        legend.position = c(0.874, 0.11),
        legend.direction = "vertical",
        legend.key.height = unit(0.01, "cm"),
        legend.key.width = unit(0.01, "cm"),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill="white",color="black")) +
  scale_colour_manual(name='Location', values=c('1 Meter North'='mediumorchid4', 'Underneath Carcass'='royalblue3', 'Soil Control'='mediumseagreen')) + 
  geom_hline(yintercept=mean(genera_vol$Rhodospirillaceae_genus)) +
  geom_hline(yintercept=mean(genera_vol$Rhodospirillaceae_genus)+2*(sd(genera_vol$Rhodospirillaceae_genus)), linetype="dashed") +
  geom_hline(yintercept=mean(genera_vol$Rhodospirillaceae_genus)-2*(sd(genera_vol$Rhodospirillaceae_genus)), linetype="dashed")
p19

# k__Bacteria;p__Actinobacteria;c__Acidimicrobiia;o__Acidimicrobiales;f__EB1017;g__
my_title <- expression(paste(italic("EB1017 genus")))

p20 = ggplot(genera_vol, aes(x=pmi_months, y=EB1017_genus, color=sample_detail)) +
  scale_y_continuous(breaks=pretty_breaks()) +
  stat_summary(geom="errorbar", fun.data=mean_se, width=0, orientation="x") +
  stat_summary(geom="line", fun.data=mean_se, orientation="x", size=1) +
  stat_summary(geom="point", fun.data=mean_se, orientation="x", size=1) +
  labs(x = "", y = "Relative Abundance", title = my_title) +
  theme_q2r() + 
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=9),
        plot.title = element_text(size=11, hjust = 0.5, face = "bold"),
        strip.text = element_text(size=6.5),
        legend.text = element_text(size=7),
        legend.title = element_text(size=8, hjust = 0.5, face = "bold"),
        legend.position = c(0.874, 0.11),
        legend.direction = "vertical",
        legend.key.height = unit(0.01, "cm"),
        legend.key.width = unit(0.01, "cm"),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill="white",color="black")) +
  scale_colour_manual(name='Location', values=c('1 Meter North'='mediumorchid4', 'Underneath Carcass'='royalblue3', 'Soil Control'='mediumseagreen')) + 
  geom_hline(yintercept=mean(genera_vol$EB1017_genus)) +
  geom_hline(yintercept=mean(genera_vol$EB1017_genus)+2*(sd(genera_vol$EB1017_genus)), linetype="dashed") +
  geom_hline(yintercept=mean(genera_vol$EB1017_genus)-2*(sd(genera_vol$EB1017_genus)), linetype="dashed")
p20

# k__Bacteria;p__Verrucomicrobia;c__[Spartobacteria];o__[Chthoniobacterales];f__[Chthoniobacteraceae];g__Chthoniobacter
my_title <- expression(paste(italic("Chthoniobacter")))

p21 = ggplot(genera_vol, aes(x=pmi_months, y=Chthoniobacter, color=sample_detail)) +
  scale_y_continuous(breaks=pretty_breaks()) +
  stat_summary(geom="errorbar", fun.data=mean_se, width=0, orientation="x") +
  stat_summary(geom="line", fun.data=mean_se, orientation="x", size=1) +
  stat_summary(geom="point", fun.data=mean_se, orientation="x", size=1) +
  labs(x = "", y = "", title = my_title) +
  theme_q2r() + 
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=9),
        plot.title = element_text(size=11, hjust = 0.5, face = "bold"),
        strip.text = element_text(size=6.5),
        legend.text = element_text(size=7),
        legend.title = element_text(size=8, hjust = 0.5, face = "bold"),
        legend.position = c(0.874, 0.11),
        legend.direction = "vertical",
        legend.key.height = unit(0.01, "cm"),
        legend.key.width = unit(0.01, "cm"),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill="white",color="black")) +
  scale_colour_manual(name='Location', values=c('1 Meter North'='mediumorchid4', 'Underneath Carcass'='royalblue3', 'Soil Control'='mediumseagreen')) + 
  geom_hline(yintercept=mean(genera_vol$Chthoniobacter)) +
  geom_hline(yintercept=mean(genera_vol$Chthoniobacter)+2*(sd(genera_vol$Chthoniobacter)), linetype="dashed") +
  geom_hline(yintercept=mean(genera_vol$Chthoniobacter)-2*(sd(genera_vol$Chthoniobacter)), linetype="dashed")
p21

# k__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhizobiales;f__Hyphomicrobiaceae;g__Devosia
my_title <- expression(paste(italic("Devosia")))

p22 = ggplot(genera_vol, aes(x=pmi_months, y=Devosia, color=sample_detail)) +
  scale_y_continuous(breaks=pretty_breaks()) +
  stat_summary(geom="errorbar", fun.data=mean_se, width=0, orientation="x") +
  stat_summary(geom="line", fun.data=mean_se, orientation="x", size=1) +
  stat_summary(geom="point", fun.data=mean_se, orientation="x", size=1) +
  labs(x = "Time (Months)", y = "", title = my_title) +
  theme_q2r() + 
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=9),
        plot.title = element_text(size=11, hjust = 0.5, face = "bold"),
        strip.text = element_text(size=6.5),
        legend.text = element_text(size=7),
        legend.title = element_text(size=8, hjust = 0.5, face = "bold"),
        legend.position = c(0.858, 0.818),
        legend.direction = "vertical",
        legend.key.height = unit(0.01, "cm"),
        legend.key.width = unit(0.01, "cm"),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill="white",color="black")) +
  scale_colour_manual(name='Location', values=c('1 Meter North'='mediumorchid4', 'Underneath Carcass'='royalblue3', 'Soil Control'='mediumseagreen')) + 
  geom_hline(yintercept=mean(genera_vol$Devosia)) +
  geom_hline(yintercept=mean(genera_vol$Devosia)+2*(sd(genera_vol$Devosia)), linetype="dashed") +
  geom_hline(yintercept=mean(genera_vol$Devosia)-2*(sd(genera_vol$Devosia)), linetype="dashed")
p22

#Multiplot
p23 = plot_grid(p20+theme(legend.position="none"),p21+theme(legend.position="none"),p19+theme(legend.position="none"),p22, labels=c("b.", "c.", "d.", "e."), ncol = 2, nrow = 2)
p23
save_plot("genera_long_multiplot.png", p23, base_height = 5.5, base_width = 9,dpi=300)

## Feature stats barplots
genera_stats = read.table("../qiime2_2020_2_analysis/longitudinal/all/genera_lme/feature_stats.tsv", header = TRUE, sep = "\t")

# importance
p24 = ggplot(genera_stats, aes(x=reorder(id,importance), y=importance)) +
  geom_bar(stat="summary",fill="gray80", color="black") +
  coord_flip() +
  labs(x = "Genera", y = "Importance") +
  geom_text(aes(label=id), y= 0.0015, hjust=0, color="black", size=2.5) +
  theme_q2r() + 
  theme(axis.text.y=element_blank(),
        axis.title=element_text(size=9),
        plot.title = element_text(size=11, hjust = 0.5, face = "bold"),
        strip.text = element_text(size=6.5),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill="white",color="black"))
p24

# cumulative avg change
p25 = ggplot(genera_stats, aes(x=reorder(id,importance), y=Cumulative_Avg_Change)) +
  geom_bar(stat="identity", fill="gray75", color="black") +
  geom_hline(yintercept = 0, linetype="solid", color = "black", size=0.5) +
  coord_flip() +
  labs(x = "", y = "Cumulative Average Change") +
  theme_q2r() + 
  theme(axis.text.y=element_blank(),
        axis.title=element_text(size=9),
        plot.title = element_text(size=11, hjust = 0.5, face = "bold"),
        strip.text = element_text(size=6.5),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill="white",color="black"))
p25 

# mean
p26 = ggplot(genera_stats, aes(x=reorder(id,importance), y=Global_Mean)) +
  geom_bar(stat="summary", fill="gray75", color="black") +
  coord_flip() +
  labs(x = "", y = "Global Mean") +
  theme_q2r() + 
  theme(axis.text.y=element_blank(),
        axis.title=element_text(size=9),
        plot.title = element_text(size=11, hjust = 0.5, face = "bold"),
        strip.text = element_text(size=6.5),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill="white",color="black"))
p26

#Multiplot
p27 = plot_grid(p24,p25,p26, labels=c("a.", "", ""), hjust = -0.2, ncol = 3, nrow = 1)
p27
save_plot("genera_stats_long_multiplot.png", p27, base_height = 3, base_width = 9,dpi=300)

# Multiplot combined the long and the barplot stats
p28 = plot_grid(p27,p23, ncol = 1, nrow = 2)
p28
save_plot("genera_feature_multiplot.png", p28, base_height = 8.5, base_width = 9,dpi=300)

#### FINAL multiplot of pcoa and long
# Beta diversity
# Unweighted
p29=unweighted$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  left_join(metadata) %>%
  left_join(shannon_data) %>%
  ggplot(aes(x=-(PC1), y=-(PC2), fill=sample_detail, size=pmi_months)) +
  geom_point(colour="black",pch=21) +
  stat_ellipse(type="t", aes(color = sample_detail)) +
  guides(size = FALSE,
         fill = FALSE,
         color = FALSE) +
  labs(fill = "Location",
       size = "Time (Months)",
       x = paste("PC1: ", round(100*unweighted$data$ProportionExplained[1]), "%"), 
       y = paste("PC2: ", round(100*unweighted$data$ProportionExplained[2]), "%"), 
       title = "Unweighted UniFrac PCoA") +
  theme_bw() +
  theme(axis.text=element_text(size=6.5),
        axis.title=element_text(size=8),
        plot.title = element_text(size=9, hjust = 0.5, face = "bold"),
        legend.position = c(0.855,0.105),
        legend.box = "horizontal",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(color = "black"),
        legend.background = element_rect(fill="white",color="black"),
        legend.key.height = unit(0.01, "cm"),
        legend.key.width = unit(0.01, "cm")) +
  scale_fill_manual(name='Location', values=c('1 Meter North'='mediumorchid4', 'Underneath Carcass'='royalblue3', 'Soil Control'='mediumseagreen')) +
  scale_colour_manual(name='Location', values=c('1 Meter North'='mediumorchid4', 'Underneath Carcass'='royalblue3', 'Soil Control'='mediumseagreen'))
p29

# Beta diversity
# Weighted
p30=weighted$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  left_join(metadata) %>%
  left_join(shannon_data) %>%
  ggplot(aes(x=PC1, y=PC2, fill=sample_detail, size=pmi_months)) +
  geom_point(colour="black",pch=21) +
  stat_ellipse(type="t", aes(color = sample_detail)) +
  guides(size = guide_legend(override.aes = list(linetype = 0)),
         fill = guide_legend(override.aes = list(linetype = 0,size=2)),
         color = FALSE) +
  labs(fill = "Location",
       size = "Time (Months)",
       x = paste("PC1: ", round(100*weighted$data$ProportionExplained[1]), "%"), 
       y = paste("PC2: ", round(100*weighted$data$ProportionExplained[2]), "%"), 
       title = "Weighted UniFrac PCoA") +
  scale_size_continuous(breaks=c(0,0.5,1,6,12,24,72,84,120),
                        limits=c(0,120)) +
  theme_bw() +
  theme(axis.text=element_text(size=6.5),
        axis.title=element_text(size=8),
        plot.title = element_text(size=9, hjust = 0.5, face = "bold"),
        legend.text = element_text(size=7),
        legend.title = element_text(size=8, face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(color = "black"),
        legend.background = element_rect(fill="white",color="white"),
        legend.key.height = unit(0.4, "cm"),
        legend.key.width = unit(0.01, "cm")) +
  scale_fill_manual(name='Location', values=c('1 Meter North'='mediumorchid4', 'Underneath Carcass'='royalblue3', 'Soil Control'='mediumseagreen')) +
  scale_colour_manual(name='Location', values=c('1 Meter North'='mediumorchid4', 'Underneath Carcass'='royalblue3', 'Soil Control'='mediumseagreen'))
p30

p31 = plot_grid(p29,p30+theme(legend.position="none"), labels=c("a.", "b."), ncol = 2, nrow = 1)
p31
legend = get_legend(p30)
p32 = plot_grid(p31, legend, rel_widths = c(4.5, 0.8))
p32
save_plot("beta_multiplot.png", p32, base_height = 5.5, base_width = 7,dpi=300)

## Beta Diversity Longitudinal 
# Unweighted UniFrac
pcoa_vol_unwunifrac = read.table("../qiime2_2020_2_analysis/longitudinal/all/pcoa-vol-unwunifrac.tsv", header = TRUE, sep = "\t", row.names = 1)

p33 = ggplot(pcoa_vol_unwunifrac, aes(x=pmi_months, y=-(Axis1), color=sample_detail)) +
  scale_y_continuous(breaks=pretty_breaks()) +
  stat_summary(geom="errorbar", fun.data=mean_se, width=0, orientation="x") +
  stat_summary(geom="line", fun.data=mean_se, orientation="x", size=1) +
  stat_summary(geom="point", fun.data=mean_se, orientation="x", size=1) +
  labs(x = "Time (Months)", y = "Distance", title = "Unweighted UniFrac PC1") +
  theme_q2r() + 
  theme(axis.text=element_text(size=6.5),
        axis.title.x = element_blank(),
        axis.title=element_text(size=8),
        plot.title = element_text(size=9, hjust = 0.5, face = "bold"),
        strip.text = element_text(size=6.5),
        legend.text = element_text(size=7),
        legend.title = element_text(size=8, hjust = 0.5, face = "bold"),
        legend.position = c(0.874, 0.11),
        legend.direction = "vertical",
        legend.key.height = unit(0.01, "cm"),
        legend.key.width = unit(0.01, "cm"),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill="white",color="black")) +
  scale_colour_manual(name='Location', values=c('1 Meter North'='mediumorchid4', 'Underneath Carcass'='royalblue3', 'Soil Control'='mediumseagreen')) + 
  geom_hline(yintercept=mean(pcoa_vol_unwunifrac$Axis1)) +
  geom_hline(yintercept=mean(pcoa_vol_unwunifrac$Axis1)+2*(sd(pcoa_vol_unwunifrac$Axis1)), linetype="dashed") +
  geom_hline(yintercept=mean(pcoa_vol_unwunifrac$Axis1)-2*(sd(pcoa_vol_unwunifrac$Axis1)), linetype="dashed")
p33

# Weighted UniFrac
pcoa_vol_wunifrac = read.table("../qiime2_2020_2_analysis/longitudinal/all/pcoa-vol-wunifrac.tsv", header = TRUE, sep = "\t", row.names = 1)

p34 = ggplot(pcoa_vol_wunifrac, aes(x=pmi_months, y=Axis1, color=sample_detail)) +
  scale_y_continuous(breaks=pretty_breaks()) +
  stat_summary(geom="errorbar", fun.data=mean_se, width=0, orientation="x") +
  stat_summary(geom="line", fun.data=mean_se, orientation="x", size=1) +
  stat_summary(geom="point", fun.data=mean_se, orientation="x", size=1) +
  labs(x = "Time (Months)", y = "Distance", title = "Weighted UniFrac PC1") +
  theme_q2r() + 
  theme(axis.text=element_text(size=6.5),
        axis.title=element_text(size=8),
        plot.title = element_text(size=9, hjust = 0.5, face = "bold"),
        strip.text = element_text(size=6.5),
        legend.text = element_text(size=7),
        legend.title = element_text(size=8, hjust = 0.5, face = "bold"),
        legend.position = c(0.905, 0.16),
        legend.direction = "vertical",
        legend.key.height = unit(0.01, "cm"),
        legend.key.width = unit(0.01, "cm"),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill="white",color="black")) +
  scale_colour_manual(name='Location', values=c('1 Meter North'='mediumorchid4', 'Underneath Carcass'='royalblue3', 'Soil Control'='mediumseagreen')) + 
  geom_hline(yintercept=mean(pcoa_vol_wunifrac$Axis1)) +
  geom_hline(yintercept=mean(pcoa_vol_wunifrac$Axis1)+2*(sd(pcoa_vol_wunifrac$Axis1)), linetype="dashed") +
  geom_hline(yintercept=mean(pcoa_vol_wunifrac$Axis1)-2*(sd(pcoa_vol_wunifrac$Axis1)), linetype="dashed")
p34

#Multiplot
p35 = plot_grid(p33+theme(legend.position="none"),p34+theme(legend.position="none"), labels=c("c.", "d."), ncol = 1, nrow = 2)
p35

#Multiplot of PCoA and Long
p36 = plot_grid(p32,p35, ncol = 1, nrow = 2)
p36
save_plot("beta_all_multiplot.png", p36, base_height = 5.5, base_width = 7,dpi=300)