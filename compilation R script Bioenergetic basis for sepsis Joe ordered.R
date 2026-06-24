
# Paper "Bioenergetic basis for sepsis""
# April 30, 2026
# By Hessel Peters-Sengers, Joe Butler, Elisa Jentho, and Miguel Soares


# to do:
# 1. not in order of manuscript
# 2. data should be loaded on top, now its at the bottom, with different variable names


library(stringr)
library(ggplot2)
library(eulerr)
library(dplyr)
library(ggrepel)

########## Data for Figure 1F (Venn Euler) ##########


results_df2_conc <- read.csv(
  "https://raw.githubusercontent.com/hesselps/Bioenergetic-basis-of-sepsis/main/Data/data_fig1F.csv"
)

forvenn2 <- data.frame( MOUSE = results_df2_conc$signif_mouse, 
                        HUMAN = results_df2_conc$signif_human )
fitv3 <- euler(forvenn2, shape="ellipse")
col1 <- "orange"
col2 <- "darkolivegreen3"

# plot venn
ve <- plot(fitv3, 
           fills = list(fill=c(col1,col2), alpha=0.5),
           quantities = list(cex=c(1), alpha=0.8, adjust_labels=T), 
           #adjust_labels=T, 
           legend = list(labels = c("Mouse",
                                    "Human"),side="top"))
ve


ggsave("Fig_1F.pdf",ve, width = 5, height = 5)




########## Data for Figure exdat3d  ##########

results_df3 <- read.csv(
  "https://raw.githubusercontent.com/hesselps/Bioenergetic-basis-of-sepsis/main/Data/data_fig1g_exdat3d_suppinfo4.csv"
)

results_df3$class <- sub("[\\(0-9].*", "", results_df3$lipid )
results_df3$class[results_df3$class=="HEX"] <- "HEXCER"

results_df3$class <- as.factor(results_df3$class)

results_df3$lipid_class <- factor(results_df3$class, levels = c("CER",
                                                            "DG",
                                                            "FA",
                                                            "HEXCER",
                                                            "LPC",
                                                            "LPE",
                                                            "PC",
                                                            "PE",
                                                            "PI",
                                                            "SM",
                                                            "ST",
                                                            "TG"))


fit <- lm(results_df3$T_stat_human~results_df3$T_stat_mouse)
summary(fit)


# Extract r and p-value
r <- cor(results_df3$T_stat_mouse, results_df3$T_stat_human, use = "complete.obs")
p <- summary(fit)$coefficients[2,4]

# Format label
p_label <- ifelse(p < 0.0001, "P < 0.0001", paste0("p = ", round(p, 4)))
ann_label <- paste0("R = ", round(r, 2), "\n", p_label)


corrplot2 <- ggplot(data=results_df3, aes(x=T_stat_mouse , y=T_stat_human))+
  geom_point(shape = 21, alpha = 0.7, col="grey20", size = 2, aes( fill=lipid_class))+
  #facet_wrap(.~lipid_class)+
  geom_smooth(method=lm)+
  xlab("Association with sepsis in mouse (t-statistic)")+
  ylab("Association with sepsis in human (t-statistic)")+
  scale_fill_manual(name="Class",
                    values = c("#FF0000", "#FF7600", "#FFEB00", "#9DFF00", "#27FF00", "#00FF4E",
                               "#00FFC4", "#00C4FF", "#004EFF", "#2700FF", "#9D00FF", "#FF00EB"  ), drop=F )+
  annotate("text", 
             x = -Inf, y = Inf,          # bottom-right corner
             hjust = -0.1, vjust = 1.5,  # nudge inward
             label = ann_label,
             size = 4)+
  theme_bw()

corrplot2

ggsave("fig_exdat3d.pdf", corrplot2, width = 10, height = 8)




########## data_fig3d_suppinfo4 ##########


pvalues <- read.csv(
  "https://raw.githubusercontent.com/hesselps/Bioenergetic-basis-of-sepsis/main/Data/data_fig3d_suppinfo4_corrmatrix.csv"
)



corrplot2 <- ggplot(data=results_df3, aes(x=T_stat_mouse, y=T_stat_human))+
  geom_point(shape = 21, alpha = 0.7, col="grey20", size = 2, aes( fill=lipid_class))+
  facet_wrap(.~lipid_class)+
  geom_text(
    data = pvalues,
    aes(x = -6.5, y = 9, label = paste0("P = ", signif(p_value, 2), "\n r = ", signif(cor, 2)) ),
    inherit.aes = FALSE
  ) +
  theme_bw()+
  geom_smooth(method=lm)+
  xlab("Association with sepsis in mouse (t-statistic)")+
  ylab("Association with sepsis in human (t-statistic)")+
  ggtitle("Correlations within Lipid Class ")+
  scale_fill_manual(name="Class",
                    values = c("#FF0000", "#FF7600", "#FFEB00", "#9DFF00", "#27FF00", "#00FF4E",
                               "#00FFC4", "#00C4FF", "#004EFF", "#2700FF", "#9D00FF", "#FF00EB"  ), drop=F )

corrplot2

ggsave("Fig_data_fig3d_suppinfo4.pdf", corrplot2, width = 10, height = 8)


########## data_fig1I, data_ext_fig10_corr, data_ext_fig10_scatter ##########

library(ggrepel)

top_all <- read.csv(
  "https://raw.githubusercontent.com/hesselps/Bioenergetic-basis-of-sepsis/main/Data/data_fig1i.csv"
)


top25msel <- top_all[top_all$label,] 

corvolc1 <- ggplot(data=top_all)+
  geom_jitter(aes(x=pearson_r, y=-log10(p_value), fill = lipid_class), shape = 21, alpha = 0.7, col="grey20", size = 2,width = 0.02)+
  theme_bw()+
  xlab("Correlation with SOFA (Pearson's r)")+
  ylab("-log10(p)")+
  geom_text_repel(
    data = top25msel,
    aes(label = lipid, x=pearson_r, y=-log10(p_value)),                  
    size = 3,
    max.overlaps = Inf,
    min.segment.length = 0,
    box.padding = 1.1,
    segment.size  = 0.3,
    segment.color = "grey50",
    #box.padding = 0.4
  )+
  scale_y_continuous(limits = c(-2, 12))+
  scale_x_continuous(limits = c(-0.6, 0.6))+
  scale_fill_manual(values = c("#FF0000", "#FF7600", "#FFEB00", "#9DFF00", "#27FF00", "#00FF4E",
                               "#00FFC4", "#00C4FF", "#004EFF", "#2700FF", "#9D00FF", "#FF00EB" ,"#FF0076"), drop=F )+
  labs(fill = "Class")


print(corvolc1)


ggsave("Fig_1i.pdf", corvolc1, width = 7, height = 5)

# #save data for Nature
# head(top_all)
# top_all$pearson_r <- top_all$r
# top_all$p_value <- top_all$p
# top_all$lipid_class <- top_all$class2
# setwd("D:/.shortcut-targets-by-id/10J9GdrycmUMcd6u97iLZyH8Psxxr964x/Lipidomics Mouse/Datasets Figures")
# write.csv(top_all[,c("lipid","pearson_r", "p_value", "lipid_class"  )], "data_fig1I.csv", row.names = F)
# 
# write.csv(top_all[top_all$lipid %in% lipid_species_vector,c("lipid","pearson_r", "p_value", "lipid_class"  )], "data_ext_fig10_corr.csv", row.names = F)
# 
# df_aggregated2_10 <- df_aggregated2[, c("SOFAtot", lipid_species_vector), drop = FALSE]
# write.csv(df_aggregated2_10, "data_ext_fig10_scatter.csv", row.names = F)


########## ########## data_extended_fig5A ##########

t_test_pvalue_log2FC <- function(row) {
  # Extract "Gr1" and "Gr2" columns as numeric vectors
  group1 <- as.numeric(row[grepl("Gr1", colnames(df_aggregated))])
  group2 <- as.numeric(row[grepl("Gr2", colnames(df_aggregated))])
  
  # Perform t-test
  test_result <- t.test(group2, group1 )
  
  # Calculate log2 fold change (mean of Gr2 - mean of Gr1)
  mean_gr1 <- mean(group1, na.rm = TRUE)
  mean_gr2 <- mean(group2, na.rm = TRUE)
  log2FC <- log2(mean_gr2 / mean_gr1)
  
  # Extract p-value and return it along with log2FC
  p_value <- test_result$p.value
  
  # Return as a named vector
  return(c(p_value = p_value, log2FC = log2FC))
}

# Apply the function row-wise to get p-values and log2FCs
results <- t(apply(df_aggregated, 1, t_test_pvalue_log2FC))

# Convert to a dataframe for easy viewing
results_df <- as.data.frame(results)

# Display the results
head(results_df)

results_df$lipid <- df_aggregated$MARS_lipids_converted
results_df$t_stat_mars <- df_aggregated$t_statistics
results_df$class <- df_aggregated$class

library(ggrepel)

results_df$padj <- p.adjust(results_df$p_value,  method="BH")
results_df$neg_log10_pvalue <- -log10(results_df$padj)
results_df$color <- ifelse(results_df$padj < 0.05 & results_df$log2FC > 0, "upregulated",
                           ifelse(results_df$padj< 0.05 & results_df$log2FC < 0, "downregulated", "non-significant"))


top_biomarkers <- results_df[results_df$lipid %in% c("PS(36:1)", "PS(38:6)", "FA19:1","FA22:5") |
                               results_df$lipid %in% c( "TG(61:2)", "LPC(19:0)","PE(O-38:6)",  "CER(42:0;O2)",  "PC(36:0)"),]

top_biomarkers

# non-sign classes
results_df$class2 <- ifelse(results_df$padj>0.05, "non-sig", results_df$class)
results_df$class2 <- ifelse(results_df$class=="LPG"|
                              results_df$class=="LPI"|
                              results_df$class=="LPS"|
                              results_df$class=="PA"|
                              results_df$class=="PG"|
                              results_df$class=="PS", "Other", results_df$class2)

results_df$class2 <- factor(results_df$class2, levels = c("CER",
                                                          "DG",
                                                          "FA",
                                                          "HEXCER",
                                                          "LPC",
                                                          "LPE",
                                                          "PC",
                                                          "PE",
                                                          "PI",
                                                          "SM",
                                                          "ST",
                                                          "TG",
                                                          "Other",
                                                          "non-sig"))

# Create the volcano plot
volcano_plot <- ggplot(results_df, aes(x = log2FC, y = neg_log10_pvalue)) +
  geom_point(aes(fill = class2), shape = 21, alpha = 0.7, col="grey20", size = 2) + # Color points based on significance
  scale_fill_manual(values = c("#FF0000", "#FF7600", "#FFEB00", "#9DFF00", "#27FF00", "#00FF4E",
                               "#00FFC4", "#00C4FF", "#004EFF", "#2700FF", "#9D00FF", "#FF00EB" ,"#FF0076" ,"grey60"))+
  #scale_color_manual(values = c("upregulated" = "red", "downregulated" = "blue", "non-significant" = "black")) + # Assign colors
  labs(x = "Log2 Fold Change", y = "-log10(BH adj p-value)", title = "CAP-sepsis vs Controls") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank()) +
  
  # Add labels to the top 20 upregulated and downregulated biomarkers
  #geom_text_repel(data = top_biomarkers, aes(label = lipid), size = 3, box.padding = 0.5, max.overlaps = Inf)
  geom_text_repel(aes(label = lipid),
                  data          = subset(top_biomarkers, log2FC > 0),
                  nudge_x       = 3 - subset(top_biomarkers, log2FC> 0)$log2FC,
                  #segment.size  = 0.2,
                  segment.color = "grey50",
                  direction     = "y",
                  hjust         = 0, size=3.5
  ) +
  geom_text_repel(aes(label = lipid),
                  data          =  subset(top_biomarkers, log2FC <0),
                  nudge_x       = -3.2 - subset(top_biomarkers, log2FC< 0)$log2FC,
                  #segment.size  = 0.2,
                  segment.color = "grey50",
                  direction     = "y",
                  hjust         = 1, size=3.5
  ) 

# Display the plot
print(volcano_plot)

head(results_df)
results_df$lipid_class <- results_df$class2
setwd("D:/.shortcut-targets-by-id/10J9GdrycmUMcd6u97iLZyH8Psxxr964x/Lipidomics Mouse/Datasets Figures")
write.csv(results_df[,c("lipid","log2FC", "neg_log10_pvalue","lipid_class" )], "data_extended_fig5A.csv", row.names = F)

###############################################
 ########## exdat3b ##########################
###############################################



volc1dat <- read.csv(
  "https://raw.githubusercontent.com/hesselps/Bioenergetic-basis-of-sepsis/main/Data/data_exdat3b.csv"
)


volc1dat$class2 <- volc1dat$lipid_class
volc1dat$class2 <- factor(volc1dat$class2, levels = c("CER",
                                                      "DG",
                                                      "FA",
                                                      "HEXCER",
                                                      "LPC",
                                                      "LPE",
                                                      "PC",
                                                      "PE",
                                                      "PI",
                                                      "SM",
                                                      "ST",
                                                      "TG",
                                                      "Other",
                                                      "non-sig"))



volc1dat$p <- 10^(-volc1dat$minus_log10_pvalue)
volc1dat$FC <- volc1dat$log2FC

top_upregulated <- volc1dat %>%
  filter(FC < 0 & p < 0.006 & class2=="LPC") %>%
  arrange(p) %>%
  head(8)

top_downregulated <- volc1dat %>%
  filter(FC < 0 & p < 0.006 & class2=="TG") %>%
  arrange(p) %>%
  head(11)

top_other_pos <- volc1dat %>%
  filter(FC > 0 & p < 0.01) %>%
  arrange(p) %>%
  head(10)

top_other_neg <- volc1dat %>%
  filter(FC < 0 & p < 0.01 & class2!="LPC" & class2!="TG") %>%
  arrange(p) %>%
  head(5)


top_biomarkers <- rbind(top_upregulated, top_downregulated, top_other_pos, top_other_neg)

top_biomarkers <- volc1dat[volc1dat$lipid %in% c("PS(36:1)", "PS(38:6)", "FA(19:1)","FA(22:5)") |
                             volc1dat$lipid %in% c( "TG(61:2)", "LPC(19:0)","PE(O-38:6)",  "Cer(42:0)",  "PC(36:0)"), ]
top_biomarkers <- top_biomarkers[!(top_biomarkers$class2=="non-sig"),]
top_biomarkers
table(volc1dat$class2)


volcano_plot_ED3_b <- ggplot(volc1dat, aes(x = FC, y = -log10(p)  )) +
  geom_point(aes(fill = class2), shape = 21, alpha = 0.7, col="grey20", size = 2) + # Color points based on significance
  
  
  scale_fill_manual(values = c("#FF0000", "#FF7600", "#FFEB00", "#9DFF00", "#27FF00", "#00FF4E",
                               "#00FFC4", "#00C4FF", "#004EFF", "#2700FF", "#9D00FF", "#FF00EB" ,"#FF0076" ,"grey60"), drop=F )+
  
  
  #scale_color_manual(values = c("upregulated" = "red", "downregulated" = "blue", "non-significant" = "black")) + # Assign colors
  labs(x = "Log2 Fold Change", y = "-log10(p-value)", title = "CLP v Ctrl") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank()) +
  
  # Add labels to the top 20 upregulated and downregulated biomarkers
  #geom_text_repel(data = top_biomarkers, aes(label = lipid), size = 3, box.padding = 0.5, max.overlaps = Inf)
  geom_text_repel(aes(label = lipid),
                  data          = subset(top_biomarkers, FC > 0),
                  nudge_x       = 5 - subset(top_biomarkers, FC> 0)$FC,
                  nudge_y       = -1.5 -log10(subset(top_biomarkers, FC> 0)$p),
                  #segment.size  = 0.2,
                  segment.color = "grey50",
                  direction     = "y",
                  hjust         = 0, size=3.5
  ) +
  geom_text_repel(aes(label = lipid),
                  data          =  subset(top_biomarkers, FC <0),
                  nudge_x       = -10 - subset(top_biomarkers, FC< 0)$FC,
                  
                  #segment.size  = 0.2,
                  segment.color = "grey50",
                  direction     = "y",
                  hjust         = 0, size=3.5
  ) 

volcano_plot_ED3_b 


ggsave("exdat3b.pdf", volcano_plot_ED3_b, width = 5, height = 6)


##### Supplementary_Data_3a (faceted)


volc1dat$class <- sub("[\\(0-9].*", "", volc1dat$lipid )
volc1dat$class[volc1dat$class=="Hex"] <- "HexCer"
volc1dat$class[volc1dat$class=="NATau"] <- "FA"

### remove Other

volc1dat_ss <- volc1dat[volc1dat$class %in%  c("Cer","DG",
                                             "FA",
                                             "HexCer",
                                             "LPC",
                                             "LPE",
                                             "PC",
                                             "PE",
                                             "PI",
                                             "SM",
                                             "ST",
                                             "TG"),]

volc1dat_ss$class3 <- volc1dat_ss$class
volc1dat_ss$class3[volc1dat_ss$lipid_class=="non-sig"] <- "non-sig"


volc1dat_ss$class3 <- factor(volc1dat_ss$class3, levels = c("Cer",
                                                      "DG",
                                                      "FA",
                                                      "HexCer",
                                                      "LPC",
                                                      "LPE",
                                                      "PC",
                                                      "PE",
                                                      "PI",
                                                      "SM",
                                                      "TG",
                                                      "non-sig"))


labels <- c("Cer(41:0)", "Cer(43:0)", "FA(16:1)", "NATau(18:1)",  "NATau(22:6)",
            "HexCer(34:2)",  "HexCer(42:1)", "LPC(20:2)",  "LPC(20:1)",  "LPC(22:3/0:0)",
            "PC(32:1)", "PE(40:8)", "PE(34:3)", "PE(P-34:2)", "PI(32:0)", "PI(36:1)", "TG(58:0)",
            "TG(O-54:4)","TG(O-52:3)" )

volc1dat_ss$tolab <- volc1dat_ss$lipid %in% labels & volc1dat_ss$minus_log10_pvalue > 2


volcano_plot_Supp3_a <- ggplot(volc1dat_ss, aes(x = FC, y = -log10(p)  )) +
  geom_point(aes(fill = class3), shape = 21, alpha = 0.7, col="grey20", size = 2) + # Color points based on significance
  
  
  scale_fill_manual(values = c("#FF0000", "#FF7600", "#FFEB00", "#9DFF00", "#27FF00", "#00FF4E",
                               "#00FFC4", "#00C4FF", "#004EFF", "#2700FF", "#FF0076" ,"grey60"), drop=F )+
  
  
  #scale_color_manual(values = c("upregulated" = "red", "downregulated" = "blue", "non-significant" = "black")) + # Assign colors
  labs(x = "Log2 Fold Change", y = "-log10(p-value)", title = "CLP v Ctrl") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank()) +
  
  # Add labels to the top 20 upregulated and downregulated biomarkers
  #geom_text_repel(data = top_biomarkers, aes(label = lipid), size = 3, box.padding = 0.5, max.overlaps = Inf)
  geom_text_repel(aes(label = lipid),
                  data          = subset(volc1dat_ss, tolab),
                  #segment.size  = 0.2,
                  min.segment.length = 0,
                  max.overlaps = Inf,
                  segment.color = "grey50",
                  direction     = "y",
                  hjust         = 0, size=3.5
  ) +
  
  facet_wrap(.~class)+
  scale_x_continuous(
    breaks = c(-5, 0, 5)
  ) 

volcano_plot_Supp3_a

ggsave("Supplementary_Data_3a.pdf", volcano_plot_Supp3_a, width = 10, height = 8)

########## data_exdat3c ##########



volc3cdat <- read.csv(
  "https://raw.githubusercontent.com/hesselps/Bioenergetic-basis-of-sepsis/main/Data/data_exdat3c.csv"
)

# t_test_pvalue_log2FC <- function(row) {
#   # Extract "Gr1" and "Gr2" columns as numeric vectors
#   group1 <- as.numeric(row[grepl("Gr1", colnames(df_aggregated))])
#   group2 <- as.numeric(row[grepl("Gr2", colnames(df_aggregated))])
#   
#   # Perform t-test
#   test_result <- t.test(group2, group1 )
#   
#   # Calculate log2 fold change (mean of Gr2 - mean of Gr1)
#   mean_gr1 <- mean(group1, na.rm = TRUE)
#   mean_gr2 <- mean(group2, na.rm = TRUE)
#   log2FC <- log2(mean_gr2 / mean_gr1)
#   
#   # Extract p-value and return it along with log2FC
#   p_value <- test_result$p.value
#   
#   # Return as a named vector
#   return(c(p_value = p_value, log2FC = log2FC))
# }
# 
# # Apply the function row-wise to get p-values and log2FCs
# results <- t(apply(df_aggregated, 1, t_test_pvalue_log2FC))

# Convert to a dataframe for easy viewing
results_df <- as.data.frame(volc3cdat)

# Display the results
head(results_df)

# results_df$lipid <- df_aggregated$MARS_lipids_converted
# results_df$t_stat_mars <- df_aggregated$t_statistics
# results_df$class <- df_aggregated$class

library(ggrepel)


top_biomarkers <- results_df[results_df$lipid %in% c("PS(36:1)", "PS(38:6)", "FA19:1","FA22:5") |
                               results_df$lipid %in% c( "TG(61:2)", "LPC(19:0)","PE(O-38:6)",  "CER(42:0;O2)",  "PC(36:0)"),]


results_df$class2 <- as.factor(results_df$lipid_class)
results_df$class2 <- factor(results_df$class2, levels = c("CER",
                                                          "DG",
                                                          "FA",
                                                          "HEXCER",
                                                          "LPC",
                                                          "LPE",
                                                          "PC",
                                                          "PE",
                                                          "PI",
                                                          "SM",
                                                          "ST",
                                                          "TG",
                                                          "Other",
                                                          "non-sig"))

# Create the volcano plot
volcano_plot <- ggplot(results_df, aes(x = log2FC, y = neg_log10_pvalue)) +
  geom_point(aes(fill = class2), shape = 21, alpha = 0.7, col="grey20", size = 2) + # Color points based on significance
  scale_fill_manual(values = c("#FF0000", "#FF7600", "#FFEB00", "#9DFF00", "#27FF00", "#00FF4E",
                               "#00FFC4", "#00C4FF", "#004EFF", "#2700FF", "#9D00FF", "#FF00EB" ,"#FF0076" ,"grey60"))+
  #scale_color_manual(values = c("upregulated" = "red", "downregulated" = "blue", "non-significant" = "black")) + # Assign colors
  labs(x = "Log2 Fold Change", y = "-log10(BH adj p-value)", title = "CAP-sepsis vs Controls") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank()) +
  
  # Add labels to the top 20 upregulated and downregulated biomarkers
  #geom_text_repel(data = top_biomarkers, aes(label = lipid), size = 3, box.padding = 0.5, max.overlaps = Inf)
  geom_text_repel(aes(label = lipid),
                  data          = subset(top_biomarkers, log2FC > 0),
                  nudge_x       = 3 - subset(top_biomarkers, log2FC> 0)$log2FC,
                  #segment.size  = 0.2,
                  segment.color = "grey50",
                  direction     = "y",
                  hjust         = 0, size=3.5
  ) +
  geom_text_repel(aes(label = lipid),
                  data          =  subset(top_biomarkers, log2FC <0),
                  nudge_x       = -3.2 - subset(top_biomarkers, log2FC< 0)$log2FC,
                  #segment.size  = 0.2,
                  segment.color = "grey50",
                  direction     = "y",
                  hjust         = 1, size=3.5
  ) 

# Display the plot
print(volcano_plot)



ggsave("exdat3c.pdf", volcano_plot, width = 5, height = 6)


#### faceted Supplementary_Data_3_b

results_df$class <- sub("[\\(0-9].*", "", results_df$lipid )
results_df$class[results_df$class=="HEX"] <- "HEXCER"


### remove Other

results_df_ss <- results_df[results_df$class %in%  c("CER","DG",
                                               "FA",
                                               "HEXCER",
                                               "LPC",
                                               "LPE",
                                               "PC",
                                               "PE",
                                               "PI",
                                               "SM",
                                               "ST",
                                               "TG"),]

results_df_ss$class3 <- results_df_ss$class
results_df_ss$class3[results_df_ss$lipid_class=="non-sig"] <- "non-sig"


results_df_ss$class3 <- factor(results_df_ss$class3, levels = c("CER",
                                                            "DG",
                                                            "FA",
                                                            "HEXCER",
                                                            "LPC",
                                                            "LPE",
                                                            "PC",
                                                            "PE",
                                                            "PI",
                                                            "SM",
                                                            "ST",
                                                            "TG",
                                                            "non-sig"))

labels <-  c("PS(36:1)", "PS(38:6)", "FA19:1","FA22:5", "TG(61:2)", "LPC(19:0)","PE(O-38:6)",  "CER(42:0;O2)",  "PC(36:0)")

results_df_ss$tolab <- results_df_ss$lipid %in% labels & results_df_ss$neg_log10_pvalue > 2


volcano_plot_Supp3_b <- ggplot(results_df_ss, aes(x = log2FC, y = neg_log10_pvalue  )) +
  geom_point(aes(fill = class3), shape = 21, alpha = 0.7, col="grey20", size = 2) + # Color points based on significance
  
  
  scale_fill_manual(values = c("#FF0000", "#FF7600", "#FFEB00", "#9DFF00", "#27FF00", "#00FF4E",
                               "#00FFC4", "#00C4FF", "#004EFF", "#2700FF","#FF00EB" , "#FF0076" ,"grey60"), drop=F )+
  
  
  #scale_color_manual(values = c("upregulated" = "red", "downregulated" = "blue", "non-significant" = "black")) + # Assign colors
  labs(x = "Log2 Fold Change", y = "-log10(p-value)", title = "CAP v Ctrl") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank()) +
  
  # Add labels to the top 20 upregulated and downregulated biomarkers
  #geom_text_repel(data = top_biomarkers, aes(label = lipid), size = 3, box.padding = 0.5, max.overlaps = Inf)
  geom_text_repel(aes(label = lipid),
                  data          = subset(results_df_ss, tolab),
                  #segment.size  = 0.2,
                  min.segment.length = 0,
                  max.overlaps = Inf,
                  segment.color = "grey50",
                  direction     = "y",
                  hjust         = 0, size=3.5
  ) +
  
  facet_wrap(.~class)+
  scale_x_continuous(
    breaks = c(-5, 0, 5)
  ) 

volcano_plot_Supp3_b

ggsave("Supplementary_Data_3b.pdf", volcano_plot_Supp3_b, width = 10, height = 8)

##################  data_exdat3f


volc2dat <- read.csv(
  "https://raw.githubusercontent.com/hesselps/Bioenergetic-basis-of-sepsis/main/Data/data_exdat3f.csv"
)

volc2dat$class_col <- volc2dat$lipid_class

volc2dat$class_col <- factor(volc2dat$class_col, levels = c("CER",
                                                            "DG",
                                                            "FA",
                                                            "HEXCER",
                                                            "LPC",
                                                            "LPE",
                                                            "PC",
                                                            "PE",
                                                            "PI",
                                                            "SM",
                                                            "ST",
                                                            "TG",
                                                            "Other",
                                                            "non-sig"))


volc2dat$p <- 10^(-volc2dat$neg_log10_p_value)


#### MAIN (non facetted)
labels <- c("ST(18:2)","PC(40:5)", "SM(40:2;O2)", "ST(20:4)","SM(42:1;O2)")

top_biomarkers <- volc2dat[volc2dat$lipid_name %in% labels,]

volcano_plot_1F <- ggplot(volc2dat, aes(x = log_2_foldchange, y = neg_log10_p_value  )) +
  geom_vline(xintercept = 0, col="gray70")+
  geom_point(aes(fill = class_col), shape = 21, alpha = 0.7, col="grey20", size = 2) + # Color points based on significance

  #scale_fill_manual(values = c("#FF0000", "#FF7600", "#FFEB00", "#27FF00", "#00FF4E",
  #                             "#00FFC4", "#00C4FF", "#004EFF", "#2700FF", "#FF00EB" ,"#FF0076" ,"grey60"), drop=T )+
  
  scale_fill_manual(values = c("#FF0000", "#FF7600", "#FFEB00", "#9DFF00", "#27FF00", "#00FF4E",
                               "#00FFC4", "#00C4FF", "#004EFF", "#2700FF", "#9D00FF", "#FF00EB", "#FF0076", "grey50"), drop=F )+
  
  
  #scale_color_manual(values = c("upregulated" = "red", "downregulated" = "blue", "non-significant" = "black")) + # Assign colors
  labs(x = "Log2 Fold Change", y = "-log10(p-value)", title = "AdipoQD/D v fl/fl") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank())+
 
  geom_text_repel(aes(label = lipid_name),
                  data          =  subset(top_biomarkers, log_2_foldchange <0),
                  nudge_x       = -4 - subset(top_biomarkers, log_2_foldchange< 0)$log_2_foldchange,
                  nudge_y = 1,
                  #segment.size  = 0.2,
                  segment.color = "grey50",
                  direction     = "y",
                  hjust         = 0, size=2.5
  ) +
  scale_y_continuous(limits = c(0,4))

volcano_plot_1F

ggsave("exdat3f.pdf", volcano_plot_1F, width = 5, height = 6)





#### suppinfo 4b and 4c #######


results_df2_2 <- read.csv(
  "https://raw.githubusercontent.com/hesselps/Bioenergetic-basis-of-sepsis/main/Data/data_suppinfo4b_4c.csv"
)


### remove disconcordant overlap lipids

results_df2_3 <- results_df2_2[ results_df2_2$overlap | 
                                  (results_df2_2$signif_mouse & !results_df2_2$signif_human) |
                                  (!results_df2_2$signif_mouse & results_df2_2$signif_human), ]

forvenn3 <- data.frame( MOUSE = results_df2_3$signif_mouse , 
                        HUMAN = results_df2_3$signif_human )
fitv3 <- euler(forvenn3, shape="ellipse")
col1 <- "darkorange"
col2 <- "darkolivegreen3"

# plot venn
ve <- plot(fitv3, 
           fills = list(fill=c(col1,col2), alpha=0.5),
           quantities = list(cex=c(1), alpha=0.8, adjust_labels=T), 
           #adjust_labels=T, 
           legend = list(labels = c("Mouse",
                                    "Human"),side="top"))
ve


###### 4c #######

labels <- c(
  "ST(20:4)",
  "FA25:1",
  "FA20:0",
  "PC(34:1)",
  "CER(41:2;O2)",
  "PC(40:5)",
  "PC(36:5)",
  "PC(40:1)",
  "PC(36:2)",
  "LPE(20:3)",
  "FA26:1",
  "PC(42:5)",
  "PI(38:2)",
  "LPE(18:1)",
  "PC(34:0)",
  "SM(33:1;O2)",
  "LPE(18:3)",
  "PI(36:1)",
  "ST(18:2)",
  "SM(40:2;O2)",
  "SM(32:2;O2)",
  "PC(38:2)",
  "SM(42:1;O2)",
  "SM(41:1;O2)",
  "PC(38:3)",
  "SM(40:1;O2)"
)

lipid_labels <- results_df2_2[results_df2_2$lipid %in% labels,]

corrplot <- ggplot(data=results_df2_2, aes(x=T_stat_mouse_sepsis_AdipoQD_D_v_fl_fl, y=T_stat_human_sepsis))+
  geom_hline(yintercept = 0, colour = "grey50")+
  geom_vline(xintercept = 0, colour = "grey50")+
  geom_point( aes( fill = overlap), shape = 21, alpha = 0.7, col="grey20", size = 2)+
  # geom_text_repel(aes(label=label2), size=3, box.padding = 2, max.overlaps = Inf,
  #                 segment.color = "lightgrey",  # Light grey lines
  #                 segment.size = 0.3)+
  scale_y_continuous(limits = c(-15,10))+
  scale_x_continuous(limits = c(-15,7))+
  geom_text_repel(aes(label = lipid),
                  data          =  lipid_labels,
                  nudge_x       = -15 - lipid_labels$T_stat_mouse_sepsis_AdipoQD_D_v_fl_fl,
                  segment.size  = 0.3,
                  segment.color = "grey80",
                  direction     = "y",
                  hjust         = 0, size=3
  )+
  theme_bw()+
  #geom_smooth(method=lm)+
  xlab("Association with ATLG-ko in mouse (t-statistic)")+
  ylab("Association with sepsis in human (t-statistic)")+
  #ggtitle("Correlation, P < 0.0001 ")+
  scale_fill_manual(name="Class",
                    values = c("grey50", "#EBDD9C" ), drop=F )


corrplot

ggsave("suppinfo4c.pdf", corrplot, width = 5, height = 6)


#############################################################################
### suppinfo5 ( facetted )
#############################################################################

volc2dat <- read.csv(
  "https://raw.githubusercontent.com/hesselps/Bioenergetic-basis-of-sepsis/main/Data/data_suppinfo5.csv"
)

volc2dat$lipid_class_fill <- factor(volc2dat$lipid_class_fill, levels = c("CER",
                                                                          "DG",
                                                                          "FA",
                                                                          "HEXCER",
                                                                          "LPC",
                                                                          "LPE",
                                                                          "PC",
                                                                          "PE",
                                                                          "PI",
                                                                          "SM",
                                                                          "ST",
                                                                          "TG",
                                                                          "non-sig"))


volc2dat_labeled <- volc2dat[volc2dat$label, ]

volcano_plot_1F_FACET <- ggplot(volc2dat, aes(x = log_2_foldchange, y = neg_log10_p_value  )) +
  geom_vline(xintercept = 0, col="gray70")+
  geom_point(aes(fill = lipid_class_fill ), shape = 21, alpha = 0.7, col="grey20", size = 2) + # Color points based on significance
  facet_wrap(.~lipid_class)+
  #scale_fill_manual(values = c("#FF0000", "#FF7600", "#FFEB00", "#27FF00", "#00FF4E",
  #                             "#00FFC4", "#00C4FF", "#004EFF", "#2700FF", "#FF00EB" ,"#FF0076" ,"grey60"), drop=T )+
  
  scale_fill_manual(values = c("#FF0000", "#FF7600", "#FFEB00", "#9DFF00", "#27FF00", "#00FF4E",
                               "#00FFC4", "#00C4FF", "#004EFF", "#2700FF", "#9D00FF", "#FF00EB" , "grey50"), drop=F )+
  
  
  #scale_color_manual(values = c("upregulated" = "red", "downregulated" = "blue", "non-significant" = "black")) + # Assign colors
  labs(x = "Log2 Fold Change", y = "-log10(p-value)", title = "AdipoQD/D v fl/fl") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank())+
  geom_text_repel(aes(label = lipid_name),
                  min.segment.length   = 0,
                  box.padding = 0.4,
                  data          = volc2dat_labeled,
                  #nudge_y       = 4.5 - 2*volc2dat_labeled$neg_log10_p_value,
                  #nudge_y       = 0.2,
                  #segment.size  = 0.2,
                  segment.color = "grey50",
                  size=2.5
  ) +
  scale_y_continuous(limits = c(0,5))

volcano_plot_1F_FACET

ggsave("suppinfo5.pdf", volcano_plot_1F_FACET, width = 5, height = 6)
