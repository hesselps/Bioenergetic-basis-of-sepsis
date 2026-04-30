# Paper "Bioenergetic basis for sepsis""
# April 30, 2026
# By Hessel Peters-Sengers, Joe Butler, and Elisa Jentho


# to do:
# 1. not in order of manuscript
# 2. data should be loaded on top, now its at the bottom, with different variable names


########## data_extended_fig5A ##########

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



########## data_extended_fig5B ##########



top_others <- results_df %>%
  filter(log2FC > 0 &neg_log10_pvalue>10) %>%
  arrange(p_value)

top_others2 <- results_df %>%
  filter(log2FC < 0 &neg_log10_pvalue>10) %>%
  arrange(p_value)

# Combine the top 20 upregulated and downregulated into one dataframe
top_biomarkers <- rbind(top_others, top_others2)

volcano_plot_facet <- ggplot(results_df[!(results_df$class2=="Other"),], aes(x = log2FC, y = neg_log10_pvalue)) +
  geom_point(aes(fill = class), shape = 21, alpha = 0.7, col="grey20", size = 2) + # Color points based on significance
  scale_fill_manual(values = c("#FF0000", "#FF7600", "#FFEB00", "#9DFF00", "#27FF00", "#00FF4E",
                               "#00FFC4", "#00C4FF", "#004EFF", "#2700FF", "#9D00FF", "#FF00EB" ))+
  #scale_color_manual(values = c("upregulated" = "red", "downregulated" = "blue", "non-significant" = "black")) + # Assign colors
  # Additional layer to overwrite points with -log10 adjusted p-value below 0.05 as grey
  geom_point(data = subset(results_df[!(results_df$class2=="Other"),], neg_log10_pvalue < -log10(0.05)),
             shape = 21, fill = "grey", color = "grey20", size = 2, alpha = 0.7) +
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
                  hjust         = 0, size=2.5
  ) +
  geom_text_repel(aes(label = lipid),
                  data          =  subset(top_biomarkers, log2FC <0),
                  nudge_x       = -3.2 - subset(top_biomarkers, log2FC< 0)$log2FC,
                  #segment.size  = 0.2,
                  segment.color = "grey50",
                  direction     = "y",
                  hjust         = 1, size=2.5
  ) +
  facet_wrap(.~class)

print(volcano_plot_facet)

setwd("D:/Elisa Jena/")
ggsave("volcano_plot_facetted.pdf",volcano_plot_facet, width = 10, height = 10)

head(results_df)
results_df$lipid_class <- results_df$class2
setwd("D:/.shortcut-targets-by-id/10J9GdrycmUMcd6u97iLZyH8Psxxr964x/Lipidomics Mouse/Datasets Figures")
write.csv(results_df[!(results_df$class2=="Other"),c("lipid","log2FC", "neg_log10_pvalue","lipid_class")], "data_extended_fig5B.csv", row.names = F)



########## data_extended_fig4D ##########



volc1dat <- read.csv("J:/.shortcut-targets-by-id/10J9GdrycmUMcd6u97iLZyH8Psxxr964x/Lipidomics Mouse/Analyses in progress by H and J/Batch1_lipid_FC.tsv", sep="\t")
plot(volc1dat$FC, -log10(volc1dat$p))

volc1dat$metabolites <- sub("\\|.*", "", volc1dat$`Metabolite.name`)
volc1dat$metabolites2 <- sub("\\;.*", "", volc1dat$metabolites)
volc1dat$metabolites3 <- gsub(" ", "", volc1dat$metabolites2)
#volc1dat$class <- gsub("[0-9:]", "", volc1dat$metabolites3)

volc1dat$metabolites4 <- gsub(" ", "(", volc1dat$metabolites2)
volc1dat$metabolites4 <- paste0( volc1dat$metabolites4, ")")

table(volc1dat$group_manual)

volc1dat$class2 <- "Other"

volc1dat$class2[volc1dat$group_manual=="Ceramides"] <- "CER" 
volc1dat$class2[grepl("Fatty", volc1dat$group_manual)] <- "FA" 
volc1dat$class2[grepl("Diradylglycerols", volc1dat$group_manual)] <- "DG" 
volc1dat$class2[grepl("Triradylcglycerols", volc1dat$group_manual)] <- "TG" 
volc1dat$class2[grepl("Glycerophosphocholines", volc1dat$group_manual)] <- "PC" 
volc1dat$class2[grepl("Glycerophosphoethanolamines", volc1dat$group_manual)] <- "PE" 
volc1dat$class2[grepl("Glycerophosphoglycerols", volc1dat$group_manual)] <- "PG" 
volc1dat$class2[grepl("Glycerophosphoinositols", volc1dat$group_manual)] <- "PI" 
volc1dat$class2[grepl("Sphingomyelins", volc1dat$group_manual)] <- "SM" 
volc1dat$class2[grepl("St", volc1dat$group_manual)] <- "ST" 
volc1dat$class2[grepl("Neutral glycosphingolipids|Phosphosphingolipids", volc1dat$group_manual)] <- "SM" 
volc1dat$class2[grepl("LPC", volc1dat$Metabolite.name)] <- "LPC" 
volc1dat$class2[grepl("LPE", volc1dat$Metabolite.name)] <- "LPE" 
volc1dat$class2[grepl("HexCer", volc1dat$Metabolite.name)] <- "HEXCER" 

table(volc1dat$class2)

volc1dat$class2[volc1dat$p>0.05] <- "non-sig"

table(volc1dat$class2)

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


volc1dat$lipid <- volc1dat$metabolites4


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
# top_others <- results_df %>%
#   filter(log2FC > 0 &neg_log10_pvalue>7) %>%
#   arrange(p_value) %>%
#   head(10)
# 
# top_others2 <- results_df %>%
#   filter(log2FC < 0 &neg_log10_pvalue>7) %>%
#   arrange(p_value) %>%
#   head(5)

# Combine the top 20 upregulated and downregulated into one dataframe
#top_biomarkers <- rbind(top_upregulated, top_downregulated, top_others, top_others2)

top_biomarkers <- rbind(top_upregulated, top_downregulated, top_other_pos, top_other_neg)

top_biomarkers <- volc1dat[volc1dat$lipid %in% c("PS(36:1)", "PS(38:6)", "FA(19:1)","FA(22:5)") |
                             volc1dat$lipid %in% c( "TG(61:2)", "LPC(19:0)","PE(O-38:6)",  "Cer(42:0)",  "PC(36:0)"), ]
top_biomarkers <- top_biomarkers[!(top_biomarkers$class2=="non-sig"),]
top_biomarkers
table(volc1dat$class2)


volcano_plot_1B <- ggplot(volc1dat, aes(x = FC, y = -log10(p)  )) +
  geom_point(aes(fill = class2), shape = 21, alpha = 0.7, col="grey20", size = 2) + # Color points based on significance
  
  #scale_fill_manual(values = c("#FF0000", "#FF7600", "#FFEB00", "#27FF00", "#00FF4E",
  #                             "#00FFC4", "#00C4FF", "#004EFF", "#2700FF", "#FF00EB" ,"#FF0076" ,"grey60"), drop=T )+
  
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

volcano_plot_1B 
setwd("D:/Elisa Jena")
ggsave("volcano_plot_1Bnew.pdf",volcano_plot_1B, width = 5, height = 6)

head(volc1dat)
volc1dat$log2FC <- volc1dat$FC
volc1dat$minus_log10_pvalue <- -log10(volc1dat$p)
setwd("D:/.shortcut-targets-by-id/10J9GdrycmUMcd6u97iLZyH8Psxxr964x/Lipidomics Mouse/Datasets Figures")
write.csv(volc1dat[,c("lipid","log2FC", "minus_log10_pvalue" )], "data_extended_fig4D.csv", row.names = F)



########## data_extended_fig4E ##########



volc1dat$class2 <- "Other"
volc1dat$class2[volc1dat$group_manual=="Ceramides"] <- "CER" 
volc1dat$class2[grepl("Fatty", volc1dat$group_manual)] <- "FA" 
volc1dat$class2[grepl("Diradylglycerols", volc1dat$group_manual)] <- "DG" 
volc1dat$class2[grepl("Triradylcglycerols", volc1dat$group_manual)] <- "TG" 
volc1dat$class2[grepl("Glycerophosphocholines", volc1dat$group_manual)] <- "PC" 
volc1dat$class2[grepl("Glycerophosphoethanolamines", volc1dat$group_manual)] <- "PE" 
#volc1dat$class2[grepl("Glycerophosphoglycerols", volc1dat$group_manual)] <- "PG" 
volc1dat$class2[grepl("Glycerophosphoinositols", volc1dat$group_manual)] <- "PI" 
volc1dat$class2[grepl("Sphingomyelins", volc1dat$group_manual)] <- "SM" 
volc1dat$class2[grepl("St", volc1dat$group_manual)] <- "ST" 
volc1dat$class2[grepl("Neutral glycosphingolipids|Phosphosphingolipids", volc1dat$group_manual)] <- "SM" 
volc1dat$class2[grepl("LPC", volc1dat$Metabolite.name)] <- "LPC" 
volc1dat$class2[grepl("LPE", volc1dat$Metabolite.name)] <- "LPE" 
volc1dat$class2[grepl("HexCer", volc1dat$Metabolite.name)] <- "HEXCER" 

table(volc1dat$class2,exclude=NULL)

volc1dat <- volc1dat[volc1dat$class2!="Other",]

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
                                                      "TG"))

volc1dat_labeled <- volc1dat %>%
  filter(p < 0.01) %>%                                # Retain points with p < 0.01
  group_by(class2) %>%                                # Group by facet
  slice_min(order_by = p, n = 3) %>%                  # Select top 5 based on p-value
  ungroup()

volc1dat_labeled <- as.data.frame(volc1dat_labeled)

volc1dat$fill_class <- ifelse(volc1dat$p > 0.05, "nonsignificant", as.character(volc1dat$class2))

volc1dat$fill_class <- factor(volc1dat$fill_class, levels = c("CER",
                                                              "DG",
                                                              "FA",
                                                              "HEXCER",
                                                              "LPC",
                                                              "LPE",
                                                              "PC",
                                                              "PE",
                                                              "PI",
                                                              "SM",
                                                              "TG",
                                                              "nonsignificant"))

table(volc1dat$fill_class)

volcano_plot_1B_FACET <- ggplot(volc1dat, aes(x = FC, y = -log10(p)  )) +
  geom_vline(xintercept = 0, col="gray70")+
  geom_point(aes(fill = fill_class), shape = 21, alpha = 0.7, col="grey20", size = 2) + # Color points based on significance
  facet_wrap(.~class2)+
  #scale_fill_manual(values = c("#FF0000", "#FF7600", "#FFEB00", "#27FF00", "#00FF4E",
  #                             "#00FFC4", "#00C4FF", "#004EFF", "#2700FF", "#FF00EB" ,"#FF0076" ,"grey60"), drop=T )+
  
  scale_fill_manual(values = c("#FF0000", "#FF7600", "#FFEB00", "#9DFF00", "#27FF00", "#00FF4E",
                               "#00FFC4", "#00C4FF", "#004EFF", "#2700FF", "#FF00EB" ,"lightgray" ), drop=F )+
  
  
  #scale_color_manual(values = c("upregulated" = "red", "downregulated" = "blue", "non-significant" = "black")) + # Assign colors
  labs(x = "Log2 Fold Change", y = "-log10(p-value)",  title = "CLP v Ctrl") +
  scale_x_continuous(breaks=c(-5, 0, 5))+
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank())+
  geom_text_repel(
    data = volc1dat_labeled,
    aes(label = lipid),                  # Replace 'gene_id' with the column containing gene names
    size = 3,
    max.overlaps = Inf,
    min.segment.length = 0,
    box.padding = 0.4
  )


volcano_plot_1B_FACET

ggsave("J:/.shortcut-targets-by-id/10J9GdrycmUMcd6u97iLZyH8Psxxr964x/Lipidomics Mouse/Analyses in progress by H and J/volcano_plot_1B_FACET_labs2.pdf",volcano_plot_1B_FACET, width = 11, height = 9)

volc1dat$log2FC <- volc1dat$FC
volc1dat$minus_log10_pvalue <- -log10(volc1dat$p)
volc1dat$lipid_class <- volc1dat$class2
setwd("D:/.shortcut-targets-by-id/10J9GdrycmUMcd6u97iLZyH8Psxxr964x/Lipidomics Mouse/Datasets Figures")
write.csv(volc1dat[,c("lipid","log2FC", "minus_log10_pvalue","lipid_class")], "data_extended_fig4E.csv", row.names = F)



########## data_extended_fig6A ##########



results_df2$class2 <-  results_df$class
results_df2$class2 <- ifelse(results_df2$class=="LPG"|
                               results_df2$class=="LPI"|
                               results_df2$class=="LPS"|
                               results_df2$class=="PA"|
                               results_df2$class=="PG"|
                               results_df2$class=="PS", "Other", results_df2$class2)


results_df2$class2 <- factor(results_df2$class2, levels = c("CER",
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
                                                            "Other"))


#scatter

plot(results_df2$t_stat_elisa, results_df2$t_stat_mars)

fit <- lm(results_df2$t_stat_mars~results_df2$t_stat_elisa)
abline(fit)
summary(fit)

b  <- coef(fit)["x"]
r_from_b <- b * sd(df$x, na.rm = TRUE) / sd(df$y, na.rm = TRUE)


corrplot <- ggplot(data=results_df2, aes(x=t_stat_elisa, y=t_stat_mars))+
  geom_point( aes( fill =class2), shape = 21, alpha = 0.7, col="grey20", size = 2)+
  theme_bw()+
  geom_smooth(method=lm)+
  xlab("Association with sepsis in mouse (t-statistic)")+
  ylab("Association with sepsis in human (t-statistic)")+
  ggtitle("Correlation, P < 0.0001 ")+
  scale_fill_manual(name="Class",
                    values = c("#FF0000", "#FF7600", "#FFEB00", "#9DFF00", "#27FF00", "#00FF4E",
                               "#00FFC4", "#00C4FF", "#004EFF", "#2700FF", "#9D00FF", "#FF00EB" ,"#FF0076" ), drop=F )


corrplot

ggsave("J:/.shortcut-targets-by-id/10J9GdrycmUMcd6u97iLZyH8Psxxr964x/Lipidomics Mouse/Analyses in progress by H and J/corrplot_new.pdf", corrplot, width = 7, height = 5)

head(results_df2)
results_df2$T_stat_human <- results_df2$t_stat_mars
results_df2$T_stat_mouse <- results_df2$t_stat_elisa
results_df2$lipid_class <- results_df2$class2
setwd("D:/.shortcut-targets-by-id/10J9GdrycmUMcd6u97iLZyH8Psxxr964x/Lipidomics Mouse/Datasets Figures")
write.csv(results_df2[,c("lipid","T_stat_human", "T_stat_mouse","lipid_class"   )], "data_extended_fig6A.csv", row.names = F)


########## data_fig1G ##########



results_df3 <- results_df2[results_df2$class2!="Other",]

results_df3$class2 <- factor(results_df3$class2, levels = c("CER",
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

corrplot2 <- ggplot(data=results_df3, aes(x=t_stat_elisa, y=t_stat_mars))+
  geom_point(shape = 21, alpha = 0.7, col="grey20", size = 2, aes( fill=class2))+
  facet_wrap(.~class2)+
  theme_bw()+
  geom_smooth(method=lm)+
  xlab("Association with sepsis in mouse (t-statistic)")+
  ylab("Association with sepsis in human (t-statistic)")+
  ggtitle("Correlations within Lipid Class ")+
  scale_fill_manual(name="Class",
                    values = c("#FF0000", "#FF7600", "#FFEB00", "#9DFF00", "#27FF00", "#00FF4E",
                               "#00FFC4", "#00C4FF", "#004EFF", "#2700FF", "#9D00FF", "#FF00EB"  ), drop=F )

corrplot2

ggsave("J:/.shortcut-targets-by-id/10J9GdrycmUMcd6u97iLZyH8Psxxr964x/Lipidomics Mouse/Analyses in progress by H and J/corrplot_new_facet_class.pdf", corrplot2, width = 10, height = 8)

#save data for Nature
head(results_df3)
results_df3$T_stat_human <- results_df3$t_stat_mars
results_df3$T_stat_mouse <- results_df3$t_stat_elisa
setwd("D:/.shortcut-targets-by-id/10J9GdrycmUMcd6u97iLZyH8Psxxr964x/Lipidomics Mouse/Datasets Figures")
write.csv(results_df3[,c("lipid","T_stat_human", "T_stat_mouse" )], "data_fig1G.csv", row.names = F)



########## data_fig1G_corr ##########



pvalues <- results_df3  %>%
  group_by(class2) %>%
  summarise(
    p_value = summary(lm(t_stat_mars ~ t_stat_elisa))$coefficients[2, 4],  # Extract p-value for x
    cor = cor(t_stat_mars, t_stat_elisa)  # Extract cor
  )

corrplot2 <- ggplot(data=results_df3, aes(x=t_stat_elisa, y=t_stat_mars))+
  geom_point(shape = 21, alpha = 0.7, col="grey20", size = 2, aes( fill=class2))+
  facet_wrap(.~class2)+
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

ggsave("J:/.shortcut-targets-by-id/10J9GdrycmUMcd6u97iLZyH8Psxxr964x/Lipidomics Mouse/Analyses in progress by H and J/corrplot_new_facet_class_pvals_and_r.pdf", corrplot2, width = 10, height = 8)
ggsave("I:/.shortcut-targets-by-id/10J9GdrycmUMcd6u97iLZyH8Psxxr964x/Lipidomics Mouse/Analyses in progress by H and J/corrplot_new_facet_class_pvals_and_r.pdf", corrplot2, width = 10, height = 8)
#save data Nature
pvalues$lipid_class <- pvalues$class2
setwd("D:/.shortcut-targets-by-id/10J9GdrycmUMcd6u97iLZyH8Psxxr964x/Lipidomics Mouse/Datasets Figures")
write.csv(pvalues[,c("lipid_class", "p_value", "cor")], "data_fig1G_corr.csv", row.names = F)


########## data_fig1F ##########

forvenn2 <- data.frame( MOUSE = results_df2_conc$sig3, 
                        HUMAN = results_df2_conc$sig3h )
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

setwd("D:/Elisa Jena/")
ggsave("mouse_human_venn_concordant.pdf",ve, width = 5, height = 5)

#save data for Nature
results_df2_conc$signif_mouse <- results_df2_conc$sig3
results_df2_conc$signif_human <- results_df2_conc$sig3h
setwd("D:/.shortcut-targets-by-id/10J9GdrycmUMcd6u97iLZyH8Psxxr964x/Lipidomics Mouse/Datasets Figures")
write.csv(results_df2_conc[,c("lipid","signif_mouse", "signif_human" )], "data_fig1F.csv", row.names = F)



########## data_ext_fig9A_fig9B ##########



label3 <- results_df2_2[results_df2_2$class3==T & results_df2_2$t_stat_mars>0 & results_df2_2$t_stat_elisa>0,]

corrplot <- ggplot(data=results_df2_2, aes(x=t_stat_elisa, y=t_stat_mars))+
  geom_hline(yintercept = 0, colour = "grey50")+
  geom_vline(xintercept = 0, colour = "grey50")+
  geom_point( aes( fill =class3), shape = 21, alpha = 0.7, col="grey20", size = 2)+
  # geom_text_repel(aes(label=label2), size=3, box.padding = 2, max.overlaps = Inf,
  #                 segment.color = "lightgrey",  # Light grey lines
  #                 segment.size = 0.3)+
  scale_y_continuous(limits = c(-15,10))+
  scale_x_continuous(limits = c(-15,7))+
  geom_text_repel(aes(label = lipid),
                  data          =  label3,
                  nudge_x       = -15 - label3$t_stat_elisa,
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

setwd("D:/.shortcut-targets-by-id/10J9GdrycmUMcd6u97iLZyH8Psxxr964x/Lipidomics Mouse/Analyses in progress by H and J")
ggsave("corrplot_atgl_vs_MARS_181_from_VENN_updated2.pdf", corrplot, width = 7, height = 5)

#save data for Nature
results_df2_2$signif_mouse <- results_df2_2$sig3
results_df2_2$signif_human <- results_df2_2$sig3h
results_df2_2$overlap <- results_df2_2$class3
results_df2_2$T_stat_human_sepsis <- results_df2_2$t_stat_mars
results_df2_2$T_stat_mouse_sepsis_AdipoQD_D_v_fl_fl <- results_df2_2$t_stat_elisa
setwd("D:/.shortcut-targets-by-id/10J9GdrycmUMcd6u97iLZyH8Psxxr964x/Lipidomics Mouse/Datasets Figures")
write.csv(results_df2_2[!is.na(results_df2_2$t_stat_elisa) ,c("lipid","signif_mouse", "signif_human","overlap", "T_stat_human_sepsis", 
                                                              "T_stat_mouse_sepsis_AdipoQD_D_v_fl_fl")], "data_ext_fig9A_fig9B.csv", row.names = F)


########## data_fig1I, data_ext_fig10_corr, data_ext_fig10_scatter ##########

top25msel <- top_all[top_all$lipid %in% vec,] # vec now not known

corvolc1 <- ggplot(data=top_all)+
  geom_jitter(aes(x=r, y=-log10(p), fill = class2), shape = 21, alpha = 0.7, col="grey20", size = 2,width = 0.02)+
  theme_bw()+
  xlab("Correlation with SOFA (Pearson's r)")+
  ylab("-log10(p)")+
  geom_text_repel(
    data = top25msel,
    aes(label = lipid, x=r, y=-log10(p)),                  
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

#save data for Nature
head(top_all)
top_all$pearson_r <- top_all$r
top_all$p_value <- top_all$p
top_all$lipid_class <- top_all$class2
setwd("D:/.shortcut-targets-by-id/10J9GdrycmUMcd6u97iLZyH8Psxxr964x/Lipidomics Mouse/Datasets Figures")
write.csv(top_all[,c("lipid","pearson_r", "p_value", "lipid_class"  )], "data_fig1I.csv", row.names = F)

write.csv(top_all[top_all$lipid %in% lipid_species_vector,c("lipid","pearson_r", "p_value", "lipid_class"  )], "data_ext_fig10_corr.csv", row.names = F)

df_aggregated2_10 <- df_aggregated2[, c("SOFAtot", lipid_species_vector), drop = FALSE]
write.csv(df_aggregated2_10, "data_ext_fig10_scatter.csv", row.names = F)



########## data_ext_fig7B ##########



top_biomarkers <- volc2dat[volc2dat$ID %in% c("pos_2533","neg_8662","pos_4668"),]

table(volc2dat$class2)

library(ggrepel)

volcano_plot_1F <- ggplot(volc2dat, aes(x = FC, y = -log10(p)  )) +
  geom_point(aes(fill = class2), shape = 21, alpha = 0.7, col="grey20", size = 2) + # Color points based on significance
  
  #scale_fill_manual(values = c("#FF0000", "#FF7600", "#FFEB00", "#27FF00", "#00FF4E",
  #                             "#00FFC4", "#00C4FF", "#004EFF", "#2700FF", "#FF00EB" ,"#FF0076" ,"grey60"), drop=T )+
  
  scale_fill_manual(values = c("#FF0000", "#FF7600", "#FFEB00", "#9DFF00", "#27FF00", "#00FF4E",
                               "#00FFC4", "#00C4FF", "#004EFF", "#2700FF", "#9D00FF", "#FF00EB" ,"#FF0076" ,"grey60"), drop=F )+
  
  
  #scale_color_manual(values = c("upregulated" = "red", "downregulated" = "blue", "non-significant" = "black")) + # Assign colors
  labs(x = "Log2 Fold Change", y = "-log10(p-value)", title = "AdipoQD/D v fl/fl") +
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
                  hjust         = 0, size=2.5
  ) +
  geom_text_repel(aes(label = lipid),
                  data          =  subset(top_biomarkers, FC <0),
                  nudge_x       = -4 - subset(top_biomarkers, FC< 0)$FC,
                  
                  #segment.size  = 0.2,
                  segment.color = "grey50",
                  direction     = "y",
                  hjust         = 0, size=2.5
  ) 


print(volcano_plot_1F)

#save data for Nature
head(volc2dat)
volc2dat2 <- volc2dat
volc2dat2$lipid_name <- volc2dat2$lipid
volc2dat2$neg_log10_p_value <- -log10(volc2dat2$p)
volc2dat2$lipid_class <- volc2dat2$class2
volc2dat2$log_2_foldchange <- volc2dat2$FC
setwd("D:/.shortcut-targets-by-id/10J9GdrycmUMcd6u97iLZyH8Psxxr964x/Lipidomics Mouse/Datasets Figures")
write.csv(volc2dat2[,c("lipid_name","log_2_foldchange", "neg_log10_p_value", "lipid_class"  )], "data_ext_fig7B.csv", row.names = F)



########## data_ext_fig7C_labelled ##########


# first read in dataset from Fig 7B, this is called "volc2dat"


volc2dat$class2 <- "Other"

volc2dat$class2[volc2dat$group_manual=="Ceramides"] <- "CER" 
volc2dat$class2[grepl("Fatty", volc2dat$group_manual)] <- "FA" 
volc2dat$class2[grepl("Diradylglycerols", volc2dat$group_manual)] <- "DG" 
volc2dat$class2[grepl("Triradylcglycerols", volc2dat$group_manual)] <- "TG" 
volc2dat$class2[grepl("Glycerophosphocholines", volc2dat$group_manual)] <- "PC" 
volc2dat$class2[grepl("Glycerophosphoethanolamines", volc2dat$group_manual)] <- "PE" 
#volc2dat$class2[grepl("Glycerophosphoglycerols", volc2dat$group_manual)] <- "PG" 
volc2dat$class2[grepl("Glycerophosphoinositols", volc2dat$group_manual)] <- "PI" 
volc2dat$class2[grepl("Sphingomyelins", volc2dat$group_manual)] <- "SM" 
volc2dat$class2[grepl("St", volc2dat$group_manual)] <- "ST" 
volc2dat$class2[grepl("Neutral glycosphingolipids|Phosphosphingolipids", volc2dat$group_manual)] <- "SM" 
volc2dat$class2[grepl("LPC", volc2dat$Metabolite.name)] <- "LPC" 
volc2dat$class2[grepl("LPE", volc2dat$Metabolite.name)] <- "LPE" 
volc2dat$class2[grepl("HexCer", volc2dat$Metabolite.name)] <- "HEXCER" 

table(volc2dat$class2,exclude=NULL)

volc2dat <- volc2dat[volc2dat$class2!="Other",]

volc2dat$class2 <- factor(volc2dat$class2, levels = c("CER",
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



volc2dat$class_col <- as.character(volc2dat$class2)
volc2dat$class_col[volc2dat$p>0.05] <- "non-sig"

table(volc2dat$class_col)

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
                                                            "non-sig"))


volc2dat$lipid <- volc2dat$metabolites4

volc2dat_labeled <- volc2dat %>%
  filter(p < 0.05) %>%                                # Retain points with p < 0.01
  group_by(class_col) %>%                                # Group by facet
  slice_min(order_by = p, n = 5) %>%                  # Select top 5 based on p-value
  ungroup()

volc2dat_labeled <- as.data.frame(volc2dat_labeled)

volcano_plot_1F_FACET <- ggplot(volc2dat, aes(x = FC, y = -log10(p)  )) +
  geom_vline(xintercept = 0, col="gray70")+
  geom_point(aes(fill = class_col), shape = 21, alpha = 0.7, col="grey20", size = 2) + # Color points based on significance
  facet_wrap(.~class2)+
  #scale_fill_manual(values = c("#FF0000", "#FF7600", "#FFEB00", "#27FF00", "#00FF4E",
  #                             "#00FFC4", "#00C4FF", "#004EFF", "#2700FF", "#FF00EB" ,"#FF0076" ,"grey60"), drop=T )+
  
  scale_fill_manual(values = c("#FF0000", "#FF7600", "#FFEB00", "#9DFF00", "#27FF00", "#00FF4E",
                               "#00FFC4", "#00C4FF", "#004EFF", "#2700FF", "#9D00FF", "#FF00EB" , "grey50"), drop=F )+
  
  
  #scale_color_manual(values = c("upregulated" = "red", "downregulated" = "blue", "non-significant" = "black")) + # Assign colors
  labs(x = "Log2 Fold Change", y = "-log10(p-value)", title = "AdipoQD/D v fl/fl") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank())+
  geom_text_repel(aes(label = `Elisa lipid dataset_converted`),
                  force_pull   = 0,
                  data          = volc2dat_labeled,
                  nudge_y       = 4 + log10(volc2dat_labeled$p),
                  angle        = 90,
                  #nudge_y       = 0.2,
                  #segment.size  = 0.2,
                  segment.color = "grey50",
                  direction     = "x",
                  hjust         = 0, size=2.5
  ) +
  scale_y_continuous(limits = c(0,5))

volcano_plot_1F_FACET

#save data for Nature
head(volc2dat_labeled)
volc2dat3 <- volc2dat_labeled
volc2dat3$lipid_name <- volc2dat3$lipid
volc2dat3$neg_log10_p_value <- -log10(volc2dat3$p)
volc2dat3$lipid_class <- volc2dat3$class_col
volc2dat3$log_2_foldchange <- volc2dat3$FC
setwd("D:/.shortcut-targets-by-id/10J9GdrycmUMcd6u97iLZyH8Psxxr964x/Lipidomics Mouse/Datasets Figures")
write.csv(volc2dat3[,c("lipid_name","log_2_foldchange", "neg_log10_p_value", "lipid_class"  )], "data_ext_fig7C_labelled.csv", row.names = F)





