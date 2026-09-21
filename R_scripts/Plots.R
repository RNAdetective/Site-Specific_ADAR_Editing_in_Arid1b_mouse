
Packages <- c('dplyr','tidyr','hrbrthemes','tidyverse','ggpubr','reshape2','circlize','ggsci',
              'ggplot2','ggstatsplot','ggsignif','cowplot','rlang','ComplexHeatmap','RColorBrewer',
              'Cairo','patchwork','ggVennDiagram','clusterProfiler','enrichplot','org.Mm.eg.db')
invisible(lapply(Packages,library, character.only=TRUE))


#######load inputfiles
setwd("C:/Users/labuser/Desktop/ASD-miRNA/NewPlots/FigureInputFiles")

filenames <- gsub("\\.csv$","", list.files(pattern="\\.csv$"))

for(i in filenames){
  assign(i, read.csv(paste(i, ".csv", sep="")))
}
###########VennDiagrams and heatmap########################################

# List of Coordinates
CoordinateList <- list(HT_drug = HTdrug_annotated$Coordinate, HT_veh = HTvehicle_annotated$Coordinate, 
                       WT_drug=WTdrug_annotated$Coordinate, WT_veh=WTvehicle_annotated$Coordinate)
list_count <- sizes <- sapply(CoordinateList, length)
list_names <- paste0(names(CoordinateList), "(", list_count, ")")

# Venn diagram
########Overall Overlapping targets#########
VennDiagram <- ggVennDiagram(CoordinateList, color = "black", lwd = 0.8, lty = 1, label = "both", category.names = list_names, set_size = 6) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") + theme(legend.position = "none") 
VennDiagram <- VennDiagram + scale_x_continuous(expand = expansion(mult = .3))


##############Coding regions and non coding and non-coding exonic, non_coding intron targets######################
######nonCoding Regulatory Regions###############
HTdrugNonCoding <- filter(HTdrug_annotated,Func.wgEncodeGencodeBasicVM16 == "UTR3" |Func.wgEncodeGencodeBasicVM16 == "UTR5" | Func.wgEncodeGencodeBasicVM16 == "intronic"| Func.wgEncodeGencodeBasicVM16 == "downstream" |Func.wgEncodeGencodeBasicVM16 == "intergenic")
HTvehNonCoding <- filter(HTvehicle_annotated,Func.wgEncodeGencodeBasicVM16 == "UTR3" |Func.wgEncodeGencodeBasicVM16 == "UTR5"| Func.wgEncodeGencodeBasicVM16 == "intronic"| Func.wgEncodeGencodeBasicVM16 == "downstream" |Func.wgEncodeGencodeBasicVM16 == "intergenic")
WTdrugNonCoding <- filter(WTdrug_annotated,Func.wgEncodeGencodeBasicVM16 == "UTR3" |Func.wgEncodeGencodeBasicVM16 == "UTR5" | Func.wgEncodeGencodeBasicVM16 == "intronic"| Func.wgEncodeGencodeBasicVM16 == "downstream" |Func.wgEncodeGencodeBasicVM16 == "intergenic")
WTvehNonCoding <- filter(WTvehicle_annotated,Func.wgEncodeGencodeBasicVM16 == "UTR3" |Func.wgEncodeGencodeBasicVM16 == "UTR5" | Func.wgEncodeGencodeBasicVM16 == "intronic"| Func.wgEncodeGencodeBasicVM16 == "downstream" |Func.wgEncodeGencodeBasicVM16 == "intergenic")

noncodingRegionList <- list(HT_drug = HTdrugNonCoding$Coordinate, HT_veh = HTvehNonCoding$Coordinate, 
                   WT_drug=WTdrugNonCoding$Coordinate, WT_veh=WTvehNonCoding$Coordinate)


list_count <- sizes <- sapply(noncodingRegionList, length)
list_names <- paste0(names(noncodingRegionList), "(", list_count, ")")



NonCodingVennDiagram <- ggVennDiagram(noncodingRegionList, color = "black", lwd = 0.8, lty = 1, category.names = list_names, set_size = 6) + 
  scale_fill_gradient(low = "lightpink", high = "hotpink") + theme(legend.position = "none") 

NonCodingVennDiagram<- NonCodingVennDiagram  + scale_x_continuous(expand = expansion(mult = .3))
##############Coding################
HTdrugCoding <- filter(HTdrug_annotated,Func.wgEncodeGencodeBasicVM16 == "exonic" )
HTvehCoding <- filter(HTvehicle_annotated,Func.wgEncodeGencodeBasicVM16 == "exonic")
WTdrugCoding <- filter(WTdrug_annotated,Func.wgEncodeGencodeBasicVM16 == "exonic" )
WTvehCoding <- filter(WTvehicle_annotated,Func.wgEncodeGencodeBasicVM16 == "exonic" )

CodingRegionList <- list(HT_drug = HTdrugCoding$Coordinate, HT_veh = HTvehCoding$Coordinate, 
                         WT_drug=WTdrugCoding$Coordinate, WT_veh=WTvehCoding$Coordinate)

list_count <- sizes <- sapply(CodingRegionList, length)
list_names <- paste0(names(CodingRegionList), "(", list_count, ")")



CodingVennDiagram <- ggVennDiagram(CodingRegionList, color = "black", lwd = 0.8, lty = 1, category.names = list_names, set_size = 6) + 
  scale_fill_gradient(low = "lightgreen", high = "green4") + theme(legend.position = "none") 

CodingVennDiagram<- CodingVennDiagram + scale_x_continuous(expand = expansion(mult = .3))
#########################ncRNA intronic and exonic Regions#############

HTdrugncRNA <- filter(HTdrug_annotated,Func.wgEncodeGencodeBasicVM16 == "ncRNA_exonic" |Func.wgEncodeGencodeBasicVM16 == "ncRNA_intronic" )
HTvehncRNA <- filter(HTvehicle_annotated,Func.wgEncodeGencodeBasicVM16 == "ncRNA_exonic"| Func.wgEncodeGencodeBasicVM16 == "ncRNA_intronic")
WTdrugncRNA <- filter(WTdrug_annotated,Func.wgEncodeGencodeBasicVM16 == "ncRNA_exonic"| Func.wgEncodeGencodeBasicVM16 == "ncRNA_intronic" )
WTvehncRNA <- filter(WTvehicle_annotated,Func.wgEncodeGencodeBasicVM16 == "ncRNA_exonic"| Func.wgEncodeGencodeBasicVM16 == "ncRNA_intronic" )

ncRNAList <- list(HT_drug = HTdrugncRNA$Coordinate, HT_veh = HTvehncRNA$Coordinate, 
                  WT_drug=WTdrugncRNA$Coordinate, WT_veh=WTvehncRNA$Coordinate)


list_count <- sizes <- sapply(ncRNAList, length)
list_names <- paste0(names(ncRNAList), "(", list_count, ")")



ncRNAVennDiagram <- ggVennDiagram(ncRNAList, color = "black", lwd = 0.8, lty = 1, category.names = list_names, set_size = 6) + 
  scale_fill_gradient(low = "orange", high = "darkorange") + theme(legend.position = "none") 

ncRNAVennDiagram <-ncRNAVennDiagram + scale_x_continuous(expand = expansion(mult = .3))

Venndia  <- combine_plots(
  list(VennDiagram, NonCodingVennDiagram, CodingVennDiagram, ncRNAVennDiagram),
  plotgrid.args = list(nrow = 2)) & theme(plot.tag = element_text(face="bold"))


ggsave(filename = "VennDiagram.pdf",
       plot = Venndia,
       scale = 0.8,
       width = 22,
       height = 12,
       device = cairo_pdf,
       units = "in",
       dpi=1200,
       bg="white")

################HeatMap#####################
WholeSampleEditedSites <- rbind(HTdrug_annotated, HTvehicle_annotated, WTdrug_annotated, WTvehicle_annotated)


WholeSampleEditedSites <- data.frame(Condition=WholeSampleEditedSites$Condition, Mean.Editing.Rate= WholeSampleEditedSites$Mean.editing.rate,
                                     Coordinate= paste0(WholeSampleEditedSites$Coordinate,":",WholeSampleEditedSites$Gene.wgEncodeGencodeBasicVM16 ))


EditedSites <- WholeSampleEditedSites%>% pivot_wider(
  names_from = Condition, values_from = Mean.Editing.Rate
)


col_fun = colorRamp2(c(0,0.1,0.3,0.6,0.9), c("white","lightgreen", "pink", "skyblue", "red"))



m <- as.matrix(EditedSites[,-1])
m <- m%>% replace(is.na(.),0)
rownames(m)<- EditedSites$Coordinate


pdf(file="C:/Users/labuser/Desktop/ASD-miRNA/NewPlots/HeatMap.pdf")


HeatMap <- ComplexHeatmap::pheatmap(m, border_color = "black", cellwidth =6,fontsize_row = 8,fontsize_col = 6,
                                    cluster_cols = F, cluster_rows = T,col=col_fun,
                                    cellheight =15,  heatmap_legend_param = list(title="Editing Rate", title_Position="centertop", legend_direction="horizontal"))

draw(HeatMap)
dev.off()

###########################SharedSitesPlot##########################

Drug_HT_WT_sharedSites$Dot_color <- ifelse(Drug_HT_WT_sharedSites$MeanDiff > 0,">HT", ">WT")


p1 <- ggplot(Drug_HT_WT_sharedSites, aes(x=MeanDiff, y =Gene.wgEncodeGencodeBasicVM16.x, color= Dot_color)) +
  geom_vline(aes(xintercept = 0), size = .25, linetype = "dashed") +
  geom_errorbarh(aes(xmax =Upper_95_CI, xmin = Lower_95_CI), size = .5, height = .2, color = "gray50") +
  geom_point(size = 2) + scale_color_brewer(palette = "Set1")+
  theme_minimal() +
  theme(axis.text.x = element_text(vjust = 0.5, hjust=0.5),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold"))+ labs(color="Mean Editing Rate")+
  ylab("Gene") + 
  xlab("Mean_Editing_HT - Mean_Editing_WT") + labs(title ="Fluoxetine" ) + scale_x_continuous(expand = expansion(add=0.05))

p1

Vehicle_HT_WT_sharedSites$Dot_color <- ifelse(Vehicle_HT_WT_sharedSites$MeanDiff > 0,">HT", ">WT")


p2 <- ggplot(Vehicle_HT_WT_sharedSites, aes(x=MeanDiff, y = Gene.wgEncodeGencodeBasicVM16.x, color= Dot_color)) +
  geom_vline(aes(xintercept = 0), size = .25, linetype = "dashed") +
  geom_errorbarh(aes(xmax = Upper_95_CI, xmin = Lower_95_CI), size = .5, height = .2, color = "gray50") +
  geom_point(size = 2) + scale_color_brewer(palette = "Set1")+
  theme_minimal() +
  theme(axis.text.x = element_text(vjust = 0.5, hjust=0.5),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold"))+ labs(color="Mean Editing Rate")+
  ylab("Gene") + 
  xlab("Mean_Editing_HT - Mean_Editing_WT") + labs(title ="Vehicle" ) + scale_x_continuous(expand = expansion(mult = c(-0.15,0.15)))


p2

Shared_sites_plot <- combine_plots(list(p1,p2), plotgrid.args = list(nrow=1))
Shared_sites_plot

ggsave(filename = "SharedSites.pdf",
       plot = Shared_sites_plot,
       scale = 0.8,
       width = 11,
       height = 8,
       device = cairo_pdf,
       units = "in",
       dpi=1200,
       bg="white")

##############GO plots for Shared Sites########################
ego1 <- enrichGO(gene =Drug_HT_WT_sharedSites$Gene.wgEncodeGencodeBasicVM16.x,
                 OrgDb = org.Mm.eg.db,
                 keyType = 'SYMBOL',
                 ont = "BP",
                 pAdjustMethod = "BH", 
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05,
                 maxGSSize = 50,
                 readable = TRUE) 

ego2 <- enrichGO(gene =Vehicle_HT_WT_sharedSites$Gene.wgEncodeGencodeBasicVM16.x,
                 OrgDb = org.Mm.eg.db,
                 keyType = 'SYMBOL',
                 ont = "BP",
                 pAdjustMethod = "BH", 
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05, 
                 maxGSSize = 50,
                 readable = TRUE) 


p3<- cnetplot(ego1,node_label="all",
              color_edge ='black' , color_item = 'steelblue',color_category = 'firebrick') +
  theme(legend.position = "none")

p4<- cnetplot(ego2,node_label="all",color_edge ='black' , color_item = 'steelblue',color_category = 'firebrick')+
  theme(legend.position = "none")

Shared_sites_pathway_plot <- combine_plots(list(p3,p4), plotgrid.args = list(nrow=1))
Shared_sites_pathway_plot
ggsave(filename = "SharedSitespathway.pdf",
       plot = Shared_sites_pathway_plot,
       scale = 0.8,
       width = 10,
       height = 6,
       device = cairo_pdf,
       units = "in",
       dpi=1200,
       bg="white")
##############Unique Sites###########################
All_Unique_Sites_data <- rbind(HT_drug_uniqueSites,HT_vehicle_uniqueSites, WT_drug_uniqueSites, WT_vehicle_uniqueSites)
violinPlot <- ggviolin(All_Unique_Sites_data, x="Condition", y="Mean.editing.rate", fill = "Condition", repel = TRUE,
                       draw_quantiles = 0.5, palette = c("#00AFBA","#E7B800", "#FC4E07", "#EF1"), 
                       add="jitter", label = "Gene.wgEncodeGencodeBasicVM16", 
                       font.label = list(size=6,color="black"), error.bar=TRUE)+
  theme(axis.title = element_text(face = "bold")) + xlab("Condition") + ylab("Mean.Editing.Rate")
 
violinPlot <-violinPlot + theme(legend.position = "none")

violinPlot

ggsave(filename = "uniquesitesViolinPlot.pdf",
       plot =violinPlot ,
       scale = 0.8,
       width = 6,
       height = 6,
       device = cairo_pdf,
       units = "in",
       dpi=1200,
       bg="white")
##########################miRNA-mRNA chordDiagram#############################
#########Prepare input files
veh_edited <- veh_sharedEdited[, c("Gene","miRNA","Folding_energy")]

veh_unedited <- veh_sharedUnedited[, c("Gene","miRNA","Folding_energy")]

Drug_edited <- Drug_sharedEdited[, c("Gene","miRNA","Folding_energy")]

Drug_unedited <- Drug_SharedUnedited[, c("Gene","miRNA","Folding_energy")]
breaks <- c(-1,-5,-10,-15,-20,-25,-30)
coul <- colorRamp2(breaks = breaks, colors = c("white","blue", "yellow","green","pink","purple", "red"))
##############veh edited and unedited
###set number of cols for plot
par(mfrow = c(1, 1), mar = c(0, 0, 0, 0))


chordDiagram(veh_edited,  col = coul , annotationTrack = "grid",fontsize(100),
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(veh_edited))))))


circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA) # here set bg.border to NA is important

chordDiagram(veh_unedited,  col = coul , annotationTrack = "grid",fontsize(100),
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(veh_unedited))))))

circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA) 
##########Drug edited unedited############
par(mfrow = c(1, 1), mar = c(0, 0, 0, 0))

chordDiagram(Drug_edited,  col = coul , annotationTrack = "grid",fontsize(100),
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(Drug_edited))))))

circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA) 

chordDiagram(Drug_unedited,  col = coul , annotationTrack = "grid",fontsize(100),
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(Drug_unedited))))))

circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA) 

lgd_FoldingEnergy = Legend(at = c(breaks), col_fun =coul , 
                           title_position = "topcenter", title = "Folding energy", direction = "horizontal")


draw(
  lgd_FoldingEnergy,
  x = unit(0.5, "npc"),
  y = unit(0.02, "npc"),
  just = c("center", "bottom")
)

################GO plots for DEGs and Differentially edited genes################
library(forcats)
set.seed(123)

egoDEgs <- enrichGO(gene =DEGsvsDEdi$Differntially_Expressed_Genes_drug_veh,
                 OrgDb = org.Mm.eg.db,
                 keyType = 'SYMBOL',
                 ont = "ALL",
                 pAdjustMethod = "BH", 
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05, 
                 maxGSSize = 50,
                 readable = TRUE)
########To view GO results##############
egoDEgs <- egoDEgs@result
write.csv(egoDEgs, "egoDEGs.csv")

P1 <- ggplot(data = egoDEgs, aes(x=fct_rev(fct_reorder(Description, ONTOLOGY)),y= -log10(p.adjust), fill  = ONTOLOGY, group = ONTOLOGY)) + 
  geom_col(width=0.8, position = "dodge")+ labs(x=NULL) + 
  scale_fill_manual(values = c("#9e66ab","#CCEDB4","#41B7C4"))
P1 <- P1+ theme_classic()  

P1 <- P1+  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+ theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1))
P1


ggsave(
  filename = "DEGsGoplot.pdf",
  plot = P1,
  scale = 0.8,
  device = cairo_pdf,
  width = 6,
  height = 6,
  units = "in",
  bg="white",
  dpi = 1200
)

egoDEedited <- enrichGO(gene =DEGsvsDEdi$DifferentiallyEditedwithmiRNA,
                    OrgDb = org.Mm.eg.db,
                    keyType = 'SYMBOL',
                    ont = "ALL",
                    pAdjustMethod = "BH", 
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05, 
                    maxGSSize = 50,
                    readable = TRUE)




P2 <- barplot(
  egoDEedited,
  x = "pvalue",
  showCategory = 15,
  split = "ONTOLOGY"
) +
  aes(fill = ONTOLOGY) +
  scale_fill_brewer(palette = "Blues") +
  enrichplot::autofacet(by = "row", scales = "free")


P2 <- P2 + scale_y_discrete()+ theme(legend.position = "none")+ labs(x="pvalue") 

P2

###########To view GO results#################
egoDEedited <- egoDEedited@result
write.csv(egoDEedited, "egoDEedited.csv")

ggsave(
  filename = "DEditedGoplot.pdf",
  plot = P2,
  scale = 0.8,
  device = cairo_pdf,
  width =14,
  height = 14,
  units = "in",
  bg="white",
  dpi = 1200
)
#######################ADAR Gene Plot################################
p1 <- ggbetweenstats(
  data = ADAR1_N_counts,
  x = Condition,
  y = Normalized_counts,
  type = "nonparametric",
  results.subtitle = TRUE,
  pairwise.display = "all",
  pairwise.comparison= TRUE,
  p.adjust.method = "BH",
  , pairwise.annotation = "p.value", 
  conf.level = 0.95,
  title = "Adar1",
  package = "ggsci",
  palette = "nrc_npg")
p1

p2 <- ggbetweenstats(
  data = ADAR2_N_counts,
  x = Condition,
  y = Normalized_counts,
  type = "nonparametric",
  results.subtitle = TRUE, 
  p.adjust.method = "BH",
  pairwise.comparison=TRUE,
  pairwise.display = "all",
  pairwise.annotation = "p.value",
  conf.level = 0.95,
  title = "Adar2",
  package = "ggsci",
  palette = "nrc_npg"
)
p2

p3 <- ggbetweenstats(
  data = ADAR3_N_counts,
  x = Condition,
  y = Normalized_counts,
  type = "nonparametric",
  results.subtitle = TRUE, 
  p.adjust.method = "BH",
  pairwise.comparison=TRUE,
  pairwise.display = "all",
  pairwise.annotation = "p.value",
  conf.level = 0.95,
  title = "Adar3",
  package = "ggsci",
  palette = "nrc_npg"
)
p3

ADAR_gene_plot <- combine_plots(
  list(p1, p2,p3),
  plotgrid.args = list(nrow = 1)
  
  
)
ADAR_gene_plot

ggsave(
  filename = "ADARGenePlot.pdf",
  plot = ADAR_gene_plot,
  scale = 0.6,
  device = cairo_pdf,
  width =25,
  height = 15,
  units = "in",
  bg="white",
  dpi = 1200
)


##################Supplemental Figures######################
##########Global Editing Rate######################################
WholeSampleStatistics <- rbind(HTdrug_annotated, HTvehicle_annotated, WTdrug_annotated, WTvehicle_annotated)


GlobalEditingRate <- ggbetweenstats(
  data = WholeSampleStatistics,
  x = Condition,
  y = Mean.editing.rate,
  type = "nonparametric",
  p.adjust.method = "BH",
  pairwise.comparison=TRUE,
  pairwise.display = "all",
  pairwise.annotation = "p.value",
  results.subtitle = FALSE, 
  title = "",
  package = "ggsci",
  palette = "nrc_npg"
)
GlobalEditingRate

ggsave(
  filename = "GlobalEditingRate.pdf",
  plot =GlobalEditingRate ,
  scale = 0.6,
  device = cairo_pdf,
  width =20,
  height = 15,
  units = "in",
  bg="white",
  dpi = 1200
)

#################Barplots for unique GO####################
ego3 <- enrichGO(gene =HT_drug_uniqueSites$Gene.wgEncodeGencodeBasicVM16,
                                   OrgDb = org.Mm.eg.db,
                                   keyType = 'SYMBOL',
                                   ont = "ALL",
                                   pAdjustMethod = "BH",
                                   pvalueCutoff = 0.05,
                                   qvalueCutoff = 0.05, 
                                   maxGSSize = 50,
                                   readable = TRUE) 




ego4 <- enrichGO(gene =HT_vehicle_uniqueSites$Gene.wgEncodeGencodeBasicVM16,
                                     OrgDb = org.Mm.eg.db,
                                   keyType = 'SYMBOL',
                                    ont = "ALL",
                                    pAdjustMethod = "BH",
                                   pvalueCutoff = 0.05,
                                   qvalueCutoff = 0.05, 
                                   maxGSSize = 50,
                                   readable = TRUE) 
 ego5 <- enrichGO(gene =WT_drug_uniqueSites$Gene.wgEncodeGencodeBasicVM16,
                                    OrgDb = org.Mm.eg.db,
                                    keyType = 'SYMBOL',
                                    ont = "ALL",
                                    pAdjustMethod = "BH", 
                                    pvalueCutoff = 0.05,
                                    qvalueCutoff = 0.05, 
                                    maxGSSize = 50,
                                    readable = TRUE) 
 ego6 <- enrichGO(gene =WT_vehicle_uniqueSites$Gene.wgEncodeGencodeBasicVM16,
                                   OrgDb = org.Mm.eg.db,
                                    keyType = 'SYMBOL',
                                    ont = "ALL",
                                     pAdjustMethod = "BH", 
                                    pvalueCutoff = 0.05,
                                    qvalueCutoff = 0.05, 
                                    maxGSSize = 50,
                                    readable = TRUE) 
######################Plot##################################################
 
 P3 <- barplot(
   ego3,
   x = "pvalue",
   showCategory = 15,
   split = "ONTOLOGY"
 ) +
   aes(fill = ONTOLOGY) +
   scale_fill_brewer(palette = "Blues") +
   enrichplot::autofacet(by = "row", scales = "free")
 
 
 P3 <- P3 + scale_y_discrete()+ theme(legend.position = "none")+ labs(x="pvalue")
 
 P3
 
 P4 <- barplot(
   ego4,
   x = "pvalue",
   showCategory = 15,
   split = "ONTOLOGY"
 ) +
   aes(fill = ONTOLOGY) +
   scale_fill_brewer(palette = "Greens") +
   enrichplot::autofacet(by = "row", scales = "free")
 
 
 P4 <- P4 + scale_y_discrete()+ theme(legend.position = "none")+ labs(x="pvalue")
 
 P4 
 
 P5 <- barplot(
   ego5,
   x = "pvalue",
   showCategory = 15,
   split = "ONTOLOGY"
 ) +
   aes(fill = ONTOLOGY) +
   scale_fill_brewer(palette = "Reds") +
   enrichplot::autofacet(by = "row", scales = "free")
 
 
 P5 <- P5 + scale_y_discrete()+ theme(legend.position = "none")+ labs(x="pvalue")
 
 P5 

 
 P6 <- barplot(
   ego6,
   x = "pvalue",
   showCategory = 15,
   split = "ONTOLOGY"
 ) +
   aes(fill = ONTOLOGY) +
   scale_fill_brewer(palette = "OrRd") +
   enrichplot::autofacet(by = "row", scales = "free")
 
 
 P6 <- P6 + scale_y_discrete()+ theme(legend.position = "none")+ labs(x="pvalue")
 
 P6 
 UniqueSitesGOplot <- combine_plots(
   list(P3,P4,P5,P6),
   plotgrid.args = list(ncol = 2), annotation.args = list(tag_levels="A"))
   
UniqueSitesGOplot
ggsave(
  filename = "UniqueSitesGOplot.tiff",
  plot =UniqueSitesGOplot,
  scale = 1,
  device = "tiff",
  width =22,
  height = 15,
  units = "in",
  bg="white",
  dpi = 1200
) 
 
 ##################To view Go results#########################
 HTdrugUniqueGO <- ego3@result
 HTvehUniqueGO <- ego4@result
 WTdrugUniqueGO <- ego5@result
 WTvehUniqueGO <- ego6@result
############################ADAR1 Isoforms######################
 ADARp110 <-filter(ADARTransCounts,Transcript == "Adarp110")
 ADARp150 <-filter(ADARTransCounts,Transcript == "Adarp150")
 
 p1 <- ggbetweenstats(
     data = ADARp150,
   x = Condition,
      y =Counts,
      type = "nonparametric",
      results.subtitle = FALSE,
   pairwise.display = "all",
   pairwise.comparison= TRUE,
   p.adjust.method = "BH",
   pairwise.annotation = "p.value", 
   conf.level = 0.95,
   title = "Adarp150",
   point.args = list(alpha = 0.8, size =4, stroke = 0)) +ggplot2::scale_color_manual(values = c("darkred", "#E69F00", "#56B4E9", "purple"))
 
 
 
 p1 <- p1 + theme(axis.title.y.right = element_blank(), 
          axis.text.y.right = element_blank(), 
          axis.ticks.y.right = element_blank(),
          axis.text.x = element_text(vjust = 0, hjust = 0.2)) 
 
 p1
 
 p2 <- ggbetweenstats(
    data = ADARp110,
    x = Condition,
    y = Counts,
   type = "nonparametric",
   results.subtitle = FALSE, 
   p.adjust.method = "BH",
    pairwise.comparison=TRUE,
   pairwise.display = "all",
    pairwise.annotation = "p.value",
    conf.level = 0.95,
   title = "Adarp110",
   point.args = list(alpha = 0.8, size =4, stroke = 0)) +ggplot2::scale_color_manual(values = c("darkred", "#E69F00", "#56B4E9", "purple"))
 
 
 
 p2 <- p2 +theme(axis.title.y.right = element_blank(), 
         axis.text.y.right = element_blank(), 
         axis.ticks.y.right = element_blank(),
         axis.text.x = element_text(vjust = 0, hjust = 0.2)) 
 p2
 
 ADAR1isoformPlots <- combine_plots(list(p1,p2))
 
 ADAR1isoformPlots

 ggsave(
   filename = "ADAR1isoformPlots.pdf",
   plot =ADAR1isoformPlots,
   scale = 0.8,
   device = cairo_pdf(),
   width =15,
   height = 12,
   units = "in",
   bg="white",
   dpi = 1200
 ) 

#############################################
 
 
 
 
 
 
 