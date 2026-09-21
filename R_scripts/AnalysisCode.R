########Libraries#####################
library(dplyr)
library(readr)
library(ggplot2)
library(stringr)
library(tidyr)
library(tidyverse)
library(Biostrings)
library(devtools)
#==============================Import data===========================================================
setwd("~/Desktop/ARID1B_ASD/frequency_counts")
filenames <- gsub("\\.csv$","",list.files(pattern = "\\.csv$"))

for (i in filenames) { 
  assign(i,read.csv(paste(i,".csv", sep="")))
  
}

HTveh1<- data.frame(SRR19308802All, Run= rep("SRR19308802",nrow(SRR19308802All)), Condition = rep("HTveh",nrow(SRR19308802All)))
names(HTveh1) <- c("CHROM","POS","REF","TOTAL","A","C","G","T","Run","Condition")
HTveh2<- data.frame(SRR19308803All, Run= rep("SRR19308803",nrow(SRR19308803All)), Condition = rep("HTveh",nrow(SRR19308803All)))
names(HTveh2) <- c("CHROM","POS","REF","TOTAL","A","C","G","T","Run","Condition")
HTveh3<- data.frame(SRR19308804All, Run= rep("SRR19308804",nrow(SRR19308804All)), Condition = rep("HTveh",nrow(SRR19308804All)))
names(HTveh3) <- c("CHROM","POS","REF","TOTAL","A","C","G","T","Run","Condition")
HTveh4<- data.frame(SRR19308805All, Run= rep("SRR19308805",nrow(SRR19308805All)), Condition = rep("HTveh",nrow(SRR19308805All)))
names(HTveh4) <- c("CHROM","POS","REF","TOTAL","A","C","G","T","Run","Condition")
#Combine HTvehicle data
HTvehicle <- rbind(HTveh1,HTveh2,HTveh3,HTveh4)
#HTdrug
HTdrug1<- data.frame(SRR19308806All, Run= rep("SRR19308806",nrow(SRR19308806All)), Condition = rep("HTdrug",nrow(SRR19308806All)))
names(HTdrug1) <- c("CHROM","POS","REF","TOTAL","A","C","G","T","Run","Condition")
HTdrug2<- data.frame(SRR19308807All, Run= rep("SRR19308807",nrow(SRR19308807All)), Condition = rep("HTdrug",nrow(SRR19308807All)))
names(HTdrug2) <- c("CHROM","POS","REF","TOTAL","A","C","G","T","Run","Condition")
HTdrug3<- data.frame(SRR19308808All, Run= rep("SRR19308808",nrow(SRR19308808All)), Condition = rep("HTdrug",nrow(SRR19308808All)))
names(HTdrug3) <- c("CHROM","POS","REF","TOTAL","A","C","G","T","Run","Condition")
HTdrug4<- data.frame(SRR19308809All, Run= rep("SRR19308809",nrow(SRR19308809All)), Condition = rep("HTdrug",nrow(SRR19308809All)))
names(HTdrug4) <- c("CHROM","POS","REF","TOTAL","A","C","G","T","Run","Condition")

#Combine HTdrug data
HTdrug <- rbind(HTdrug1,HTdrug2,HTdrug3,HTdrug4)
# WT veh
WTveh1<- data.frame(SRR19308794All, Run= rep("SRR19308794",nrow(SRR19308794All)), Condition = rep("WTveh",nrow(SRR19308794All)))
names(WTveh1) <- c("CHROM","POS","REF","TOTAL","A","C","G","T","Run","Condition")
WTveh2<- data.frame(SRR19308795All, Run= rep("SRR19308795",nrow(SRR19308795All)), Condition = rep("WTveh",nrow(SRR19308795All)))
names(WTveh2) <- c("CHROM","POS","REF","TOTAL","A","C","G","T","Run","Condition")
WTveh3<- data.frame(SRR19308796All, Run= rep("SRR19308796",nrow(SRR19308796All)), Condition = rep("WTveh",nrow(SRR19308796All)))
names(WTveh3) <- c("CHROM","POS","REF","TOTAL","A","C","G","T","Run","Condition")
WTveh4<- data.frame(SRR19308797All, Run= rep("SRR19308797",nrow(SRR19308797All)), Condition = rep("WTveh",nrow(SRR19308797All)))
names(WTveh4) <- c("CHROM","POS","REF","TOTAL","A","C","G","T","Run","Condition")
#Combine WTvehicle data
WTvehicle <- rbind(WTveh1,WTveh2,WTveh3,WTveh4)
#WTdrug

WTdrug1<- data.frame(SRR19308798All, Run= rep("SRR19308798",nrow(SRR19308798All)), Condition = rep("WTdrug",nrow(SRR19308798All)))
names(WTdrug1) <- c("CHROM","POS","REF","TOTAL","A","C","G","T","Run","Condition")
WTdrug2<- data.frame(SRR19308799All, Run= rep("SRR19308799",nrow(SRR19308799All)), Condition = rep("WTdrug",nrow(SRR19308799All)))
names(WTdrug2) <- c("CHROM","POS","REF","TOTAL","A","C","G","T","Run","Condition")
WTdrug3<- data.frame(SRR19308800All, Run= rep("SRR19308800",nrow(SRR19308800All)), Condition = rep("WTdrug",nrow(SRR19308800All)))
names(WTdrug3) <- c("CHROM","POS","REF","TOTAL","A","C","G","T","Run","Condition")
WTdrug4<- data.frame(SRR19308801All, Run= rep("SRR19308801",nrow(SRR19308801All)), Condition = rep("WTdrug",nrow(SRR19308801All)))
names(WTdrug4) <- c("CHROM","POS","REF","TOTAL","A","C","G","T","Run","Condition")
#Combine WTdrug data
WTdrug <- rbind(WTdrug1,WTdrug2,WTdrug3,WTdrug4)
###########Filter to keep only A and T in REF allele############################
HTvehicle <- filter(HTvehicle, REF =="A" | REF == "T")
HTdrug <- filter(HTdrug, REF=="A"| REF=="T")
WTvehicle <- filter(WTvehicle, REF=="A"| REF=="T")
WTdrug <- filter(WTdrug, REF=="A"| REF=="T")

#Calculate editing frequency for each site: distinguishes between A->G and T->C
for (n in 1:nrow(HTvehicle)){
  if (HTvehicle$REF[n] == "A"){
    HTvehicle$Editing.rate[n] = HTvehicle$G[n]/HTvehicle$TOTAL[n]
  }
  else{
    HTvehicle$Editing.rate[n] = HTvehicle$C[n]/HTvehicle$TOTAL[n]
  }
} 

for (n in 1:nrow(HTdrug)){
  if (HTdrug$REF[n] == "A"){
    HTdrug$Editing.rate[n] = HTdrug$G[n]/HTdrug$TOTAL[n]
  }
  else{
    HTdrug$Editing.rate[n] = HTdrug$C[n]/HTdrug$TOTAL[n]
  }
} 


for (n in 1:nrow(WTvehicle)){
  if (WTvehicle$REF[n] == "A"){
    WTvehicle$Editing.rate[n] = WTvehicle$G[n]/WTvehicle$TOTAL[n]
  }
  else{
    WTvehicle$Editing.rate[n] = WTvehicle$C[n]/WTvehicle$TOTAL[n]
  }
} 

for (n in 1:nrow(WTdrug)){
  if (WTdrug$REF[n] == "A"){
    WTdrug$Editing.rate[n] = WTdrug$G[n]/WTdrug$TOTAL[n]
  }
  else{
    WTdrug$Editing.rate[n] = WTdrug$C[n]/WTdrug$TOTAL[n]
  }
} 
############save samples Tables with no filters#####################
FileList <- list(HTvehicle=HTvehicle, HTdrug=HTdrug, WTvehicle=WTvehicle,WTdrug=WTdrug)
setwd("~/Desktop/ARID1B_ASD/WholeSampleFiles/CountsWithNoFilters")
for (i in names(FileList)) {write.csv(FileList[[i]], paste0(i,".csv"))
  
}

#===================Filter ADAR sites===========================================

HTvehicle <- filter(HTvehicle, Editing.rate <= .99 & Editing.rate >= 0.01)
HTvehicle <- filter(HTvehicle, Editing.rate <= .49 | Editing.rate  >= 0.51)

HTdrug <- filter(HTdrug, Editing.rate <= .99 & Editing.rate >= 0.01)
HTdrug <- filter(HTdrug, Editing.rate <= .49 | Editing.rate >= 0.51)

WTvehicle <- filter(WTvehicle, Editing.rate <= .99 & Editing.rate >= 0.01)
WTvehicle<- filter(WTvehicle, Editing.rate <= .49 | Editing.rate >= 0.51)

WTdrug <- filter(WTdrug, Editing.rate <= .99 & Editing.rate >= 0.01)
WTdrug <- filter(WTdrug, Editing.rate <= .49 | Editing.rate >= 0.51)

#=============Calculate mean editing/sample/Coordinate=================================
HTvehicle$CHROM <-paste0("chr", HTvehicle$CHROM)
HTvehicle$Coordinate <- paste0(HTvehicle$CHROM, ":", HTvehicle$POS)
HTdrug$CHROM<-paste0("chr", HTdrug$CHROM)
HTdrug$Coordinate <- paste0(HTdrug$CHROM, ":", HTdrug$POS)
WTvehicle$CHROM <-paste0("chr", WTvehicle$CHROM)
WTvehicle$Coordinate <- paste0(WTvehicle$CHROM, ":", WTvehicle$POS)
WTdrug$CHROM<-paste0("chr", WTdrug$CHROM)
WTdrug$Coordinate <- paste0(WTdrug$CHROM, ":", WTdrug$POS)



HTvehicle_avg <- HTvehicle %>% group_by(Coordinate) %>%summarise(Mean.editing.rate = mean(Editing.rate))

HTvehicle<- left_join(HTvehicle_avg, 
                      HTvehicle %>% group_by(Coordinate) %>% 
                        summarise_at(vars(-group_cols()), .funs = ~paste(unique(.), collapse ="_")) %>% 
                        ungroup()
)

HTdrug_avg <- HTdrug %>% group_by(Coordinate) %>%summarise(Mean.editing.rate = mean(Editing.rate))

HTdrug<- left_join(HTdrug_avg, 
                   HTdrug %>% group_by(Coordinate) %>% 
                     summarise_at(vars(-group_cols()), .funs = ~paste(unique(.), collapse ="_")) %>% 
                     ungroup()
)

WTvehicle_avg <- WTvehicle %>% group_by(Coordinate) %>%summarise(Mean.editing.rate = mean(Editing.rate))

WTvehicle<- left_join(WTvehicle_avg, 
                      WTvehicle %>% group_by(Coordinate) %>% 
                        summarise_at(vars(-group_cols()), .funs = ~paste(unique(.), collapse ="_")) %>% 
                        ungroup()
)
WTdrug_avg <- WTdrug %>% group_by(Coordinate) %>%summarise(Mean.editing.rate = mean(Editing.rate))

WTdrug<- left_join(WTdrug_avg, 
                   WTdrug %>% group_by(Coordinate) %>% 
                     summarise_at(vars(-group_cols()), .funs = ~paste(unique(.), collapse ="_")) %>% 
                     ungroup()
)

################Filter to include sites with edited in >=50% samples
HTvehicle <- HTvehicle%>% mutate(No.of.Edited.Samples=str_count(Run, "SRR"))
HTdrug <- HTdrug%>% mutate(No.of.Edited.Samples=str_count(Run, "SRR"))
WTvehicle <- WTvehicle%>% mutate(No.of.Edited.Samples=str_count(Run, "SRR"))
WTdrug <- WTdrug%>% mutate(No.of.Edited.Samples=str_count(Run, "SRR"))

HTvehicle <- filter(HTvehicle, No.of.Edited.Samples > 1)
HTdrug <- filter(HTdrug, No.of.Edited.Samples > 1)
WTvehicle <- filter(WTvehicle, No.of.Edited.Samples > 1)
WTdrug <- filter(WTdrug, No.of.Edited.Samples > 1)

#=================savedata========================================================
setwd("~/Desktop/ARID1B_ASD/WholeSampleFiles")
FileList <- list(HTvehicle=HTvehicle, HTdrug=HTdrug, WTvehicle=WTvehicle,WTdrug=WTdrug)

for (i in names(FileList)) {write.csv(FileList[[i]], paste0(i,".csv"))
  
}


#=========================Import REDIportal varinat data============================
REDIportal <- read.delim("REDIportal_mm10.txt")
#================ filter for known genetic variations from dbsnp========================

REDIportal <- filter(REDIportal, dbsnp == "-")

#==================Add site annotations from REDIportal========================
REDIportal <- REDIportal[,c(1,2,5,9,10,11)]
REDIportal$Coordinate <- paste0(REDIportal$Region, ":", REDIportal$Position)
REDIportal <- REDIportal %>% relocate(Coordinate, 1)
##########to check duplicates in REDIportal File #################
length(unique(REDIportal$Coordinate)) == nrow(REDIportal)
print(REDIportal$Coordinate[duplicated(REDIportal$Coordinate)])
#=======================Annotate with REDIportal=======================
HTvehicle_annotated <- left_join(HTvehicle, REDIportal, by="Coordinate")
HTvehicle_annotated <- na.omit(HTvehicle_annotated)
write.csv(HTvehicle_annotated, "HTvehicle_annotated.csv")
HTdrug_annotated <- left_join(HTdrug, REDIportal, by="Coordinate")
HTdrug_annotated <- na.omit(HTdrug_annotated)
write.csv(HTdrug_annotated, "HTdrug_annotated.csv")

WTvehicle_annotated <- left_join(WTvehicle, REDIportal, by="Coordinate")
WTvehicle_annotated <- na.omit(WTvehicle_annotated)
write.csv(WTvehicle_annotated, "WTvehicle_annotated.csv")

WTdrug_annotated <- left_join(WTdrug, REDIportal, by="Coordinate")
WTdrug_annotated <- na.omit(WTdrug_annotated)
write.csv(WTdrug_annotated, "WTdrug_annotated.csv")
##################Combine  samples for unique sites#################
All_sample_data_but_HTdrug <- rbind(HTvehicle_annotated, WTdrug_annotated, WTvehicle_annotated)
All_sample_data_but_HTveh <- rbind(HTdrug_annotated, WTdrug_annotated, WTvehicle_annotated)
All_sample_data_but_WTdrug <- rbind(HTdrug_annotated,HTvehicle_annotated,WTvehicle_annotated)
All_sample_data_but_WTveh <- rbind(HTdrug_annotated,HTvehicle_annotated ,WTdrug_annotated)
#========== Shared sites between HT and WT=================================

###Rename editing rate column##########
names(HTdrug_annotated)[2] <-"Mean.editing.rate.HTD"
names(HTvehicle_annotated)[2] <-"Mean.editing.rate.HTV"
names(WTdrug_annotated)[2] <- "Mean.editing.rate.WTD"
names(WTvehicle_annotated) [2] <-"Mean.editing.rate.WTV"
###Shared sites between HT/drug and WT/drug group
Drug_HT_WT_shared <- left_join(HTdrug_annotated, WTdrug_annotated, by="Coordinate")
Drug_HT_WT_shared <-Drug_HT_WT_shared %>% drop_na(Mean.editing.rate.WTD)

##make a function to calculate standard deviation
Var <- function(x){
  values<- as.numeric(unlist(strsplit(x, "_")))
  var(value)
}

###Drug shared

Drug_HT_WT_shared$var.HTD <- sapply(Drug_HT_WT_shared$Editing.rate.x, Var)
Drug_HT_WT_shared$var.WTD <- sapply(Drug_HT_WT_shared$Editing.rate.y, Var)

Drug_HT_WT_shared <- Drug_HT_WT_shared%>% mutate(MeanDiff= Mean.editing.rate.HTD- Mean.editing.rate.WTD,
                                                 SE= sqrt(var.HTD) / length(Mean.editing.rate.HTD)+
                                                   (var.WTD)/ (length(Mean.editing.rate.WTD)),t_critical= qt(0.975, df = 2 - 1),
                                                 Upper_95_CI = MeanDiff + t_critical * SE,
                                                 Lower_95_CI = MeanDiff - t_critical * SE
)


###Filter CI
##if upper x lower < 0, CI limit has zero
Drug_HT_WT_shared <- Drug_HT_WT_shared %>% mutate(upperxlower= Upper_95_CI*Lower_95_CI, Has_zero= upperxlower < 0)
Drug_HT_WT_shared <- filter(Drug_HT_WT_shared, Has_zero=="FALSE")
setwd("~/Desktop/ARID1B_ASD/Differential_Editing")
write.csv(Drug_HT_WT_shared, "Drug_HT_WT_sharedSites.csv")

###Shared sites between HT/veh and WT/veh group
Vehicle_HT_WT_shared <- left_join(HTvehicle_annotated, WTvehicle_annotated, by="Coordinate")
Vehicle_HT_WT_shared <- Vehicle_HT_WT_shared %>% drop_na(Mean.editing.rate.WTV)


Vehicle_HT_WT_shared$var.HTV<- sapply(Vehicle_HT_WT_shared$Editing.rate.x, Var)
Vehicle_HT_WT_shared$var.WTV <- sapply(Vehicle_HT_WT_shared$Editing.rate.y, Var)
Vehicle_HT_WT_shared <- Vehicle_HT_WT_shared%>% mutate(MeanDiff= Mean.editing.rate.HTV - Mean.editing.rate.WTV,
                                                       SE= sqrt(var.HTV) / length(Mean.editing.rate.HTV)+
                                                         (var.WTV)/ (length(Mean.editing.rate.WTV)),t_critical= qt(0.975, df = 2 - 1),
                                                       Upper_95_CI = MeanDiff + t_critical * SE,
                                                       Lower_95_CI = MeanDiff - t_critical * SE
)

###Filter CI

##if upper x lower < 0, CI limit has zero

Vehicle_HT_WT_shared <- Vehicle_HT_WT_shared %>% mutate(upperxlower= Upper_95_CI*Lower_95_CI, Has_zero= upperxlower < 0)
Vehicle_HT_WT_shared <- filter(Vehicle_HT_WT_shared, Has_zero=="FALSE")
write.csv(Vehicle_HT_WT_shared, "Vehicle_HT_WT_sharedSites.csv")

############Unique sites####################
#Find unique sites for HT drug group
HT_drug_unique <- anti_join(HTdrug_annotated,All_sample_data_but_HTdrug, by="Coordinate")
HT_drug_unique <- na.omit(HT_drug_unique)


#Find Unique sites for HT vehicle group
HT_vehicle_unique <- anti_join(HTvehicle_annotated, All_sample_data_but_HTveh, by="Coordinate")
HT_vehicle_unique <- na.omit(HT_vehicle_unique)
#Save results
write.csv(HT_drug_unique, "HT_drug_uniqueSites.csv")
write.csv(HT_vehicle_unique, "HT_vehicle_uniqueSites.csv")

#Find WT unique sites for drug group
WT_drug_unique <- anti_join(WTdrug_annotated,All_sample_data_but_WTdrug, by="Coordinate")
WT_drug_unique <- na.omit(WT_drug_unique)


#Find Unique sites for WT vehicle group
WT_vehicle_unique <- anti_join(WTvehicle_annotated, All_sample_data_but_WTveh, by="Coordinate")
WT_vehicle_unique <- na.omit(WT_vehicle_unique)
#Save results
write.csv(WT_drug_unique, "WT_drug_uniqueSites.csv")
write.csv(WT_vehicle_unique, "WT_vehicle_uniqueSites.csv")

#############################################################mRNA-miRNA Interaction Analaysis#################

###Extract FASTA seq. for miRNA binding analysis

###Load ref geome/create function to extract fasta
#Load ref_genome
ref_genome <- readDNAStringSet("mm10.fa", format = "fasta")
ref_genome

#Define window size (+- 100bp from edited position)
window_size <- 100
extract_FASTA <- function(chr, position,ref_genome) {  
  start_pos <- position - window_size
  end_pos <- position + window_size
  start_pos <- max(start_pos, 1)
  end_pos <- min(end_pos, length(ref_genome[[chr]])) 
  seq <- ref_genome[[chr]][start_pos:end_pos] 
  return(as.character(seq))
}

#### Drug-treated-shared (only sites filtered for CI)
##Make a new data frame with Chr, Gene, and POS column
Drug_HT_WT_shared_FilteredCI <- Drug_HT_WT_sharedSites
Vehicle_HT_WT_shared_FilteredCI <- Vehicle_HT_WT_sharedSites

Fasta_Drug_HT_WT_shared <- data.frame(Drug_HT_WT_shared_FilteredCI$CHROM.x, Drug_HT_WT_shared_FilteredCI$POS.x, Drug_HT_WT_shared_FilteredCI$Gene.wgEncodeGencodeBasicVM16.x)
########change colnames for clarity################
names(Fasta_Drug_HT_WT_shared)[1]<- "chromosome"
names(Fasta_Drug_HT_WT_shared)[2]<- "position"
names(Fasta_Drug_HT_WT_shared)[3]<- "Gene"
##Extract FASTA
Fasta_Drug_HT_WT_shared$sequence.fasta <- mapply(extract_FASTA,Fasta_Drug_HT_WT_shared$chromosome , Fasta_Drug_HT_WT_shared$position, MoreArgs = list(ref_genome = ref_genome))

#### Vehicle-treated-group
Fasta_Vehicle_HT_WT_shared <- data.frame(Vehicle_HT_WT_shared_FilteredCI$CHROM.x, Vehicle_HT_WT_shared_FilteredCI$POS.x, Vehicle_HT_WT_shared_FilteredCI$Gene.wgEncodeGencodeBasicVM16.x)
names(Fasta_Vehicle_HT_WT_shared)[1]<- "chromosome"
names(Fasta_Vehicle_HT_WT_shared)[2]<- "position"
names(Fasta_Vehicle_HT_WT_shared)[3]<- "Gene"

##Extract FASTA
Fasta_Vehicle_HT_WT_shared$sequence.fasta <- mapply(extract_FASTA,Fasta_Vehicle_HT_WT_shared$chromosome , Fasta_Vehicle_HT_WT_shared$position, MoreArgs = list(ref_genome = ref_genome))
###Change A>G or T>C at edited position for edited version of  FASTA seq
#Edited FASTA (Drug_treated_group)
Fasta_Drug_HT_WT_shared$edited_sequence <- ifelse(
  nchar(Fasta_Drug_HT_WT_shared$sequence.fasta) >= 101,
  paste0(
    substr(Fasta_Drug_HT_WT_shared$sequence.fasta, 1, 100),
    ifelse(
      substr(Fasta_Drug_HT_WT_shared$sequence.fasta, 101, 101) == "A", "G",
      ifelse(
        substr(Fasta_Drug_HT_WT_shared$sequence.fasta, 101, 101) == "T", "C",
        substr(Fasta_Drug_HT_WT_shared$sequence.fasta, 101, 101)
      )
    ),
    substr(Fasta_Drug_HT_WT_shared$sequence.fasta, 102, nchar(Fasta_Drug_HT_WT_shared$sequence.fasta))
  ),
  Fasta_Drug_HT_WT_shared$sequence.fasta
)

#Edited FASTA (Vehicle_treated_group)
Fasta_Vehicle_HT_WT_shared$edited_sequence <- ifelse(
  nchar(Fasta_Vehicle_HT_WT_shared$sequence.fasta) >= 101,
  paste0(
    substr(Fasta_Vehicle_HT_WT_shared$sequence.fasta, 1, 100),
    ifelse(
      substr(Fasta_Vehicle_HT_WT_shared$sequence.fasta, 101, 101) == "A", "G",
      ifelse(
        substr(Fasta_Vehicle_HT_WT_shared$sequence.fasta, 101, 101) == "T", "C",
        substr(Fasta_Vehicle_HT_WT_shared$sequence.fasta, 101, 101)
      )
    ),
    substr(Fasta_Vehicle_HT_WT_shared$sequence.fasta, 102, nchar(Fasta_Vehicle_HT_WT_shared$sequence.fasta))
  ),
  Fasta_Vehicle_HT_WT_shared$sequence.fasta
)


########From Table to Fasta################
##Drug_group
Fasta_edited_Drug_HT_WT_shared <- data.frame(Fasta_Drug_HT_WT_shared$Gene, Fasta_Drug_HT_WT_shared$edited_sequence)
names(Fasta_edited_Drug_HT_WT_shared)[1]<- "Gene"
names(Fasta_edited_Drug_HT_WT_shared)[2]<- "edited.fasta"
Fasta_edited_Drug_HT_WT_shared = dataframe2fas(Fasta_edited_Drug_HT_WT_shared, file="Drug_shared_edited_seqlist.fasta")

Fasta_unedited_Drug_HT_WT_shared <- data.frame(Fasta_Drug_HT_WT_shared$Gene, Fasta_Drug_HT_WT_shared$sequence.fasta)
names(Fasta_unedited_Drug_HT_WT_shared)[1]<- "Gene"
names(Fasta_unedited_Drug_HT_WT_shared)[2]<- "unedited.fasta"
Fasta_unedited_Drug_HT_WT_shared = dataframe2fas(Fasta_unedited_Drug_HT_WT_shared, file="Fasta_unedited_Drug_HT_WT_shared")

write.fasta(Fasta_unedited_Drug_HT_WT_shared,"Drug_shared_unedited_seqlist.fasta")
##vehicle_group
Fasta_edited_vehicle_HT_WT_shared <- data.frame(Fasta_Vehicle_HT_WT_shared$Gene, Fasta_Vehicle_HT_WT_shared$edited_sequence)
names(Fasta_edited_vehicle_HT_WT_shared)[1]<- "Gene"
names(Fasta_edited_vehicle_HT_WT_shared)[2]<- "edited.fasta"
Fasta_edited_vehicle_HT_WT_shared = dataframe2fas(Fasta_edited_vehicle_HT_WT_shared, file="MyFasta")

write.fasta(Fasta_edited_vehicle_HT_WT_shared,"vehicle_shared_edited_seqlist.fasta")


Fasta_unedited_vehicle_HT_WT_shared <- data.frame(Fasta_Vehicle_HT_WT_shared$Gene,Fasta_Vehicle_HT_WT_shared$sequence.fasta)
names(Fasta_unedited_vehicle_HT_WT_shared)[1]<- "Gene"
names(Fasta_unedited_vehicle_HT_WT_shared)[2]<- "unedited.fasta"
Fasta_unedited_vehicle_HT_WT_shared = dataframe2fas(Fasta_unedited_vehicle_HT_WT_shared, file="Fasta_unedited_vehicle_HT_WT_shared")

write.fasta(Fasta_unedited_vehicle_HT_WT_shared,"vehicle_shared_unedited_seqlist.fasta")

####################################
#####Unique sites FASTA#######
##########GET FASTA##################
#######
names(HT_drug_unique)[20]<- "Gene"
names(HT_vehicle_unique)[20]<- "Gene"


FASTA_HTDU <- HT_drug_unique[, c("POS","CHROM","Gene")]
names(FASTA_HTDU) <- c("Position","chromosome","Gene")

FASTA_HTVU <- HT_vehicle_unique[, c("POS","CHROM","Gene")]
names(FASTA_HTVU) <- c("Position","chromosome","Gene")

##Add a new column for chromosome with prefix "chr"
FASTA_HTDU$sequence.fasta <- mapply(extract_FASTA,FASTA_HTDU$chromosome , FASTA_HTDU$Position, MoreArgs = list(ref_genome = ref_genome))
FASTA_HTVU$sequence.fasta <- mapply(extract_FASTA, FASTA_HTVU$chromosome,FASTA_HTVU$Position, MoreArgs = list(ref_genome=ref_genome))

#Edited FASTA
####HTDU
FASTA_HTDU$edited_sequence <- ifelse(nchar(FASTA_HTDU$sequence.fasta) >= 101,
                                               paste0(
                                                 substr(FASTA_HTDU$sequence.fasta, 1, 100),
                                                 ifelse(substr(FASTA_HTDU$sequence.fasta, 101, 101) == "A", "G",
                                                        ifelse(substr(FASTA_HTDU$sequence.fasta, 101, 101) == "T", "C",
                                                               substr(FASTA_HTDU$sequence.fasta, 101, 101))),
                                                 substr(FASTA_HTDU$sequence.fasta, 102, nchar(FASTA_HTDU$sequence.fasta))
                                               ),
                                               FASTA_HTDU$sequence.fasta
)

#####HTVU
FASTA_HTVU$edited_sequence <- ifelse(nchar(FASTA_HTVU$sequence.fasta) >= 101,
                                                  paste0(
                                                    substr(FASTA_HTVU$sequence.fasta, 1, 100),
                                                    ifelse(substr(FASTA_HTVU$sequence.fasta, 101, 101) == "A", "G",
                                                           ifelse(substr(FASTA_HTVU$sequence.fasta, 101, 101) == "T", "C",
                                                                  substr(FASTA_HTVU$sequence.fasta, 101, 101))),
                                                    substr(FASTA_HTVU$sequence.fasta, 102, nchar(FASTA_HTVU$sequence.fasta))
                                                  ),
                                                  FASTA_HTVU$sequence.fasta
)

#####From table to fasta
Fasta_HTDU_edited <- data.frame(FASTA_HTDU$Gene, FASTA_HTDU$edited_sequence)
names(Fasta_HTDU_edited) <- c("Gene","edited.fasta" )
Fasta_HTDU_edited= dataframe2fas(Fasta_HTDU_edited, file="Fasta_HTDU_edited.fasta")
setwd("~/Desktop/ASD_miRNA/remoteRNA22v2")
write.fasta(Fasta_HTDU_edited, "HTDU_edited.txt")

Fasta_HTDU_Unedited <- data.frame(FASTA_HTDU$Gene, FASTA_HTDU$sequence.fasta)
names(Fasta_HTDU_Unedited) <- c("Gene","Unedited.fasta" )
Fasta_HTDU_Unedited= dataframe2fas(Fasta_HTDU_Unedited, file="Fasta_HTDU_Unedited.fasta")
write.fasta(Fasta_HTDU_Unedited, "HTDU_unedited.txt")

########HTVU
Fasta_HTVU_edited <- data.frame(FASTA_HTVU$Gene, FASTA_HTVU$edited_sequence)
names(Fasta_HTVU_edited) <- c("Gene", "edited.fasta")
Fasta_HTVU_edited= dataframe2fas(Fasta_HTVU_edited, file="Fasta_HTVU_edited.fasta")
write.fasta(Fasta_HTVU_edited, "HTVU_edited.txt")
Fasta_HTVU_unedited <- data.frame(FASTA_HTVU$Gene, FASTA_HTVU$sequence.fasta)
names(Fasta_HTVU_unedited) <- c("Gene", "unedited.fasta")
Fasta_HTVU_unedited= dataframe2fas(Fasta_HTVU_unedited, file="Fasta_HTVU_unedited.txt")
write.fasta(Fasta_HTVU_unedited, "HTVU_unedited.txt")


######get fasta for wt unique
names(WT_drug_uniqueSites)[19]<- "Gene"
names(WT_vehicle_uniqueSites)[19]<- "Gene"

FASTA_WTDU <- WT_drug_uniqueSites [, c("POS","CHROM","Gene" )]
names(FASTA_WTDU) <- c("Position","chromosome","Gene")

FASTA_WTVU <- WT_vehicle_uniqueSites[, c("POS","CHROM","Gene")]
names(FASTA_WTVU) <- c("Position","chromosome","Gene")

##Extract FASTA
FASTA_WTDU$sequence.fasta <- mapply(extract_FASTA,FASTA_WTDU$chromosome , FASTA_WTDU$Position, MoreArgs = list(ref_genome = ref_genome))
FASTA_WTVU$sequence.fasta <- mapply(extract_FASTA,FASTA_WTVU$chromosome , FASTA_WTVU$Position, MoreArgs = list(ref_genome = ref_genome))



###Change A>G or T>C at edited position for edited version of  FASTA seq
#Edited FASTA (WTDU)
FASTA_WTDU$edited_sequence <- ifelse(
  nchar(FASTA_WTDU$sequence.fasta) >= 101,
  paste0(
    substr(FASTA_WTDU$sequence.fasta, 1, 100),
    ifelse(
      substr(FASTA_WTDU$sequence.fasta, 101, 101) == "A", "G",
      ifelse(
        substr(FASTA_WTDU$sequence.fasta, 101, 101) == "T", "C",
        substr(FASTA_WTDU$sequence.fasta, 101, 101)
      )
    ),
    substr(FASTA_WTDU$sequence.fasta, 102, nchar(FASTA_WTDU$sequence.fasta))
  ),
  FASTA_WTDU$sequence.fasta
)

#####WTVU

FASTA_WTVU$edited_sequence <- ifelse(
  nchar(FASTA_WTVU$sequence.fasta) >= 101,
  paste0(
    substr(FASTA_WTVU$sequence.fasta, 1, 100),
    ifelse(
      substr(FASTA_WTVU$sequence.fasta, 101, 101) == "A", "G",
      ifelse(
        substr(FASTA_WTVU$sequence.fasta, 101, 101) == "T", "C",
        substr(FASTA_WTVU$sequence.fasta, 101, 101)
      )
    ),
    substr(FASTA_WTVU$sequence.fasta, 102, nchar(FASTA_WTVU$sequence.fasta))
  ),
  FASTA_WTVU$sequence.fasta
)


#####From table to fasta
####WTD
Fasta_WTDU_edited <- data.frame(FASTA_WTDU$Gene, FASTA_WTDU$edited_sequence)
names(Fasta_WTDU_edited) <- c("Gene","edited.fasta" )
Fasta_WTDU_edited= dataframe2fas(Fasta_WTDU_edited, file="Fasta_WTDU_edited.fasta")
write.fasta(Fasta_WTDU_edited, "WTDU_edited.txt")

Fasta_WTDU_Unedited <- data.frame(FASTA_WTDU$Gene, FASTA_WTDU$sequence.fasta)
names(Fasta_WTDU_Unedited) <- c("Gene","Unedited.fasta" )
Fasta_WTDU_Unedited= dataframe2fas(Fasta_WTDU_Unedited, file="Fasta_WTDU_Unedited.fasta")
write.fasta(Fasta_WTDU_Unedited, "WTDU_unedited.txt")
########WTVU
Fasta_WTVU_edited <- data.frame(FASTA_WTVU$Gene, FASTA_WTVU$edited_sequence)
names(Fasta_WTVU_edited) <- c("Gene", "edited.fasta")
Fasta_WTVU_edited= dataframe2fas(Fasta_WTVU_edited, file="Fasta_WTVU_edited.fasta")
write.fasta(Fasta_WTVU_edited, "WTVU_edited.txt")

Fasta_WTVU_unedited <- data.frame(FASTA_WTVU$Gene, FASTA_WTVU$sequence.fasta)
names(Fasta_WTVU_unedited) <- c("Gene", "unedited.fasta")
Fasta_WTVU_unedited= dataframe2fas(Fasta_WTVU_unedited, file="Fasta_WTVU_unedited.fasta")
write.fasta(Fasta_WTVU_unedited, "WTVU_unedited.txt")

#########################END#################################
