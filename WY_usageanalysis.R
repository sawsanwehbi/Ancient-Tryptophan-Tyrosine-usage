# Load the ape package
library(ape)
library(protr)
library(seqinr)
library(tidyr)
library(tidyverse)
library(matrixStats)
library(phytools)
library(phangorn)
library(common)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(stringr)
library(gridtext)
library(grid)
library(ggpubr)
library(bio3d)
library(dplyr)
library(weights)
library(philentropy)
library(deming)
library(scales)
library(diagis)
library("msa")
library('geoR')

Bacterial_supergroups <- c('CPR', 'PVC', 'FCB', 'Terrabacteria', 'Proteobacteria')
Archaeal_supergroups <- c('Asgard', 'TACK', 'DPANN', 'Euryarchaeota')
All_supergroups <- c(Bacterial_supergroups,Archaeal_supergroups)
AA_properties <- read.csv('../Pfam Trees/AminoAcid_properties.csv', header = T)

AncientPostLUCA <- read.csv('../Tryptophan paper/AncientPostLUCA.csv', header = T)
PfamConAACNQbac <- read.csv('pfam_asr_aac_NQBac_0.4_0.35preLBCA',header=T)
PfamConAACNQarc <- read.csv('pfam_asr_aac_NQArc_0.4_0.35preLBCA',header=T)

Pfam_ConAAC <- PfamConAACNQarc
Pfam_ConAAC$median_PfamLen <- Pfam_contemp$median_PfamLen

### Calculate clan properties ####
uniqueClans <- unique(Pfam_ConAAC$clans)
Clan_ConAAfreq_Df <- data.frame(uniqueClans)
for (aa in 2:21) {
  Clan_ConAAfreq <- sapply(1:length(uniqueClans), function (i){
    sum(Pfam_ConAAC[,aa][which(Pfam_ConAAC$clans == uniqueClans[i])]*
          Pfam_ConAAC$Conserved_length[which(Pfam_ConAAC$clans == uniqueClans[i])])/
      sum(Pfam_ConAAC$Conserved_length[which(Pfam_ConAAC$clans == uniqueClans[i])])})
  Clan_ConAAfreq_Df <- cbind(Clan_ConAAfreq_Df , Clan_ConAAfreq) }
colnames(Clan_ConAAfreq_Df)[2:21] <- colnames(Pfam_ConAAC)[2:21]
Clan_ConAAfreq_Df$Clan_Conlength <- sapply(1:length(uniqueClans), function (i){
  max(Pfam_ConAAC$Conserved_length[which(Pfam_ConAAC$clans == uniqueClans[i])])})
colnames(Clan_ConAAfreq_Df)[1] <- 'Clans'
Clan_ConAAfreq_Df$gap <- sapply(1:length(uniqueClans), function (i){
  sum(Pfam_ConAAC$gap[which(Pfam_ConAAC$clans == uniqueClans[i])])})
Clan_ConAAfreq_Df$contemp_len <- sapply(1:length(uniqueClans), function (i){
  sum(Pfam_ConAAC$median_PfamLen[which(Pfam_ConAAC$clans == uniqueClans[i])])})
Clan_ConAAfreq_Df$concat_len <- sapply(1:length(uniqueClans), function (i){
  sum(Pfam_ConAAC$Conserved_length[which(Pfam_ConAAC$clans == uniqueClans[i])])})
Clan_ancestor <- vector() 
for ( i in 1:length(uniqueClans)) {
  if (length(which(grepl('LUCA', Pfam_ConAAC$ancestor[which(Pfam_ConAAC$clans %in% uniqueClans[i])]))) > 1)
  {Clan_ancestor[i] <- 'preLUCA'} 
  else if (length(which(grepl('preLUCA', Pfam_ConAAC$ancestor[which(Pfam_ConAAC$clans %in% uniqueClans[i])]))) > 0)
  { Clan_ancestor[i] <- 'preLUCA'} 
  else if (length(which(grepl('LUCA', Pfam_ConAAC$ancestor[which(Pfam_ConAAC$clans %in% uniqueClans[i])]))) == 1)
  {Clan_ancestor[i] <- 'LUCA'} 
  else if (  'LBCA' %in% Pfam_ConAAC$ancestor[which(Pfam_ConAAC$clans %in% uniqueClans[i])] &
             'LACA' %in% Pfam_ConAAC$ancestor[which(Pfam_ConAAC$clans %in% uniqueClans[i])])
  {Clan_ancestor[i] <- 'LUCA'} 
  else if (  'LBCA' %in% Pfam_ConAAC$ancestor[which(Pfam_ConAAC$clans %in% uniqueClans[i])] |
             'LACA' %in% Pfam_ConAAC$ancestor[which(Pfam_ConAAC$clans %in% uniqueClans[i])]
             & !'unclassifiable' %in% Pfam_ConAAC$ancestor[which(Pfam_ConAAC$Clans %in% uniqueClans[i])]) 
  {Clan_ancestor[i] <- 'postLUCA'} 
  else if((length(which(Pfam_ConAAC$ancestor[which(Pfam_ConAAC$clans %in% uniqueClans[i])] %in% c('post-LBCA',Bacterial_supergroups))) > 2) &
          (length(which(Pfam_ConAAC$ancestor[which(Pfam_ConAAC$clans %in% uniqueClans[i])] %in% Archaeal_supergroups)) > 1))
  {Clan_ancestor[i] <- 'unclassifiable' }
  else if((length(which(unique(Pfam_ConAAC$ancestor[which(Pfam_ConAAC$clans %in% uniqueClans[i])]) %in% c('post-LBCA',Bacterial_supergroups))) > 2) |
          (length(which(unique(Pfam_ConAAC$ancestor[which(Pfam_ConAAC$clans %in% uniqueClans[i])]) %in% Archaeal_supergroups)) > 1) 
          & !'unclassifiable' %in% Pfam_ConAAC$ancestor[which(Pfam_ConAAC$clans %in% uniqueClans[i])]) 
  {Clan_ancestor[i] <- 'postLUCA' }
  else if((length(which(unique(Pfam_ConAAC$ancestor[which(Pfam_ConAAC$clans %in% uniqueClans[i])]) %in% c('post-LBCA',Bacterial_supergroups))) == 1) |
          (length(which(unique(Pfam_ConAAC$ancestor[which(Pfam_ConAAC$clans %in% uniqueClans[i])]) %in% Archaeal_supergroups)) == 1) | 
          (length(which(unique(Pfam_ConAAC$ancestor[which(Pfam_ConAAC$clans %in% uniqueClans[i])]) %in% c('post-LBCA',Bacterial_supergroups))) == 2)
          & !'unclassifiable' %in% Pfam_ConAAC$ancestor[which(Pfam_ConAAC$clans %in% uniqueClans[i])]) 
  {Clan_ancestor[i] <- 'modern' }
  else {  Clan_ancestor[i] <- 'unclassifiable'} }
Clan_ConAAfreq_Df$Clan_ancestor <- Clan_ancestor
#write.csv(Clan_ConAAfreq_Df, 'Clan_data_ancestralAAC_NQbac.csv', row.names = F)

#### Pre-LACA and Pre-LBCA clan classification ######
postLUCAclans <- Clan_ConAAfreq_Df$Clans[which(Clan_ConAAfreq_Df$Clan_ancestor == 'postLUCA')]
LBCAclans <- vector()
LACAclans <- vector()
for ( i in 1:length(postLUCAclans)){
  if (any(Pfam_ConAAC$ancestor[which(Pfam_ConAAC$clans == postLUCAclans[i])] =='LBCA'))
  {LBCAclans[i] <- postLUCAclans[i]} 
  else if (any(Pfam_ConAAC$ancestor[which(Pfam_ConAAC$clans == postLUCAclans[i])] =='LACA'))
  {LACAclans[i] <- postLUCAclans[i]} else {next} }
LBCAclans <- LBCAclans[-which(is.na(LBCAclans))]
LACAclans <- LACAclans[-which(is.na(LACAclans))]
BroaddiversifiedLACAclans <- LACAclans[which(grepl('CL....', LACAclans))] #20
BroaddiversifiedLBCAclans <- LBCAclans[which(grepl('CL....', LBCAclans))] #156
#Broaddiversifiedclans <- c(BroaddiversifiedLBCAclans, BroaddiversifiedLACAclans) 
## broader classification of ancient postLUCA

## PostLUCA clans with >1 LACA/LBCA or 1 preLACA/preLBCA pfam; stricter classification
MultipleLACAclans <- sapply(1:length(BroaddiversifiedLACAclans), function(i){
  length(which(Pfam_ConAAC$ancestor[which(Pfam_ConAAC$clans %in% BroaddiversifiedLACAclans[i])] == 'LACA'))
})
MultipleLBCAclans <- sapply(1:length(BroaddiversifiedLBCAclans), function(i){
  length(which(Pfam_ConAAC$ancestor[which(Pfam_ConAAC$clans %in% BroaddiversifiedLBCAclans[i])] == 'LBCA'))
})
diversifiedLACAclans <- BroaddiversifiedLACAclans[which(MultipleLACAclans  > 1)] #4
diversifiedLACAclans <-  unique(c(diversifiedLACAclans ,unique(Pfam_ConAAC$clans[which(Pfam_ConAAC$pfamIDs %in% AncientPostLUCA$PFAM_IDs[which(
  AncientPostLUCA$classified_ancestor == 'preLACA')])]) )) #53
diversifiedLACAclans <- Clan_ConAAfreq_Df$Clans[which( Clan_ConAAfreq_Df$Clans %in% diversifiedLACAclans & 
                                                         Clan_ConAAfreq_Df$Clan_ancestor == 'postLUCA')] #31

diversifiedLBCAclans <- BroaddiversifiedLBCAclans[which(MultipleLBCAclans  > 1)] #77
diversifiedLBCAclans <-  unique(c(diversifiedLBCAclans , unique(Pfam_ConAAC$clans[which(Pfam_ConAAC$pfamIDs %in% AncientPostLUCA$PFAM_IDs[which(
  AncientPostLUCA$classified_ancestor == 'preLBCA')])])))#364
diversifiedLBCAclans <- Clan_ConAAfreq_Df$Clans[which( Clan_ConAAfreq_Df$Clans %in% diversifiedLBCAclans & 
                                                         Clan_ConAAfreq_Df$Clan_ancestor == 'postLUCA')] #251
diversifiedclans <- c(diversifiedLBCAclans, diversifiedLACAclans) 
# diversifiedLBCAclans and diversifiedLACAclans are preLBCA and preLACA clans

BroaddiversifiedLACAclans <- unique(c(BroaddiversifiedLACAclans, diversifiedLACAclans )) #39
BroaddiversifiedLBCAclans <- unique(c(BroaddiversifiedLBCAclans, diversifiedLBCAclans ))  #315
Broaddiversifiedclans <- c(BroaddiversifiedLBCAclans, BroaddiversifiedLACAclans) 
# Broaddiversifiedclans include diverged post-LUCA clans that do not exactly fit or pre-LACA/pre-LBCA criteria
# but we are unsure enough about their age that we can exclude them when using postLUCA AAC for usage ratios

# Modern archaea and bacteria clans
modernClans <- Clan_ConAAfreq_Df$Clans[which(Clan_ConAAfreq_Df$Clan_ancestor == 'modern')]
#modernClans <- modernClans[which(grepl('PF.....', modernClans))] ## SINGLE COPY MODERN
modernBacteriaclans <- Pfam_ConAAC$clans[which(Pfam_ConAAC$clans %in% modernClans 
                                               & Pfam_ConAAC$ancestor %in% Bacterial_supergroups | 
                                                 Pfam_ConAAC$clans %in% modernClans & Pfam_ConAAC$ancestor == 'post-LBCA')]
modernArchaeaclans <- Pfam_ConAAC$clans[which(Pfam_ConAAC$clans %in% modernClans 
                                              & Pfam_ConAAC$ancestor %in% Archaeal_supergroups)]
CommonArcBacclans <- which( modernArchaeaclans %in% modernBacteriaclans)
modernArchaeaclans  <- modernArchaeaclans[-CommonArcBacclans] #163
modernBacteriaclans <- modernBacteriaclans[-CommonArcBacclans] #2114
modernBacteriaclans <- modernBacteriaclans[-which(modernBacteriaclans %in%
                                                    unique(Pfam_ConAAC$clans[which(Pfam_ConAAC$ancestor == 'post-LBCA')]) )]#1578


### AAC per age group ####
Clan_ConAAfreq_Df <- read.csv('Clan_data_ancestralAAC_NQarc.csv')
#Clan_ConAAfreq_Df <- read.csv('Clan_data_ancestralAAC_NQbac.csv')
NoDivpostLUCA_Clans_ConAAC <- colWeightedMeans(as.matrix(Clan_ConAAfreq_Df[which(Clan_ConAAfreq_Df$Clan_ancestor == 'postLUCA' &
                                                                                   !Clan_ConAAfreq_Df$Clans %in% Broaddiversifiedclans ),2:21]),
                                               Clan_ConAAfreq_Df$Clan_Conlength[which(Clan_ConAAfreq_Df$Clan_ancestor == 'postLUCA' &
                                                                                        !Clan_ConAAfreq_Df$Clans %in% Broaddiversifiedclans)] )
allLUCA_Clans_ConAAC <- colWeightedMeans(as.matrix(Clan_ConAAfreq_Df[which(grepl('LUCA',Clan_ConAAfreq_Df$Clan_ancestor) & 
                         Clan_ConAAfreq_Df$Clan_ancestor != 'postLUCA'),2:21]),
                        Clan_ConAAfreq_Df$Clan_Conlength[which(grepl('LUCA',Clan_ConAAfreq_Df$Clan_ancestor) &
                          Clan_ConAAfreq_Df$Clan_ancestor != 'postLUCA')] )
LUCA_Clans_ConAAC <- colWeightedMeans(as.matrix(Clan_ConAAfreq_Df[which(Clan_ConAAfreq_Df$Clan_ancestor == 'LUCA'),2:21]),
                                      Clan_ConAAfreq_Df$Clan_Conlength[which(Clan_ConAAfreq_Df$Clan_ancestor == 'LUCA')] )
preLUCA_Clans_ConAAC <- colWeightedMeans(as.matrix(Clan_ConAAfreq_Df[which(Clan_ConAAfreq_Df$Clan_ancestor == 'preLUCA'),2:21]),
                                         Clan_ConAAfreq_Df$Clan_Conlength[which(Clan_ConAAfreq_Df$Clan_ancestor == 'preLUCA')] )
postLUCA_Clans_ConAAC <- colWeightedMeans(as.matrix(Clan_ConAAfreq_Df[which(Clan_ConAAfreq_Df$Clan_ancestor == 'postLUCA'),2:21]),
                                          Clan_ConAAfreq_Df$Clan_Conlength[which(Clan_ConAAfreq_Df$Clan_ancestor == 'postLUCA')] )
modern_Clans_ConAAC <- colWeightedMeans(as.matrix(Clan_ConAAfreq_Df[which(Clan_ConAAfreq_Df$Clan_ancestor == 'modern'),2:21]),
                                        Clan_ConAAfreq_Df$Clan_Conlength[which(Clan_ConAAfreq_Df$Clan_ancestor == 'modern')] )
modernArc_Clans_ConAAC <- colWeightedMeans(as.matrix(Clan_ConAAfreq_Df[which(Clan_ConAAfreq_Df$Clans %in% modernArchaeaclans ),2:21]),
                                           Clan_ConAAfreq_Df$Clan_Conlength[which( Clan_ConAAfreq_Df$Clans %in% modernArchaeaclans)] )
modernBac_Clans_ConAAC <- colWeightedMeans(as.matrix(Clan_ConAAfreq_Df[which(Clan_ConAAfreq_Df$Clans %in% modernBacteriaclans ),2:21]),
                                           Clan_ConAAfreq_Df$Clan_Conlength[which( Clan_ConAAfreq_Df$Clans %in% modernBacteriaclans)] )
DivLACA_ClanAAC <- colWeightedMeans(as.matrix(Clan_ConAAfreq_Df[which(Clan_ConAAfreq_Df$Clans %in% diversifiedLACAclans ),2:21]),
                                    Clan_ConAAfreq_Df$Clan_Conlength[which(Clan_ConAAfreq_Df$Clans %in% diversifiedLACAclans )] )
DivLBCA_ClanAAC <- colWeightedMeans(as.matrix(Clan_ConAAfreq_Df[which(Clan_ConAAfreq_Df$Clans %in% diversifiedLBCAclans ),2:21]),
                                    Clan_ConAAfreq_Df$Clan_Conlength[which(Clan_ConAAfreq_Df$Clans %in% diversifiedLBCAclans )] )
OrderedmodernBac_Clans_ConAAC <- modernBac_Clans_ConAAC[match(  AA_properties$Letter, names(modernBac_Clans_ConAAC))]
OderedmodernArc_Clans_ConAAC <- modernArc_Clans_ConAAC[match(  AA_properties$Letter, names(modernArc_Clans_ConAAC))]
## recent bacteria include diversified clans but with only 1 LBCA/LACA clan ie not old enough for pre-LBCA/LACA
RecentBacteria_Clans_ConAAC <- colWeightedMeans(as.matrix(Clan_ConAAfreq_Df[which(!Clan_ConAAfreq_Df$Clans %in% diversifiedLBCAclans 
                                                                                  & Clan_ConAAfreq_Df$Clans %in% LBCAclans ),2:21]),
                                                Clan_ConAAfreq_Df$Clan_Conlength[which(!Clan_ConAAfreq_Df$Clans %in% diversifiedLBCAclans
                                                                                       & Clan_ConAAfreq_Df$Clans %in% LBCAclans )] )
#RecentBacteria_Clans_ConAAC <- RecentBacteria_Clans_ConAAC[match(  AA_properties$Letter, names(RecentBacteria_Clans_ConAAC))]
AllLBCA_Clans_ConAAC <- colWeightedMeans(as.matrix(Clan_ConAAfreq_Df[which( Clan_ConAAfreq_Df$Clans %in% LBCAclans ),2:21]),
                                         Clan_ConAAfreq_Df$Clan_Conlength[which( Clan_ConAAfreq_Df$Clans %in% LBCAclans )] )
#AllLBCA_Clans_ConAAC  <- AllLBCA_Clans_ConAAC [match(  AA_properties$Letter, names(AllLBCA_Clans_ConAAC ))]
RecentArchaea_Clans_ConAAC <- colWeightedMeans(as.matrix(Clan_ConAAfreq_Df[which(!Clan_ConAAfreq_Df$Clans %in% diversifiedLACAclans
                                                                                 & Clan_ConAAfreq_Df$Clans %in% LACAclans ),2:21]),
                                               Clan_ConAAfreq_Df$Clan_Conlength[which(!Clan_ConAAfreq_Df$Clans %in% diversifiedLACAclans 
                                                                                      & Clan_ConAAfreq_Df$Clans %in% LACAclans )] )
AllLACA_Clans_ConAAC <- colWeightedMeans(as.matrix(Clan_ConAAfreq_Df[which( Clan_ConAAfreq_Df$Clans %in% LACAclans ),2:21]),
                                         Clan_ConAAfreq_Df$Clan_Conlength[which( Clan_ConAAfreq_Df$Clans %in% LACAclans )] )


#### AAC SE ####
DivLACA_sd <- colWeightedSds(as.matrix(Clan_ConAAfreq_Df[which(Clan_ConAAfreq_Df$Clans %in% diversifiedLACAclans ),2:21]),
                             Clan_ConAAfreq_Df$Clan_Conlength[which(Clan_ConAAfreq_Df$Clans %in% diversifiedLACAclans )] )
DivLACA_se <- DivLACA_sd/sqrt(length(diversifiedLACAclans))
DivLACA_se <- DivLACA_se[match(  AA_properties$Letter, names(DivLACA_se))]
DivLBCA_sd <- colWeightedSds(as.matrix(Clan_ConAAfreq_Df[which(Clan_ConAAfreq_Df$Clans %in% diversifiedLBCAclans ),2:21]),
                             Clan_ConAAfreq_Df$Clan_Conlength[which(Clan_ConAAfreq_Df$Clans %in% diversifiedLBCAclans )] )
DivLBCA_se <- DivLBCA_sd/sqrt(length(diversifiedLBCAclans))
DivLBCA_se <- DivLBCA_se[match(  AA_properties$Letter, names(DivLBCA_se))]
modernBac_Clanssd <- colWeightedSds(as.matrix(Clan_ConAAfreq_Df[which(Clan_ConAAfreq_Df$Clans %in% modernBacteriaclans ),2:21]),
                                    Clan_ConAAfreq_Df$Clan_Conlength[which( Clan_ConAAfreq_Df$Clans %in% modernBacteriaclans)] )
modernBac_Clansse <- modernBac_Clanssd/sqrt(length(modernBacteriaclans))
modernBac_Clansse <- modernBac_Clansse[match(  AA_properties$Letter, names(modernBac_Clansse ))]
modernArc_Clans_ConAAC <- colWeightedMeans(as.matrix(Clan_ConAAfreq_Df[which(Clan_ConAAfreq_Df$Clans %in% modernArchaeaclans ),2:21]),
                                           Clan_ConAAfreq_Df$Clan_Conlength[which( Clan_ConAAfreq_Df$Clans %in% modernArchaeaclans)] )
modernArc_Clanssd <- colWeightedSds(as.matrix(Clan_ConAAfreq_Df[which(Clan_ConAAfreq_Df$Clans %in% modernArchaeaclans ),2:21]),
                                    Clan_ConAAfreq_Df$Clan_Conlength[which( Clan_ConAAfreq_Df$Clans %in% modernArchaeaclans)] )
modernArc_Clansse <- modernArc_Clanssd/sqrt(length(modernArchaeaclans))
modernArc_Clansse <- modernArc_Clansse[match(  AA_properties$Letter, names(modernArc_Clansse))]
NoDivLUCA_sd <- colWeightedSds(as.matrix(Clan_ConAAfreq_Df[which(Clan_ConAAfreq_Df$Clan_ancestor == 'LUCA'),2:21]),
                               Clan_ConAAfreq_Df$Clan_Conlength[which(Clan_ConAAfreq_Df$Clan_ancestor == 'LUCA')] )
NoDivLUCA_sd <- NoDivLUCA_sd/sqrt(length(which(Clan_ConAAfreq_Df$Clan_ancestor == 'LUCA')))

RecentBacteria_se <- colWeightedSds(as.matrix(Clan_ConAAfreq_Df[which(!Clan_ConAAfreq_Df$Clans %in% diversifiedLBCAclans
                     & Clan_ConAAfreq_Df$Clans %in% LBCAclans  ),2:21]),
                     Clan_ConAAfreq_Df$Clan_Conlength[which(!Clan_ConAAfreq_Df$Clans %in% diversifiedLBCAclans
                     & Clan_ConAAfreq_Df$Clans %in% LBCAclans  )] )/sqrt(length(which(!Clan_ConAAfreq_Df$Clans %in% diversifiedLBCAclans
                                                                                                                                            & Clan_ConAAfreq_Df$Clans %in% LBCAclans )))
RecentBacteria_se <- RecentBacteria_se[match(  AA_properties$Letter, names(RecentBacteria_se))]
RecentArchaea_se <- colWeightedSds(as.matrix(Clan_ConAAfreq_Df[which(!Clan_ConAAfreq_Df$Clans %in% diversifiedLACAclans
                                                                     & Clan_ConAAfreq_Df$Clans %in% LACAclans  ),2:21]),
                                   Clan_ConAAfreq_Df$Clan_Conlength[which(!Clan_ConAAfreq_Df$Clans %in% diversifiedLACAclans
                                                                          & Clan_ConAAfreq_Df$Clans %in% LACAclans  )] )/sqrt(length(which(!Clan_ConAAfreq_Df$Clans %in% diversifiedLACAclans
                                                                                                                                           & Clan_ConAAfreq_Df$Clans %in% LACAclans )))
RecentArchaea_se <- RecentArchaea_se[match(  AA_properties$Letter, names(RecentArchaea_se))]
AllLBCA_se <- colWeightedSds(as.matrix(Clan_ConAAfreq_Df[which(Clan_ConAAfreq_Df$Clans %in% LBCAclans  ),2:21]),
                             Clan_ConAAfreq_Df$Clan_Conlength[which( Clan_ConAAfreq_Df$Clans %in% LBCAclans)] )/sqrt(length(which(Clan_ConAAfreq_Df$Clans %in% LBCAclans )))                                                                                                                                    
AllLBCA_se <- AllLBCA_se[match(  AA_properties$Letter, names(AllLBCA_se))]

allLUCAweightedse <- vector()
LUCAweightedse <- vector()
preLUCAweightedse <- vector()
postLUCAweightedse <- vector()
NoDivpostLUCAweightedse <- vector()
DivLACAweightedse <- vector()
DivLBCAweightedse <- vector()
RecentArchaeaweightedse <- vector()
RecentBacteriaweightedse <- vector()
ModernArchaeaweightedse <- vector()
ModernBacteriaweightedse <- vector()
Ancient4PreLACAweightedse  <- vector()
colnames(Clan_ConAAfreq_Df)
for (colnb in 2:21) {
  allLUCAweightedse[colnb] <- weighted_se(Clan_ConAAfreq_Df[,colnb][which(grepl('LUCA',Clan_ConAAfreq_Df$Clan_ancestor))],
                                          Clan_ConAAfreq_Df$Clan_Conlength[which(grepl('LUCA',Clan_ConAAfreq_Df$Clan_ancestor))])}
for (colnb in 2:21) {
  LUCAweightedse[colnb] <- weighted_se(Clan_ConAAfreq_Df[,colnb][which(Clan_ConAAfreq_Df$Clan_ancestor == 'LUCA')],
                                       Clan_ConAAfreq_Df$Clan_Conlength[which(Clan_ConAAfreq_Df$Clan_ancestor == 'LUCA')])}
for (colnb in 2:21){
  postLUCAweightedse[colnb] <- weighted_se(Clan_ConAAfreq_Df[,colnb][which(Clan_ConAAfreq_Df$Clan_ancestor == 'postLUCA')],
                                           Clan_ConAAfreq_Df$Clan_Conlength[which(Clan_ConAAfreq_Df$Clan_ancestor == 'postLUCA')])}
for (colnb in 2:21){
  preLUCAweightedse[colnb] <- weighted_se(Clan_ConAAfreq_Df[,colnb][which(Clan_ConAAfreq_Df$Clan_ancestor == 'preLUCA')],
                                          Clan_ConAAfreq_Df$Clan_Conlength[which(Clan_ConAAfreq_Df$Clan_ancestor == 'preLUCA')])}
for (colnb in 2:21){
  NoDivpostLUCAweightedse[colnb] <- weighted_se(Clan_ConAAfreq_Df[,colnb][which(Clan_ConAAfreq_Df$Clan_ancestor == 'postLUCA' &
                                                                                  !Clan_ConAAfreq_Df$Clans %in% Broaddiversifiedclans )],
                                                Clan_ConAAfreq_Df$Clan_Conlength[which(Clan_ConAAfreq_Df$Clan_ancestor == 'postLUCA' &
                                                                                         !Clan_ConAAfreq_Df$Clans %in% Broaddiversifiedclans)])}
for (colnb in 2:21){
  DivLACAweightedse[colnb] <- weighted_se(Clan_ConAAfreq_Df[,colnb][which(Clan_ConAAfreq_Df$Clans %in% diversifiedLACAclans )],
                                          Clan_ConAAfreq_Df$Clan_Conlength[which(Clan_ConAAfreq_Df$Clans %in% diversifiedLACAclans  )])}
for (colnb in 2:21){
  DivLBCAweightedse[colnb] <- weighted_se(Clan_ConAAfreq_Df[,colnb][which(Clan_ConAAfreq_Df$Clans %in% diversifiedLBCAclans )],
                                          Clan_ConAAfreq_Df$Clan_Conlength[which(Clan_ConAAfreq_Df$Clans %in% diversifiedLBCAclans  )])}
for (colnb in 2:21){
  RecentArchaeaweightedse[colnb] <- weighted_se(Clan_ConAAfreq_Df[,colnb][which(!Clan_ConAAfreq_Df$Clans %in% Broaddiversifiedclans  
                                                                                & Clan_ConAAfreq_Df$Clans %in% LACAclans )],
                                                Clan_ConAAfreq_Df$Clan_Conlength[which(!Clan_ConAAfreq_Df$Clans %in% Broaddiversifiedclans  
                                                                                       & Clan_ConAAfreq_Df$Clans %in% LACAclans  )])}
for (colnb in 2:21){
  RecentBacteriaweightedse[colnb] <- weighted_se(Clan_ConAAfreq_Df[,colnb][which(!Clan_ConAAfreq_Df$Clans %in% Broaddiversifiedclans  
                                                                                 & Clan_ConAAfreq_Df$Clans %in% LBCAclans )],
                                                 Clan_ConAAfreq_Df$Clan_Conlength[which(!Clan_ConAAfreq_Df$Clans %in% Broaddiversifiedclans  
                                                                                        & Clan_ConAAfreq_Df$Clans %in% LBCAclans  )])}
for (colnb in 2:21){
  ModernArchaeaweightedse[colnb] <- weighted_se(Clan_ConAAfreq_Df[,colnb][which(Clan_ConAAfreq_Df$Clans %in% modernArchaeaclans)],
                                                Clan_ConAAfreq_Df$Clan_Conlength[which(Clan_ConAAfreq_Df$Clans %in% modernArchaeaclans )])}
for (colnb in 2:21){
  ModernBacteriaweightedse[colnb] <- weighted_se(Clan_ConAAfreq_Df[,colnb][which(Clan_ConAAfreq_Df$Clans %in% modernBacteriaclans)],
                                                 Clan_ConAAfreq_Df$Clan_Conlength[which(Clan_ConAAfreq_Df$Clans %in% modernBacteriaclans )])}
for (colnb in 2:21){
  Ancient4PreLACAweightedse[colnb] <- weighted_se(Clan_ConAAfreq_Df[,colnb][which(Clan_ConAAfreq_Df$Clans %in% BroaddiversifiedLACAclans[which(MultipleLACAclans  > 1)] )],
                                                  Clan_ConAAfreq_Df$Clan_Conlength[which(Clan_ConAAfreq_Df$Clans %in% BroaddiversifiedLACAclans[which(MultipleLACAclans  > 1)] )])}


### Calculate usages ####
NoDivpostLUCA_Clans_ConAAC <- colWeightedMeans(as.matrix(Clan_ConAAfreq_Df[which(Clan_ConAAfreq_Df$Clan_ancestor == 'postLUCA' &
                               !Clan_ConAAfreq_Df$Clans %in% Broaddiversifiedclans ),2:21]),
                                Clan_ConAAfreq_Df$Clan_Conlength[which(Clan_ConAAfreq_Df$Clan_ancestor == 'postLUCA' &
                               !Clan_ConAAfreq_Df$Clans %in% Broaddiversifiedclans)] )
NoDiv_LUCAusage <- LUCA_Clans_ConAAC/NoDivpostLUCA_Clans_ConAAC
NoDiv_LUCAusage <- NoDiv_LUCAusage[match(  AA_properties$Letter, names(NoDiv_LUCAusage))]
NoDiv_preLUCAusage <- preLUCA_Clans_ConAAC/NoDivpostLUCA_Clans_ConAAC 
NoDiv_preLUCAusage <- NoDiv_preLUCAusage[match(  AA_properties$Letter, names(NoDiv_preLUCAusage))]

Ancient4PreLACAusage <- colWeightedMeans(as.matrix(Clan_ConAAfreq_Df[which(Clan_ConAAfreq_Df$Clans %in% BroaddiversifiedLACAclans[which(MultipleLACAclans  > 1)]),2:21]),
                                         Clan_ConAAfreq_Df$Clan_Conlength[which(Clan_ConAAfreq_Df$Clans %in% BroaddiversifiedLACAclans[which(MultipleLACAclans  > 1)] )] )/ LUCA_Clans_ConAAC
colWeightedSds(as.matrix(Clan_ConAAfreq_Df[which(Clan_ConAAfreq_Df$Clans %in% BroaddiversifiedLACAclans[which(MultipleLACAclans  > 1)]),2:21]),
               Clan_ConAAfreq_Df$Clan_Conlength[which(Clan_ConAAfreq_Df$Clans %in% BroaddiversifiedLACAclans[which(MultipleLACAclans  > 1)] )] )
0.001535329/sqrt(4) # 4 pre-laca clans se NQ.bac
0.001356431/ sqrt(4) # 4 pre-laca clans se NQ.arch


DivLBCAclan_usage <- DivLBCA_ClanAAC/ LUCA_Clans_ConAAC
DivLACAclan_usage <- DivLACA_ClanAAC/ LUCA_Clans_ConAAC
DivLACAclan_usage <- DivLACAclan_usage[match(  AA_properties$Letter, names(DivLACAclan_usage))]
DivLBCAclan_usage <- DivLBCAclan_usage[match(  AA_properties$Letter, names(DivLBCAclan_usage))]


RecentBacusage <- RecentBacteria_Clans_ConAAC/LUCA_Clans_ConAAC
RecentBacusage <- RecentBacusage[match(  AA_properties$Letter, names(RecentBacusage))]
RecentArcusage <- RecentArchaea_Clans_ConAAC/LUCA_Clans_ConAAC
RecentArcusage <- RecentArcusage[match(  AA_properties$Letter, names(RecentArcusage))]

ModernArcusage <- modernArc_Clans_ConAAC/LUCA_Clans_ConAAC
ModernArcusage <- ModernArcusage[match(  AA_properties$Letter, names(ModernArcusage))]
ModernBacusage <- modernBac_Clans_ConAAC/LUCA_Clans_ConAAC
ModernBacusage <- ModernBacusage[match(  AA_properties$Letter, names(ModernBacusage))]


#### standard errors for usages #####
LUCAclanratio_var <-  (LUCAweightedse[-1]^2)/((postLUCA_Clans_ConAAC)^2) +
  (postLUCAweightedse[-1]^2)*((LUCA_Clans_ConAAC)^2)/((postLUCA_Clans_ConAAC)^4)
LUCAclanratio_se <- sqrt(LUCAclanratio_var)
LUCAclanratio_se <- LUCAclanratio_se[match(AA_properties$Letter,names(LUCAclanratio_se))]
preLUCAclanratio_var <-  (preLUCAweightedse[-1]^2)/((postLUCA_Clans_ConAAC)^2) +
  (postLUCAweightedse[-1]^2)*((preLUCA_Clans_ConAAC)^2)/((postLUCA_Clans_ConAAC)^4)
preLUCAclanratio_se <- sqrt(preLUCAclanratio_var)
preLUCAclanratio_se <- preLUCAclanratio_se[match(AA_properties$Letter,names(preLUCAclanratio_se))]

NodivLUCAclanratio_var <-  (LUCAweightedse[-1]^2)/((NoDivpostLUCA_Clans_ConAAC)^2) +
  (NoDivpostLUCAweightedse[-1]^2)*((LUCA_Clans_ConAAC)^2)/((NoDivpostLUCA_Clans_ConAAC)^4)
NodivLUCAclanratio_se <- sqrt(NodivLUCAclanratio_var)
NodivLUCAclanratio_se <- NodivLUCAclanratio_se[match(AA_properties$Letter,names(NodivLUCAclanratio_se))]

NodivpreLUCAclanratio_var <-  (preLUCAweightedse[-1]^2)/((NoDivpostLUCA_Clans_ConAAC)^2) +
  (NoDivpostLUCAweightedse[-1]^2)*((preLUCA_Clans_ConAAC)^2)/((NoDivpostLUCA_Clans_ConAAC)^4)
NodivpreLUCAclanratio_se <- sqrt(NodivpreLUCAclanratio_var)
NodivpreLUCAclanratio_se <- NodivpreLUCAclanratio_se[match(AA_properties$Letter,names(NodivpreLUCAclanratio_se))]


DivLACAclanratio_var <-  (DivLACAweightedse[-1]^2)/((LUCA_Clans_ConAAC)^2) +
  (LUCAweightedse[-1]^2)*((DivLACA_ClanAAC)^2)/((LUCA_Clans_ConAAC)^4)
DivLACAclanratio_se <- sqrt(DivLACAclanratio_var)
DivLACAclanratio_se <- DivLACAclanratio_se[match(AA_properties$Letter,names(DivLACAclanratio_se))]
# SE on the 4 most ancient archaeal clan usages
sqrt((Ancient4PreLACAweightedse[-1]^2)/((LUCA_Clans_ConAAC)^2) +
       (LUCAweightedse[-1]^2)*((colWeightedMeans(as.matrix(Clan_ConAAfreq_Df[which(Clan_ConAAfreq_Df$Clans %in% BroaddiversifiedLACAclans[which(MultipleLACAclans  > 1)]),2:21]),
                                                 Clan_ConAAfreq_Df$Clan_Conlength[which(Clan_ConAAfreq_Df$Clans %in% BroaddiversifiedLACAclans[which(MultipleLACAclans  > 1)] )] ))^2)/((LUCA_Clans_ConAAC)^4))

DivLBCAclanratio_var <- (DivLBCAweightedse[-1]^2)/((LUCA_Clans_ConAAC)^2) +
  (LUCAweightedse[-1]^2)*((DivLBCA_ClanAAC)^2)/((LUCA_Clans_ConAAC)^4)
DivLBCAclanratio_se <- sqrt(DivLBCAclanratio_var)
DivLBCAclanratio_se <- DivLBCAclanratio_se[match(AA_properties$Letter,names(DivLBCAclanratio_se))]

sqrt( (DivLBCAweightedse[-1]^2)/((DivLACA_ClanAAC)^2) +
        (DivLACAweightedse[-1]^2)*((DivLBCA_ClanAAC)^2)/((DivLACA_ClanAAC)^4))
DivLBCALACAusage

RecentArcusage_var <- (RecentArchaeaweightedse[-1]^2)/((LUCA_Clans_ConAAC)^2) +
  (LUCAweightedse[-1]^2)*((RecentArchaea_Clans_ConAAC)^2)/((LUCA_Clans_ConAAC)^4)
RecentArcusage_se <- sqrt(RecentArcusage_var)
RecentArcusage_se <- RecentArcusage_se[match(AA_properties$Letter,names(RecentArcusage_se ))]
RecentBacusage_var <- (RecentBacteriaweightedse[-1]^2)/((LUCA_Clans_ConAAC)^2) +
  (LUCAweightedse[-1]^2)*((RecentBacteria_Clans_ConAAC)^2)/((LUCA_Clans_ConAAC)^4)
RecentBacusage_se <- sqrt(RecentBacusage_var)
RecentBacusage_se <- RecentBacusage_se[match(AA_properties$Letter,names(RecentBacusage_se ))]

ModernArcusage_var <- (ModernArchaeaweightedse[-1]^2)/((LUCA_Clans_ConAAC)^2) +
  (LUCAweightedse[-1]^2)*((modernArc_Clans_ConAAC)^2)/((LUCA_Clans_ConAAC)^4)
ModernArcusage_se <- sqrt(ModernArcusage_var)
ModernArcusage_se <- ModernArcusage_se[match(AA_properties$Letter,names(ModernArcusage_se))]
ModernBacusage_var <- (ModernBacteriaweightedse[-1]^2)/((LUCA_Clans_ConAAC)^2) +
  (LUCAweightedse[-1]^2)*((modernBac_Clans_ConAAC)^2)/((LUCA_Clans_ConAAC)^4)
ModernBacusage_se <- sqrt(ModernBacusage_var)
ModernBacusage_se <- ModernBacusage_se[match(AA_properties$Letter,names(ModernBacusage_se))]



#### Figure 4B and Supp Figure 9-10 ####
# run twice for nq.arch and nq.bac
#laca_lbca_df <- read.csv( 'laca_lbca_df_arc.csv',header = T)
laca_lbca_df <- read.csv( 'laca_lbca_df_bac.csv',header = T)

ggplot( laca_lbca_df,aes(y =as.numeric(DivLACAclan_usage), 
   x = as.numeric(DivLBCAclan_usage), label =X)) + geom_text(color='blue', size = 10) +
  geom_errorbar(aes(ymin=as.numeric(DivLACAclan_usage)-DivLACAclanratio_se, ymax=as.numeric(DivLACAclan_usage)+DivLACAclanratio_se)) +
  geom_errorbarh(aes(xmin=as.numeric(DivLBCAclan_usage)-DivLBCAclanratio_se, xmax=as.numeric(DivLBCAclan_usage)+DivLBCAclanratio_se)) +
  ylab('Pre-LACA/LUCA') + xlab('Pre-LBCA/LUCA')  +
  theme(legend.position = 'none') +
  geom_hline(yintercept = 1, linetype = 2) + geom_vline(xintercept = 1,linetype = 2) +
  theme(axis.text=element_text(size=24),axis.title=element_text(size=32,face="bold"),
        legend.text=element_text(size=25), legend.title = element_text(size=25,face="bold")) +
  coord_trans(x = "log10", y = "log10") +
  scale_x_continuous(breaks = c(0.6, 0.9 , 1.2, 1.5),
                     minor_breaks = c(0.7,0.8,1,1.1,1.3,1.4, 1.6,1.7),
                     guide = guide_axis(minor.ticks = TRUE),
                     expand = expansion(mult = c(0.4, 0.2))) +
  scale_y_continuous(breaks = c(0.6, 0.9 , 1.2, 1.5),
                     minor_breaks = c(0.7,0.8,1,1.1,1.3,1.4,1.6,1.7),
                     guide = guide_axis(minor.ticks = TRUE),
                     expand = expansion(mult = c(0.1, 0.1)))

ggplot( laca_lbca_df,aes(y = as.numeric(ModernArcusage), 
                         x = as.numeric(ModernBacusage), label =X)) + geom_text(color='blue', size = 10) +
  geom_errorbar(aes(ymin=as.numeric(ModernArcusage)-ModernArcusage_se, ymax=as.numeric(ModernArcusage)+ModernArcusage_se)) +
  geom_errorbarh(aes(xmin=as.numeric(ModernBacusage)- ModernBacusage_se, xmax=as.numeric(ModernBacusage)+ ModernBacusage_se)) +
  ylab('Modern Archaea/LUCA') + xlab('Modern Bacteria/LUCA')  +
  theme(legend.position = 'none') +  
  #xlim(0.67 ,1.83) + 
  geom_hline(yintercept = 1, linetype = 2) + geom_vline(xintercept = 1,linetype = 2) +
  theme(axis.text=element_text(size=24),axis.title=element_text(size=32,face="bold"),
        legend.text=element_text(size=25), legend.title = element_text(size=25,face="bold")) +
  coord_trans(x = "log10", y = "log10") +
  scale_x_continuous(breaks = c(0.6, 0.9 , 1.2, 1.5),
                     minor_breaks = c(0.7,0.8,1,1.1,1.3,1.4, 1.5),
                     guide = guide_axis(minor.ticks = TRUE),
                     expand = expansion(mult = c(0.4, 0.2))) +
  scale_y_continuous(breaks = c(0.6, 0.9 , 1.2, 1.5),
                     minor_breaks = c(0.7,0.8,1,1.1,1.3,1.4, 1.5),
                     guide = guide_axis(minor.ticks = TRUE),
                     expand = expansion(mult = c(0.4, 0.2)),
                     limits = c(0.6,1.5))


#### Supplementary Figure 11A-D; run twice to generate figures with nq.bac and nq.arch ####
W_sitesproblist <- readRDS("W_sitesproblistPfamIDs0.4_arc.rds")
#W_sitesproblist <- readRDS("W_sitesproblistPfamIDs0.4_bac.rds")
W_preLACA <- unlist(replace(W_sitesproblist[which(names(W_sitesproblist) %in% AncientPostLUCA$PFAM_IDs[which(AncientPostLUCA$classified_ancestor == 'preLACA')])],
                            lengths(W_sitesproblist[which(names(W_sitesproblist) %in% AncientPostLUCA$PFAM_IDs[which(AncientPostLUCA$classified_ancestor == 'preLACA')])]) == 0, 0))
names(W_preLACA ) <- str_extract(names(W_preLACA ),'.......')

W_preLBCA <- unlist(replace(W_sitesproblist[which(names(W_sitesproblist) %in% AncientPostLUCA$PFAM_IDs[which(AncientPostLUCA$classified_ancestor == 'preLBCA')])],
                            lengths(W_sitesproblist[which(names(W_sitesproblist) %in% AncientPostLUCA$PFAM_IDs[which(AncientPostLUCA$classified_ancestor == 'preLBCA')])]) == 0, 0))
length(which(W_preLBCA < 0.7))/ length(W_preLBCA) #70%

PreLACAecdf_fun <- ecdf(W_preLACA )
PreLACA_contr <- PreLACAecdf_fun(W_preLACA )*W_preLACA 
PreLBCAecdf_fun <- ecdf(W_preLBCA )
PreLBCA_contr <- PreLBCAecdf_fun(W_preLBCA )*W_preLBCA 
W_ecdf <- bind_rows(
  data.frame(val = W_preLACA, contribution= PreLACA_contr,group = "preLACA"),
  data.frame(val = W_preLBCA, contribution= PreLBCA_contr, group = "preLBCA"))

Y_sitesproblist <- readRDS("Y_sitesproblistPfamIDs0.4_arc.rds")
#Y_sitesproblist <- readRDS("Y_sitesproblistPfamIDs0.4_bac.rds")
Y_preLACA <- unlist(replace(Y_sitesproblist[which(names(Y_sitesproblist) %in% AncientPostLUCA$PFAM_IDs[which(AncientPostLUCA$classified_ancestor == 'preLACA')])],
                            lengths(Y_sitesproblist[which(names(Y_sitesproblist) %in% AncientPostLUCA$PFAM_IDs[which(AncientPostLUCA$classified_ancestor == 'preLACA')])]) == 0, 0))
names(Y_preLACA ) <- str_extract(names(Y_preLACA ),'.......')
Y_preLBCA <- unlist(replace(Y_sitesproblist[which(names(Y_sitesproblist) %in% AncientPostLUCA$PFAM_IDs[which(AncientPostLUCA$classified_ancestor == 'preLBCA')])],
                            lengths(Y_sitesproblist[which(names(Y_sitesproblist) %in% AncientPostLUCA$PFAM_IDs[which(AncientPostLUCA$classified_ancestor == 'preLBCA')])]) == 0, 0))
PreLACAecdf_fun <- ecdf(Y_preLACA )
PreLACA_contr <- PreLACAecdf_fun(Y_preLACA )*Y_preLACA 
PreLBCAecdf_fun <- ecdf(Y_preLBCA )
PreLBCA_contr <- PreLBCAecdf_fun(Y_preLBCA )*Y_preLBCA 
Y_ecdf <- bind_rows(
  data.frame(val = Y_preLACA, contribution= PreLACA_contr,group = "preLACA"),
  data.frame(val = Y_preLBCA, contribution= PreLBCA_contr, group = "preLBCA"))

grid.arrange(ggplot(W_ecdf, aes(x=val, y= contribution,color = group)) + geom_step() +
               xlab('Ancestral W confidence') +  ylab('Cumulative Contribution') +
               theme(axis.text=element_text(size=24),axis.title=element_text(size=32,face="bold"),
                     legend.text=element_text(size=25), legend.title = element_text(size=25,face="bold"),
                     legend.position = "inside", legend.position.inside = c(0.25,0.7)) ,
             ggplot(Y_ecdf, aes(x=val, y= contribution,color = group)) + geom_step() +
               xlab('Ancestral Y confidence') +  ylab('Cumulative Contribution') +
               theme(axis.text=element_text(size=24),axis.title=element_text(size=32,face="bold"),
                     legend.text=element_text(size=25), legend.title = element_text(size=25,face="bold"),
                     legend.position = "inside", legend.position.inside = c(0.25,0.7)) , ncol =2)


### Supp figure 11E ####
# RUN twice with nq.bac and nq.arch
W_freq_df <- data.frame(preLUCA_Clans_ConAAC[match(  AA_properties$Letter, names(preLUCA_Clans_ConAAC))],
                        LUCA_Clans_ConAAC[match(  AA_properties$Letter, names(LUCA_Clans_ConAAC))],
                        DivLBCA_ClanAAC[match(  AA_properties$Letter, names(DivLBCA_ClanAAC))],
                        DivLACA_ClanAAC[match(  AA_properties$Letter, names(DivLACA_ClanAAC))], 
                        RecentBacteria_Clans_ConAAC[match(  AA_properties$Letter, names(RecentBacteria_Clans_ConAAC))], 
                        RecentArchaea_Clans_ConAAC[match(  AA_properties$Letter, names(RecentArchaea_Clans_ConAAC))],
                        OrderedmodernBac_Clans_ConAAC, OderedmodernArc_Clans_ConAAC, 
                        preLUCAweightedse[-1][match(  AA_properties$Letter, names(preLUCA_Clans_ConAAC))] ,
                        LUCAweightedse[-1][match(  AA_properties$Letter, names(LUCA_Clans_ConAAC))] ,
                        DivLBCA_se,DivLACA_se,RecentBacteria_se,RecentArchaea_se,modernBac_Clansse,modernArc_Clansse)
W_freq_df <- data.frame(as.numeric(W_freq_df[19,1:8]), as.numeric(W_freq_df[19,9:16]))
rownames(W_freq_df) <- c('PreLUCA','LUCA','PreLBCA','PreLACA','LBCA','LACA','Modern Bacteria','Modern Archaea')
colnames(W_freq_df) <- c('W','SE')
W_freq_df_bac <- W_freq_df 

#W_freq_df_ArcBac <- rbind(W_freq_df_bac,W_freq_df_arc)
#W_freq_df_ArcBac$age <- rep(rownames(W_freq_df_ArcBac)[1:8],2)
#W_freq_df_ArcBac$group <- c(rep('NQ.bac',8), rep('NQ.arch',8))
W_freq_df_ArcBac <- read.csv('W_freq_df_ArcBac.csv', header = T)

ggplot(W_freq_df_ArcBac, aes(x = fct_inorder(age), y = W, fill = fct_inorder(group))) +
  geom_bar(stat='identity', position='dodge') + geom_errorbar(aes(ymin = W - SE, ymax = W + SE), width = 0.2, position=position_dodge(0.9) ) + 
  xlab(NULL) + ylab('Ancestrally reconstructed W frequency') +     theme(legend.position = 'none') +
  theme(axis.text=element_text(size=20),axis.title=element_text(size=30,face="bold"),
        legend.text=element_text(size=20), legend.title = element_text(size=20,face="bold"))

## file contains number of 100% conserved W labeled in black on top of bars
conserved_W <- read.csv('Conserved_W.csv', header = T)

