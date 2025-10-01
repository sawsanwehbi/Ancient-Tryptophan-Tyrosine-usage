
# Load the ape package
library(ape)
library(adephylo)
library(protr)
library(seqinr)
library(broom)
library(ggrepel)
library(stringr)
library(ggplot2)
library(diagis)
library(phytools)
library(dplyr)
library(matrixStats)
library(deming)
library(Biostrings)
library("thacklr")
library("pwalign")
library(protr)
library(reshape2)


setwd("~/Desktop/Tryptophan paper")
Tol_tree <- read.tree('Fixed_ToL.newick')
Distance2aaRS <- read.csv('AARSprotoenzyme_LUCA.csv', header = T)
BiosynReq <- read.csv('../PFAM Trees/BiosyntheticRequirement_AA.csv', header = T)
AA_properties <- read.csv('../PFAM Trees/AminoAcid_properties.csv', header = T)
#NoDivAA_properties <- read.csv( '../Tryptophan paper/NoDivAA_properties.csv', header= T)
NoDivAA_properties <- read.csv( '../Tryptophan paper/StrictNoDivAA_properties.csv', header= T)
colnames(NoDivAA_properties)[1] <- 'Letter'
NoDivAA_properties <- inner_join(AA_properties , NoDivAA_properties, by= 'Letter')
NQbac <- read.csv('../PFAM Trees/NQbac.csv', header = T)


##### Protozyme distance analysis ####
#Distance2aaRS <- Distance2aaRS[-22,]
#Distance2aaRS$WehbiOrder <- rank(-Distance2aaRS$luca)
#Distance2aaRS$WehbiOrder_ties <- rank(-Distance2aaRS$luca, ties.method = "min")
cor.test(Distance2aaRS$distanceofAARSfromprotozyme, Distance2aaRS$WehbiOrder, method='spearman', exact = F)
cor.test(Distance2aaRS$distanceofAARSfromprotozyme, Distance2aaRS$WehbiOrder_ties, method='spearman', exact = F)
cor.test(Distance2aaRS$distanceofAARSfromprotozyme, Distance2aaRS$Trifonov2000order_ties, method='spearman', exact = F)


gridExtra::grid.arrange(
  ggplot(Distance2aaRS, aes(x = as.numeric(WehbiOrder_ties),y = as.numeric(distanceofAARSfromprotozyme),
    label = amino.acid)) +  geom_text( size = 10) +  ylab('Distance of aaRS from protozyme') +
    xlab('Wehbi 2024 rank order') +   annotate(geom="text", x=8, y=6.7, label=paste0("Spearman's rho= -0.56"), color="red", size = 11) +
    annotate(geom="text", x=8, y=6.3, label="p = 0.005", color="red", size = 11) + 
    geom_smooth(method='lm', se =T) + theme(axis.text=element_text(size=20),
    axis.title=element_text(size=30,face="bold")) + expand_limits(x = 0)  ,
    ggplot(Distance2aaRS, aes(x = as.numeric(Trifonov2000order_ties),y = as.numeric(distanceofAARSfromprotozyme),
    label = amino.acid)) +  geom_text( size = 10) +  ylab('Distance of aaRS from protozyme') +
    xlab('Trifonov 2000 rank order') +   annotate(geom="text", x=8, y=2.5, label=paste0("Spearman's rho= -0.49"), color="red", size = 11) +
    annotate(geom="text", x=8, y=2.1, label="p = 0.019", color="red", size = 11) + geom_smooth(method='lm',se=T) +
    theme(axis.text=element_text(size=20), axis.title=element_text(size=30,face="bold"))  + 
    expand_limits(x = 0) , ncol = 2)


NoDivmweights_model <- lm( NoDivAA_properties$NoDiv_LUCAusage ~ NoDivAA_properties$MolecularWeightsDa,
weights = 1/(NoDivAA_properties$NodivLUCAclanratio_se^2))
summary(NoDivmweights_model)
NoDivmweights_model_confinterval <- broom::augment(NoDivmweights_model, interval="confidence")
NoDivmweights_modelclass1 <- lm( NoDivAA_properties$NoDiv_LUCAusage[which(NoDivAA_properties$aaRS_Class == 1 | NoDivAA_properties$aaRS_Class == '1_2' )] ~ 
               NoDivAA_properties$MolecularWeightsDa[which(NoDivAA_properties$aaRS_Class == 1 | NoDivAA_properties$aaRS_Class == '1_2')],
    weights = 1/(NoDivAA_properties$NodivLUCAclanratio_se[which(NoDivAA_properties$aaRS_Class == 1 | NoDivAA_properties$aaRS_Class == '1_2')]^2))
NoDivmweights_modelclass2 <- lm( NoDivAA_properties$NoDiv_LUCAusage[which(NoDivAA_properties$aaRS_Class == 2 | NoDivAA_properties$aaRS_Class == '1_2')] ~ 
             NoDivAA_properties$MolecularWeightsDa[which(NoDivAA_properties$aaRS_Class == 2 | NoDivAA_properties$aaRS_Class == '1_2')],
            weights = 1/(NoDivAA_properties$NodivLUCAclanratio_se[which(NoDivAA_properties$aaRS_Class == 2 | NoDivAA_properties$aaRS_Class == '1_2')]^2))
NoDivmweights_model_confintervalclass1 <- broom::augment(NoDivmweights_modelclass1, interval="confidence")
NoDivmweights_model_confintervalclass2 <- broom::augment(NoDivmweights_modelclass2, interval="confidence")
#850x760
NoDivLUCA_plot <- ggplot(NoDivAA_properties, aes(y = as.numeric(NoDiv_LUCAusage), 
  x = as.numeric(MolecularWeightsDa), label = Letter)) + 
  geom_text(aes(colour = factor(aaRS_Class)), size = 16) + 
  labs(color='aaRS Class') + 
  xlab('Molecular Weight (Da)') + ylab('LUCA clan usage')  + 
  theme(legend.position="none") +
  annotate(geom="text", y=0.6, x=100, label="Class II", color="#619CFF", size = 12) +
  annotate(geom="text", y=0.64, x=100, label="Class I", color="#F8766D", size = 12) +
  geom_errorbar(aes(ymin=as.numeric(NoDiv_LUCAusage)-NodivLUCAclanratio_se, ymax=as.numeric(NoDiv_LUCAusage)+NodivLUCAclanratio_se)) +
  theme(axis.text=element_text(size=24),axis.title=element_text(size=32,face="bold"),
        legend.text=element_text(size=14), legend.title = element_text(size=16,face="bold")) + 
  #geom_line(aes(y = predict(NoDivmweights_model)), linewidth = 1,color = 'black') +
  geom_line(data= NoDivAA_properties[which(NoDivAA_properties$aaRS_Class==1 | NoDivAA_properties$aaRS_Class == '1_2'),] ,aes(y = predict(NoDivmweights_modelclass1)), linewidth = 1,color = "#F8766D") +
  geom_line(data= NoDivAA_properties[which(NoDivAA_properties$aaRS_Class==2 | NoDivAA_properties$aaRS_Class == '1_2'),], aes(y = predict(NoDivmweights_modelclass2)), linewidth = 1,color = "#619CFF") +
  annotate(geom="text", x=120, y=0.77, label=paste0(paste0('Weighted R',supsc('2')), "= 0.52"), color= "#619CFF", size = 14) +
  annotate(geom="text", x=120, y=0.72, label="p = 0.018", color= "#619CFF", size = 14) +
  annotate(geom="text", x=160, y=1.2, label=paste0(paste0('Weighted R',supsc('2')), "= 0.57"), color= "#F8766D", size = 14) +
  annotate(geom="text", x=160, y=1.15, label="p = 0.007", color= "#F8766D", size = 14) +
  geom_ribbon(data=NoDivAA_properties[which(NoDivAA_properties$aaRS_Class==1 | NoDivAA_properties$aaRS_Class == '1_2'),] , aes(ymin=NoDivmweights_model_confintervalclass1$.lower, ymax=NoDivmweights_model_confintervalclass1$.upper), colour="#F8766D", alpha=0.3) +
  geom_ribbon(data=NoDivAA_properties[which(NoDivAA_properties$aaRS_Class==2 | NoDivAA_properties$aaRS_Class == '1_2'),] , aes(ymin=NoDivmweights_model_confintervalclass2$.lower, ymax=NoDivmweights_model_confintervalclass2$.upper), colour="#619CFF", alpha=0.3)

summary(NoDivmweights_premodel)
NoDivmweights_premodel <- lm( NoDivAA_properties$NoDiv_preLUCAusage ~ NoDivAA_properties$MolecularWeightsDa,
                           weights = 1/(NoDivAA_properties$NodivpreLUCAclanratio_se^2))
NoDivmweights_premodel_confinterval <- broom::augment(NoDivmweights_premodel, interval="confidence")
NoDivmweights_premodelclass1 <- lm( NoDivAA_properties$NoDiv_preLUCAusage[which(NoDivAA_properties$aaRS_Class == 1)] ~ 
                                      NoDivAA_properties$MolecularWeightsDa[which(NoDivAA_properties$aaRS_Class == 1)],
                                 weights = 1/(NoDivAA_properties$NodivpreLUCAclanratio_se[which(NoDivAA_properties$aaRS_Class == 1)]^2))
NoDivmweights_premodelclass2 <- lm( NoDivAA_properties$NoDiv_preLUCAusage[which(NoDivAA_properties$aaRS_Class == 2)] ~ 
                                      NoDivAA_properties$MolecularWeightsDa[which(NoDivAA_properties$aaRS_Class == 2)],
                                 weights = 1/(NoDivAA_properties$NodivpreLUCAclanratio_se[which(NoDivAA_properties$aaRS_Class == 2)]^2))
NoDivmweights_premodel_confintervalclass1 <- broom::augment(NoDivmweights_premodelclass1, interval="confidence")
NoDivmweights_premodel_confintervalclass2 <- broom::augment(NoDivmweights_premodelclass2, interval="confidence")
NoDivpreLUCA_plot <- ggplot(NoDivAA_properties, aes(y = as.numeric(NoDiv_preLUCAusage), 
  x = as.numeric(MolecularWeightsDa), label = Letter)) + 
  xlab('Molecular Weight (Da)') + ylab('preLUCA clan usage')  + 
  theme(legend.position="none") + geom_text(aes(colour = factor(aaRS_Class)), size = 16) + 
  labs(color='aaRS Class') + 
  annotate(geom="text", y=1.2, x=180, label="Class II", color="#619CFF", size = 12) +
  annotate(geom="text", y=1.15, x=180, label="Class I", color="#F8766D", size = 12) +
  #annotate(geom="text", y=1.1, x=180, label="Class I & II", color="#00BA38", size = 12) +
  geom_errorbar(aes(ymin=as.numeric(NoDiv_preLUCAusage)-NodivpreLUCAclanratio_se, ymax=as.numeric(NoDiv_preLUCAusage)+NodivpreLUCAclanratio_se)) +
  theme(axis.text=element_text(size=24),axis.title=element_text(size=32,face="bold"),
        legend.text=element_text(size=14), legend.title = element_text(size=16,face="bold")) + 
  geom_line(aes(y = predict(NoDivmweights_premodel)), linewidth = 1,color = 'black') +
  annotate(geom="text", x=140, y=0.73, label=paste0(paste0('Weighted R',supsc('2')), "= 0.37"), color= "black", size = 14) +
  annotate(geom="text", x=140, y=0.67, label="p = 0.004", color= "black", size = 14) +
  #annotate(geom="text", x=70, y=1.2, label="b)", color="black", size = 16) +
  geom_ribbon(aes(ymin=NoDivmweights_premodel_confinterval$.lower, ymax=NoDivmweights_premodel_confinterval$.upper), colour=NA, alpha=0.3)


NoDivConclan_wls_model <- lm(NoDivAA_properties$Filtered_avg_2000 ~ NoDivAA_properties$NoDiv_LUCAusage , 
                        weights = 1 / (AA_properties$Filtered_sd_2000^2))
NoDivConclan_wls_confinterval <- broom::augment(NoDivConclan_wls_model , interval="confidence")
summary(NoDivConclan_wls_model )
NoDivLUCAvsTrifonov_plot <- ggplot(NoDivAA_properties, aes(y = as.numeric(Filtered_avg_2000 ), 
x = as.numeric(NoDiv_LUCAusage), label = Letter)) + 
  ylab('Trifonov (2000) order') +
  xlab('LUCA clan usage') + 
  geom_text(color='blue', size = 12) +
  theme(axis.text=element_text(size=24),axis.title=element_text(size=32,face="bold"),
        legend.text=element_text(size=27), legend.title = element_text(size=16,face="bold")) + 
  scale_fill_brewer(palette="Dark2") +  theme(legend.position = 'none') +
  guides(color = guide_legend(override.aes = list(size = 16))) + 
  geom_line(aes(y = predict(NoDivConclan_wls_model)), linewidth = 1,color = 'black') +
  annotate(geom="text", y=26, x=0.9, label=paste0("Weighted ",paste0('R',supsc('2')), "= 0.36"), color="black", size = 14) +
  annotate(geom="text", y=23, x=0.9, label="p = 0.005", color="black", size = 14) +
  #annotate(geom="text", y=26, x=0.63, label="c)", color="black", size = 16) +
  geom_ribbon(aes(ymin=NoDivConclan_wls_confinterval$.lower, ymax=NoDivConclan_wls_confinterval$.upper), colour=NA, alpha=0.3)


NoDivLACAvsLBCA_plot <- ggplot(NoDivAA_properties, aes(y = as.numeric(DivLACAclan_usage), 
  x = as.numeric(DivLBCAclan_usage), label = Letter)) + 
  #geom_errorbar(aes(ymin=as.numeric(DivLACAclan_usage)-DivLACAclanratio_se, ymax=as.numeric(DivLACAclan_usage)+DivLACAclanratio_se)) +
  #geom_errorbarh(aes(xmin=as.numeric(DivLBCAclan_usage)-DivLBCAclanratio_se, xmax=as.numeric(DivLBCAclan_usage)+DivLBCAclanratio_se)) +
  geom_text(aes(colour = factor(middle_nt_letter)), size = 8) + 
  labs(color='Middle nucleotide') + 
  ylab('Ancient Archaea usage') + xlab('Ancient Bacteria usage')  +
  theme(legend.position = 'none') + xlim(0.65 ,1.3) + ylim(0.3,1.6) +
  #scale_fill_brewer(palette="Dark2") +
  geom_hline(yintercept = 1, linetype = 2) + geom_vline(xintercept = 1,linetype = 2) +
  #annotate(geom="text", y=1.4, x=0.7, label="d)", color="black", size = 16) +
  scale_y_continuous(breaks = c(0.6,1,1.4)) +
  theme(axis.text=element_text(size=24),axis.title=element_text(size=32,face="bold"),
  legend.text=element_text(size=25), legend.title = element_text(size=25,face="bold")) +
  #annotate(geom="text", y=0.6, x=0.8, label="Subclass a", color="#619CFF", size = 12) +
  #annotate(geom="text", y=0.5, x=0.8, label="Subclass c", color="#F8766D", size = 12) +
  #annotate(geom="text", y=0.4, x=0.8, label="Subclass b", color="#00BA38", size = 12)
  annotate(geom="point", y=DivLACAmiddlent_usage[1], x=DivLBCAmiddlent_usage[1], color="#F8766D", size =8) +
  annotate(geom="point", y=DivLACAmiddlent_usage[2], x=DivLBCAmiddlent_usage[2], color="#E76BF3", size = 8) +
  annotate(geom="point", y=DivLACAmiddlent_usage[3], x=DivLBCAmiddlent_usage[3], color="#00B6EB", size = 8) +
  annotate(geom="point", y=DivLACAmiddlent_usage[4], x=DivLBCAmiddlent_usage[4], label=rownames(DivLBCA_LACAmiddlent_usageDf)[4], color="#B79F00", size = 8) 

## 2000x1800
grid.arrange(NoDivLUCA_plot, NoDivpreLUCA_plot,NoDivLUCAvsTrifonov_plot, NoDivLACAvsLBCA_plot , ncol =2 )
# 3000x1000 or 2700x700
grid.arrange(NoDivLUCA_plot, NoDivpreLUCA_plot,NoDivLUCAvsTrifonov_plot,  ncol =3 )

ggplot(SinglevsMultiNoDivAA_properties, aes(y = as.numeric(strict_DivLACAclan_usage), 
   x = as.numeric(strict_DivLBCAclan_usage), label = Letter)) +
  geom_errorbar(aes(ymin=as.numeric(strict_DivLACAclan_usage)-strict_DivLACAclanratio_se, ymax=as.numeric(strict_DivLACAclan_usage)+strict_DivLACAclanratio_se)) +
  geom_errorbarh(aes(xmin=as.numeric(strict_DivLBCAclan_usage)-strict_DivLBCAclanratio_se, xmax=as.numeric(strict_DivLBCAclan_usage)+strict_DivLBCAclanratio_se)) +
  geom_text(color='blue', size = 10) + ylab('Ancient Archaea usage') + xlab('Ancient Bacteria usage')  +
  theme(legend.position = 'none') + xlim(0 ,1.68) + ylim(0,1.68) +
  geom_hline(yintercept = 1, linetype = 2) + geom_vline(xintercept = 1,linetype = 2) +
  #annotate(geom="text", y=1.4, x=0.7, label="d)", color="black", size = 16) +
  theme(axis.text=element_text(size=24),axis.title=element_text(size=32,face="bold"),
  legend.text=element_text(size=25), legend.title = element_text(size=25,face="bold"))

### testing differences in single vs multi laca/lbca ancient clans ####
grid.arrange(ggplot(SinglevsMultiNoDivAA_properties, aes(y = as.numeric(strict_DivLACAclan_usage), 
  x = as.numeric(DivLACAclan_usage), label = Amino.Acid)) + geom_text( size = 8) +
  geom_errorbar(aes(ymin=as.numeric(strict_DivLACAclan_usage)-strict_DivLACAclanratio_se, ymax=as.numeric(strict_DivLACAclan_usage)+strict_DivLACAclanratio_se)) +
  geom_errorbarh(aes(xmin=as.numeric(DivLACAclan_usage)-DivLACAclanratio_se, xmax=as.numeric(DivLACAclan_usage)+DivLACAclanratio_se)) +
    ylab('Multi-LACA clan usage') + xlab('Single-LACA clan usage')  +
  geom_abline(aes(slope=1, intercept = 0),linewidth = 1,color = 'darkred') + 
  geom_hline(yintercept = 1, linetype = 2) + geom_vline(xintercept = 1,linetype = 2) +
  theme(axis.text=element_text(size=24),axis.title=element_text(size=32,face="bold"),
        legend.text=element_text(size=25), legend.title = element_text(size=25,face="bold")),
  
ggplot(SinglevsMultiNoDivAA_properties, aes(y = as.numeric(strict_DivLBCAclan_usage), 
  x = as.numeric(DivLBCAclan_usage), label = Amino.Acid)) + geom_text( size = 8) +
  geom_errorbar(aes(ymin=as.numeric(strict_DivLBCAclan_usage)-strict_DivLBCAclanratio_se, ymax=as.numeric(strict_DivLBCAclan_usage)+strict_DivLBCAclanratio_se)) +
  geom_errorbarh(aes(xmin=as.numeric(DivLBCAclan_usage)-DivLBCAclanratio_se, xmax=as.numeric(DivLBCAclan_usage)+DivLBCAclanratio_se)) +
  ylab('Multi-LBCA clan usage') + xlab('Single-LBCA clan usage')  +
  geom_abline(aes(slope=1, intercept = 0),linewidth = 1,color = 'darkred') + 
  geom_hline(yintercept = 1, linetype = 2) + geom_vline(xintercept = 1,linetype = 2) +
  theme(axis.text=element_text(size=24),axis.title=element_text(size=32,face="bold"),
        legend.text=element_text(size=25), legend.title = element_text(size=25,face="bold")), ncol=2)


##### Plotting middle nucleotide usages ####
DivLBCA_LACAmiddlent_usageDf <- read.csv('DivLBCA_LACAmiddlent_usage.csv', header = T)
rownames(DivLBCA_LACAmiddlent_usageDf) <- c('Adenine','Uracil','Guanine', 'Cytosine')
ggplot(DivLBCA_LACAmiddlent_usageDf , aes(y = as.numeric(DivLACAmiddlent_usage), 
  x = as.numeric(DivLBCAmiddlent_usage), label = rownames(DivLBCA_LACAmiddlent_usageDf))) + 
  annotate(geom="text", y=DivLACAmiddlent_usage[1], x=DivLBCAmiddlent_usage[1], label=rownames(DivLBCA_LACAmiddlent_usageDf)[1], color="#F8766D", size = 12) +
  annotate(geom="text", y=DivLACAmiddlent_usage[2], x=DivLBCAmiddlent_usage[2], label=rownames(DivLBCA_LACAmiddlent_usageDf)[2], color="#E76BF3", size = 12) +
  annotate(geom="text", y=DivLACAmiddlent_usage[3], x=DivLBCAmiddlent_usage[3], label=rownames(DivLBCA_LACAmiddlent_usageDf)[3], color="#00B6EB", size = 12) +
  annotate(geom="text", y=DivLACAmiddlent_usage[4], x=DivLBCAmiddlent_usage[4], label=rownames(DivLBCA_LACAmiddlent_usageDf)[4], color="#B79F00", size = 12) +
  ylab('Ancient Archaea usage') + xlab('Ancient Bacteria usage')  + xlim(0.88 ,1.1) + ylim(0.88,1.1) +
  geom_hline(yintercept = 1, linetype = 2) + geom_vline(xintercept = 1,linetype = 2) +
  #scale_y_continuous(breaks = c(0.6,0.8,1,1.05)) + scale_x_continuous(breaks = c(0.6,0.8,1,1.05)) +
 theme(axis.text=element_text(size=24),axis.title=element_text(size=32,face="bold"),
        legend.text=element_text(size=25), legend.title = element_text(size=25,face="bold")) 


#### Updated lUCA usage vs NE and NSA ####
NoDivNE_model <- lm(BiosynReq$NoDiv_LUCAusage ~ BiosynReq$Nb_BiosyntEnzymes,
weights = 1/(BiosynReq$NoDiv_LUCAusage^2))
NoDivNE_confinterval <- broom::augment(NoDivNE_model, interval="confidence")
summary(NoDivNE_model)
NoDivLUCA_NE_plot <- ggplot(BiosynReq, aes(y = as.numeric(NoDiv_LUCAusage), 
  x = as.numeric(Nb_BiosyntEnzymes), label = Amino.Acid)) + 
  xlab('Number of biosynthetic enzymes') + ylab('LUCA clan usage')  + 
  theme(legend.position="none") + geom_text(aes(colour = factor(aaRS_Class)), size = 16) + 
  labs(color='aaRS Class') + 
  annotate(geom="text", y=1.2, x=23, label="Class II", color="#619CFF", size = 12) +
  annotate(geom="text", y=1.15, x=23, label="Class I", color="#F8766D", size = 12) +
  #annotate(geom="text", y=1.1, x=180, label="Class I & II", color="#00BA38", size = 12) +
  geom_errorbar(aes(ymin=as.numeric(NoDiv_LUCAusage)-NodivLUCAclanratio_se, ymax=as.numeric(NoDiv_LUCAusage)+NodivLUCAclanratio_se)) +
  theme(axis.text=element_text(size=24),axis.title=element_text(size=32,face="bold"),
        legend.text=element_text(size=14), legend.title = element_text(size=16,face="bold")) + 
  geom_line(aes(y = predict(NoDivNE_model)), linewidth = 1,color = 'black') +
  annotate(geom="text", x=14, y=0.73, label=paste0(paste0('Weighted R',supsc('2')), "= 0.23"), color= "black", size = 14) +
  annotate(geom="text", x=14, y=0.67, label="p = 0.03", color= "black", size = 14) +
  #annotate(geom="text", x=70, y=1.2, label="b)", color="black", size = 16) +
  geom_ribbon(aes(ymin=NoDivNE_confinterval$.lower, ymax=NoDivNE_confinterval$.upper), colour=NA, alpha=0.3)


NoDivNSA_model <- lm( NoDivAA_properties$NoDiv_LUCAusage ~ NoDivAA_properties$Nb_nonH_sidechain,
                     weights = 1/(NoDivAA_properties$NoDiv_LUCAusage^2))
NoDivNSA_confinterval <- broom::augment(NoDivNSA_model, interval="confidence")
NoDivNSA_modelclass1 <- lm( NoDivAA_properties$NoDiv_LUCAusage[which(NoDivAA_properties$aaRS_Class == 1 | NoDivAA_properties$aaRS_Class == "1_2")] ~ 
                              NoDivAA_properties$Nb_nonH_sidechain[which(NoDivAA_properties$aaRS_Class == 1 | NoDivAA_properties$aaRS_Class == "1_2")],
                              weights = 1/(NoDivAA_properties$NodivLUCAclanratio_se[which(NoDivAA_properties$aaRS_Class == 1 | NoDivAA_properties$aaRS_Class == "1_2")]^2))
NoDivNSA_modelclass2 <- lm( NoDivAA_properties$NoDiv_LUCAusage[which(NoDivAA_properties$aaRS_Class == 2 | NoDivAA_properties$aaRS_Class == "1_2")] ~ 
                              NoDivAA_properties$Nb_nonH_sidechain[which(NoDivAA_properties$aaRS_Class == 2 | NoDivAA_properties$aaRS_Class == "1_2")],
                        weights = 1/(NoDivAA_properties$NodivLUCAclanratio_se[which(NoDivAA_properties$aaRS_Class == 2 | NoDivAA_properties$aaRS_Class == "1_2")]^2))
NoDivNSA_model_confintervalclass1 <- broom::augment(NoDivNSA_modelclass1, interval="confidence")
NoDivNSA_model_confintervalclass2 <- broom::augment(NoDivNSA_modelclass2, interval="confidence")
summary(NoDivNSA_model)

NoDivLUCA_NSA_plot <- ggplot(NoDivAA_properties, aes(y = as.numeric(NoDiv_LUCAusage), 
  x = as.numeric(Nb_nonH_sidechain), label = Letter)) + 
  xlab('Number of non-H in side chain') + ylab('LUCA clan usage')  + 
  theme(legend.position="none") + scale_x_continuous(breaks = c(0,2,4,6,8,10)) + 
  geom_text(color='blue', size = 12,  position = pd) +
  #annotate(geom="text", y=0.83, x=1, label="Class II", color="#619CFF", size = 12) +
  #annotate(geom="text", y=0.78, x=1, label="Class I", color="#F8766D", size = 12) +
  #annotate(geom="text", y=1.1, x=180, label="Class I & II", color="#00BA38", size = 12) +
  geom_errorbar(aes(ymin=as.numeric(NoDiv_LUCAusage)-NodivLUCAclanratio_se, ymax=as.numeric(NoDiv_LUCAusage)+NodivLUCAclanratio_se), position = pd) +
  theme(axis.text=element_text(size=24),axis.title=element_text(size=32,face="bold"),
        legend.text=element_text(size=14), legend.title = element_text(size=16,face="bold")) + 
  geom_line(aes(y = predict(NoDivNSA_model)), linewidth = 1,color = 'black') +
  annotate(geom="text", x=4, y=0.7, label=paste0(paste0('Weighted R',supsc('2')), "= 0.68"), color= "black", size = 14) +
  annotate(geom="text", x=4, y=0.64, label="p = 7e-6", color= "black", size = 14) +
  #annotate(geom="text", x=70, y=1.2, label="b)", color="black", size = 16) +
  geom_ribbon(aes(ymin=NoDivNSA_confinterval$.lower, ymax=NoDivNSA_confinterval$.upper), colour=NA, alpha=0.3)


BiosynReqmodified <-  NoDivAA_properties[- which(NoDivAA_properties$Letter == 'N' | NoDivAA_properties$Letter == 'Q' | NoDivAA_properties$Letter == 'L' ),]
BiosynReqmodified$Nb_nonH_sidechain[2] <- 6
BiosynReqmodified$aaRS_Class[2] <- 2
BiosynReqmodified$Letter[2] <- 'J'
BiosynReqmodified$aaRS_Class[9] <- 1
ModifiedNoDivNSA_model <- lm(BiosynReqmodified$NoDiv_LUCAusage ~ BiosynReqmodified$Nb_nonH_sidechain,
                           weights = 1/(BiosynReqmodified$NodivLUCAclanratio_se^2))
ModifiedNoDivNSA_modelclass2 <- lm(BiosynReqmodified$NoDiv_LUCAusage[which(BiosynReqmodified$aaRS_Class == 2)] ~ 
             BiosynReqmodified$Nb_nonH_sidechain[which(BiosynReqmodified$aaRS_Class == 2)],
            weights = 1/(BiosynReqmodified$NodivLUCAclanratio_se[which(BiosynReqmodified$aaRS_Class == 2  )]^2))
ModifiedNoDivNSA_modelclass1 <- lm(BiosynReqmodified$NoDiv_LUCAusage[which(BiosynReqmodified$aaRS_Class == 1 )] ~ 
             BiosynReqmodified$Nb_nonH_sidechain[which(BiosynReqmodified$aaRS_Class == 1  )],
           weights = 1/(BiosynReqmodified$NodivLUCAclanratio_se[which(BiosynReqmodified$aaRS_Class == 1 )]^2))
ModifiedNoDivNSA_model_confintervalclass1 <- broom::augment(ModifiedNoDivNSA_modelclass1, interval="confidence")
ModifiedNoDivNSA_model_confintervalclass2 <- broom::augment(ModifiedNoDivNSA_modelclass2, interval="confidence")
summary(ModifiedNoDivNSA_modelclass1 )

pd <- position_dodge(width = 0.5)
ModifiedNoDivLUCA_NSA_plot <- ggplot(BiosynReqmodified, aes(y = as.numeric(NoDiv_LUCAusage), 
  x = as.numeric(Nb_nonH_sidechain), label = Letter)) + 
  xlab('Number of non-H in side chain') + ylab('LUCA clan usage')  + 
  theme(legend.position="none") + geom_text(aes(colour = factor(aaRS_Class)), size = 16, position = pd) + 
  labs(color='aaRS Class') + scale_x_continuous(breaks = c(0,2,4,6,8,10)) +
  annotate(geom="text", y=0.8, x=1, label="Class II", color="#00BFC4", size = 12) +
  annotate(geom="text", y=0.75, x=1, label="Class I", color="#F8766D", size = 12) +
  geom_errorbar(aes(ymin=as.numeric(NoDiv_LUCAusage)-NodivLUCAclanratio_se, ymax=as.numeric(NoDiv_LUCAusage)+NodivLUCAclanratio_se), position = pd) +
  theme(axis.text=element_text(size=24),axis.title=element_text(size=32,face="bold"),
        legend.text=element_text(size=14), legend.title = element_text(size=16,face="bold")) + 
  #geom_line(aes(y = predict(NoDivNSA_model)), linewidth = 1,color = 'black') +
  #annotate(geom="text", x=4, y=0.73, label=paste0(paste0('Weighted R',supsc('2')), "= 0.69"), color= "black", size = 14) +
  #annotate(geom="text", x=4, y=0.67, label="p = 5e-6", color= "black", size = 14) +
  #annotate(geom="text", x=70, y=1.2, label="b)", color="black", size = 16) +
  #geom_ribbon(aes(ymin=NoDivNSA_confinterval$.lower, ymax=NoDivNSA_confinterval$.upper), colour=NA, alpha=0.3)
  geom_line(data= BiosynReqmodified[which(BiosynReqmodified$aaRS_Class==1  ),] ,aes(y = predict(ModifiedNoDivNSA_modelclass1)), linewidth = 1,color = "#F8766D") +
  geom_line(data= BiosynReqmodified[which(BiosynReqmodified$aaRS_Class==2),], aes(y = predict(ModifiedNoDivNSA_modelclass2)), linewidth = 1,color = "#00BFC4") +
  annotate(geom="text", x=3, y=0.65, label=paste0(paste0('Weighted R',supsc('2')), "= 0.52"), color= "#00BFC4", size = 14) +
  annotate(geom="text", x=3, y=0.6, label="p = 0.02", color= "#00BFC4", size = 14) +
  annotate(geom="text", x=7, y=1.2, label=paste0(paste0('Weighted R',supsc('2')), "= 0.95"), color= "#F8766D", size = 14) +
  annotate(geom="text", x=7, y=1.15, label="p = 4e-05", color= "#F8766D", size = 14) +
  geom_ribbon(data=BiosynReqmodified[which(BiosynReqmodified$aaRS_Class==1  ),] , aes(ymin=ModifiedNoDivNSA_model_confintervalclass1$.lower, ymax=ModifiedNoDivNSA_model_confintervalclass1$.upper), colour="#F8766D", alpha=0.3) +
  geom_ribbon(data=BiosynReqmodified[which(BiosynReqmodified$aaRS_Class==2 ),] , aes(ymin=ModifiedNoDivNSA_model_confintervalclass2$.lower, ymax=ModifiedNoDivNSA_model_confintervalclass2$.upper), colour="#00BFC4", alpha=0.3)



