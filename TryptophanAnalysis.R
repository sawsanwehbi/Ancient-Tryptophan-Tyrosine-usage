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
library(aphid)

setwd("~/Desktop/Tryptophan paper")
Tol_tree <- read.tree('Fixed_ToL.newick')
Distance2aaRS <- read.csv('AARSprotoenzyme_LUCA.csv', header = T)
PF00579_function <- read.tree('../Tryptophan paper/PF00579_function_geneTree.newick')
PF00579_aln <- read.fasta('../Tryptophan paper/PF00579_aln400.fasta', as.string = TRUE)
PF00579_conbac <- read.tree('../Tryptophan paper/PF00579_conHGTdropped_bac_branches.treefile')
plot.phylo(PF00579_conbac)
longbranchclade <- extract.clade(PF00579_conbac, interactive = T)
PF00579_conbac_longbranchdropped <- drop.tip(PF00579_conbac , longbranchclade$tip.label)
PF00579_conbac_longbranchdropped <- drop.tip(PF00579_conbac_longbranchdropped, 'G000403645_PF00579.Euryarchaeota.Archaea') 

write.tree(PF00579_conbac_longbranchdropped , 'PF00579_conbac_longbranchdropped_noLongArchaea.newick')
write.tree(unroot(PF00579_conbac_longbranchdropped) , 'PF00579_conbac_longbranchdropped_noLongArchaea_unrooted.newick')

PF00579_conbac_longbranchdropped <- multi2di(PF00579_conbac_longbranchdropped )
mad_PF00579_conbac_longbranchdropped <- mad( PF00579_conbac_longbranchdropped, 'full' )
plot.phylo(mad_PF00579_conbac_longbranchdropped[[6]][[1]])

write.fasta(PF00579_aln[which(names(PF00579_aln) %in% PF00579_conbac_longbranchdropped$tip.label) ],
            names(PF00579_aln)[which(names(PF00579_aln) %in% PF00579_conbac_longbranchdropped$tip.label)]
            , '../Tryptophan paper/PF00579_aln317.fasta')


PF00579_cprdropped <- drop.tip(PF00579_conbac , PF00579_conbac$tip.label[grepl('CPR',PF00579_conbac$tip.label)])
PF00579_cprdropped <- drop.tip(PF00579_cprdropped, 'G000403645_PF00579.Euryarchaeota.Archaea') 
write.tree(unroot(PF00579_cprdropped) , '../Tryptophan paper/PF00579_conarc_CPRdropped.newick')
write.fasta(PF00579_aln[which( names(PF00579_aln) %in% PF00579_cprdropped$tip.label )],
            names(PF00579_aln)[which( names(PF00579_aln) %in%PF00579_cprdropped$tip.label)]
            , '../Tryptophan paper/PF00579_alnCPRdropped.fasta')

plot.phylo(PF00579_cpr_BacterialYiso1dropped)
BacterialYiso1 <- extract.clade(PF00579_cprdropped, interactive = T)
PF00579_cpr_BacterialYiso1dropped <- drop.tip(PF00579_cprdropped, BacterialYiso1$tip.label)
write.tree(unroot(PF00579_cpr_BacterialYiso1dropped) , '../Tryptophan paper/PF00579_conarc_CPR_BacY1dropped.newick')
write.fasta(PF00579_aln[which( names(PF00579_aln) %in% PF00579_cpr_BacterialYiso1dropped$tip.label )],
            names(PF00579_aln)[which( names(PF00579_aln) %in% PF00579_cpr_BacterialYiso1dropped$tip.label)]
            , '../Tryptophan paper/PF00579_alnCPR_BacY1dropped.fasta')



reconciled.treefiles <- read.tree('YRS_nonrecombination_reconciledarcbac.newick')
write.tree(unroot(reconciled.treefiles), '../Tryptophan paper/YRS_nonrecombination_reconciledarcbac_unrooted.newick')
write.tree(dropped_treefiles, '../Tryptophan paper/PF00579_arcbac_reconciled_HGTdropped.newick')
write.tree(midpoint(reconciled.treefiles), '../Tryptophan paper/PF00579_arcbac_reconciled_midpoint.newick')


reconciled.treefiles_unrooted <- read.tree('PF00579_reconciledunrooted.newick')
unique_topologies <- unique.multiPhylo(reconciled.treefiles_unrooted[c(1:400)])
count <- function(item, list) {
  total = 0
  for (i in 1:length(list)) {
    if (all.equal.phylo(item, list[[i]], use.edge.length = FALSE)) {
      total = total + 1
    }
  }
  return(total)
}
result <- data.frame(unique_topology = rep(0, length(unique_topologies)),
                     count = rep(0, length(unique_topologies)))
for (i in 1:length(unique_topologies)) {
  result[i, ] <- c(i, count(unique_topologies[[i]], reconciled.treefiles_unrooted ))
}
result$percentage <- ((result$count/length(reconciled.treefiles_unrooted ))*100)


classIandII <- read.fasta('../sequence.fasta', as.string = T)
classIandIIpfam <- as.list(c(substr(classIandII$`ClassI_TyrRS_NC_000913.3:c1717222-1715948`, 685,1216),
substr(classIandII$`ClassII_SerRS_NC_000913.3:939428-940720`, 91,982)))
classIandIIpfam[[1]] <- reverse(classIandIIpfam[[1]])
write.fasta(classIandIIpfam, names(classIandII), '../ClassIandII_pfamNuc.fasta')


classIandIIpfam2ndpos <-as.list(c( paste0( str_split_1(classIandIIpfam[[1]],"")[seq(2,532,2)], collapse = ""),
       paste0( str_split_1(classIandIIpfam[[2]],"")[seq(2,892,2)], collapse = "") ))
write.fasta(classIandIIpfam2ndpos, names(classIandII), '../ClassIandII_pfamNuc2ndpos.fasta')

##### structural analysis ####
pdbfam <- read.delim('../Differential-Retention-main/PDBfam.txt', header = TRUE)
PF00579_pdbfam <- pdbfam[pdbfam$Pfam_Acc == 'PF00579',]
writeLines(unique(PF00579_pdbfam$UniprotCode), 'PF00579_uniprot')
pf00579_namesdf <- read.delim('../Tryptophan paper/idmapping_2024_11_20.tsv')
pf00579_rows <- read.csv('../Tryptophan paper/PF00579_rows.csv')
pf00579_rows <- pf00579_rows[-which(duplicated(pf00579_rows$From)),]
pf00579_rows <- inner_join(pf00579_rows ,pf00579_namesdf, by ='From')
pf00579_matching <- read.delim('../Tryptophan paper/PF00579matches_Uniprot.tsv', header = FALSE)
colnames(pf00579_matching)[1] <- 'Entry'
pf00579_rows <- inner_join( pf00579_rows, pf00579_matching, by = 'Entry')
Protein3Di_seq <- read.fasta('../Tryptophan paper/PF00579protein147_3Dialn.fa')
pf00579_rows$protein3Di_sequence <- Protein3Di_seq[match(pf00579_rows$Entry , names(Protein3Di_seq))]
pf00579_rows <- pf00579_rows[-which(nchar(pf00579_rows$protein3Di_sequence) == 4),] # remove empty sequences

# get EMBL CDS or RefSeq protein accesions from inteproscan df, map to uniprot using ID mapping
# use uniprot id to download pdb
writeLines(sapply(1:nrow(OnePfamUniprotPDBnononcanonYRS), function(i){
  paste0('curl https://alphafold.ebi.ac.uk/files/AF-',OnePfamUniprotPDBnononcanonYRS$Entry[i],
         '-F1-model_v4.pdb  -o AF-',OnePfamUniprotPDBnononcanonYRS$Entry[i],'-F1.pdb')  }))

tree_3di <- read.tree('../Tryptophan paper/PF00579_3DiAF_branches.treefile')
aln_3di <- seqinr::read.fasta('../Tryptophan paper/PF00579_Trunc_3Dialn.fasta', as.string=T)
pf00579_seqdf <- read.csv('PF00579_rows.csv', header =T)
UniprotPDBidmapping <- read.csv('UniprotPDBidmapping.csv', header =T)
UniprotPDBidmapping <- inner_join(pf00579_seqdf,UniprotPDBidmapping , 'From')
UniprotPDBidmapping <- UniprotPDBidmapping[-which(duplicated(UniprotPDBidmapping$Entry)),]
UniprotPDBidmapping <- UniprotPDBidmapping[-which(duplicated(UniprotPDBidmapping$From)),]
UniprotPDBidmapping$Protein.names <- gsub('-.*','',UniprotPDBidmapping$Protein.names)
UniprotPDBidmapping$SequenceFuncNames <- paste0(UniprotPDBidmapping$accesion_names,'_', UniprotPDBidmapping$Entry,".", 
               UniprotPDBidmapping$Protein.names, '.', UniprotPDBidmapping$SequenceNames,".",  UniprotPDBidmapping$Domain)
UniprotPDBidmapping$SequenceFuncNames <- gsub('..........._PF00579', '', UniprotPDBidmapping$SequenceFuncNames)
UniprotPDBnononcanonYRS <- UniprotPDBidmapping[-which(UniprotPDBidmapping$accesion_names %in%
                              gsub('_PF00579.*','',names(BacYRSnoncanon_aln)))  ,]
UniprotPDBnononcanonYRS <- UniprotPDBnononcanonYRS[-which(grepl('CPR',UniprotPDBnononcanonYRS$SequenceFuncNames)),]

## Using I3IM11  as a ref (length 333AA) PFAM MATCH 4:287 / 124:661 IN 210 protein 3di alignment 
# 123: 653 in 209 3di aln with no 
Protein3Di_seq <- seqinr::read.fasta('../foldmason_ss.fa',as.string = T)
Truncated_alignedsequences <- substr(as.list(Protein3Di_seq),123,653)
write.fasta(as.list(Truncated_alignedsequences) , 
            UniprotPDBnononcanonYRS$SequenceFuncNames[match(names(Protein3Di_seq ),UniprotPDBnononcanonYRS$Entry)], 
            'PF00579protein210_3DialnNoG000403645.fa')
# penalty gap opening 15, gap extension 1

reconciledcon3di_tree <- read.tree('../PF00579protein210_3Dialn.fa.treefile')
reconciledcon3di_tree <- drop.tip(reconciledcon3di_tree , which(grepl('G000403645', reconciledcon3di_tree$tip.label)))
write.tree(reconciledcon3di_tree ,'../PF00579protein210_3DialnNoG000403645.treefile' )


Protein3Diaa_seq <- seqinr::read.fasta('../foldmason_aa.fa',as.string = T)
write.fasta(as.list(Protein3Diaa_seq) , 
            UniprotPDBnononcanonYRS$SequenceFuncNames[match(names(Protein3Diaa_seq ),UniprotPDBnononcanonYRS$Entry)], 
            '../Tryptophan paper/PF00579foldmason_aa.fa')

reconciled3Di_trees <- read.tree('3DIPF00579.rec_uml')
write.tree(unroot(reconciled3Di_trees), '3DIPF00579_reconciled.newick')
partition_trees <- read.tree('PF00579_partition.rec_uml')
write.tree(unroot.multiPhylo(partition_trees), 'PF00579_partition.rec_unrooted.newick')

PF00579aln_3di <- seqinr::read.fasta('PF00579protein210_3DialnNoG000403645.fa', as.string=F)
PF00579aln_3diAA <- seqinr::read.fasta('PF00579protein210_3DiAAalnNoG000403645.fa', as.string=F)
PF00579aln_3diAA <- PF00579aln_3diAA[match(names(PF00579aln_3di), names(PF00579aln_3diAA))]

nchar(PF00579aln_3di$G000363885_N6VRA9.Tryptophan.Euryarchaeota.Archaea[1])
# 1:755 AA aln, 756:1286 3Di aln
write.fasta(as.list(paste0(PF00579aln_3diAA,PF00579aln_3di)), names(PF00579aln_3diAA),
            'PF00579_AA3Dialn.fasta')

## replace aa sites in muscle alignment with 3Di
Musclealn_3Di <- list()
for ( seq in 1:209){
template_chars <-  PF00579aln_3diAA[[seq]]
fill_chars <- PF00579aln_3di[[seq]][-which(PF00579aln_3di[[seq]] == '-')]

fill_index <- 1
for (i in seq_along(template_chars)) {
  if (template_chars[i] != "-") {
    template_chars[i] <- fill_chars[fill_index]
    fill_index <- fill_index + 1
  }
}

Musclealn_3Di <- append (Musclealn_3Di, paste(template_chars, collapse = "")) }
PF00579aln_3diAAseqs <- sapply(1:209,function(x){paste(PF00579aln_3diAA[[x]], collapse = "")})
write.fasta(as.list(paste0(PF00579aln_3diAAseqs,Musclealn_3Di)), names(PF00579aln_3diAA),
            'PF00579_AA3Di_Musclealn.fasta')
write.fasta(as.list(Musclealn_3Di), names(PF00579aln_3diAA),
            'PF00579_3Di_Musclealn.fasta')

accesion_line <- gsub(',', "",toString(accesion_Df[,2]))
paste('./dataformat download gene accesion', accesion_line )
      
parttrees <- read.tree('../PF00579_Musclepartition_con.parttrees')
pf00579_AAtree <- parttrees[[1]]
pf00579_3ditree <- parttrees[[2]]
trp_3ditree <- keep.tip(parttrees[[2]], grep('Tryptophan', parttrees[[2]]$tip.label))
trp_AAtree <- keep.tip(parttrees[[1]], grep('Tryptophan', parttrees[[1]]$tip.label))


# Assuming you have a list of trees named 'trees'
library(ggtree)
cls <-list(bacteria=grep('Bacteria',pf00579_AAtree$tip.label ),
          Archaea=grep('Archaea',pf00579_AAtree$tip.label ))
pf00579_AAtree_new <- groupOTU(pf00579_AAtree,.node=cls)
pf00579_3ditree_new <- groupOTU(pf00579_3ditree,.node=cls)

  
ggtree(c(pf00579_AAtree_new ,pf00579_3ditree_new ), root.position = 3) + 
  facet_wrap(~.id, scales = "fixed") +
  geom_treescale(width = 3.0, x=0) +
  geom_tree(aes(color=group)) +
  scale_color_manual(values=c("blue","red")) +
  #geom_tiplab(size=0.8,aes(color=group), show.legend=FALSE)+
  theme(legend.position="none") 


generaxtree$tip.label[grep('G000091085', generaxtree$tip.label)] #Chlamydophila pneumoniae strain AR39. Y
generaxtree$tip.label[grep('G000011125', generaxtree$tip.label)] #Aeropyrum pernix strain K1 W
generaxtree <- read.tree('../Tryptophan paper/PF00579_Musclepartition_con3Di.newick')
is.binary.phylo(binarygeneraxtree)
binarygeneraxtree <- multi2di.phylo(generaxtree)
write.tree(binarygeneraxtree, '../Tryptophan paper/PF00579_Musclepartition_con3Di.newick')

NonCanonYRS_recomb <- read.fasta('NonCanon_BacYRS_recombination.fasta',as.string = T)
median(nchar(NonCanonYRS_recomb))


#### sample AA partition ####
for (i in 3:100){
subsample <- sample(1:755, 75, FALSE) #  453, 302 sites are 60% and 40% of aln
subsampleAA <- lapply(1:209,function(x){PF00579aln_3diAA[[x]][subsample]})
subsampleAAesqs <-  sapply(1:209,function(x){paste(subsampleAA[[x]], collapse = "")})
write.fasta(as.list(paste0(subsampleAAesqs,Musclealn_3Di)), names(PF00579aln_3diAA),
           paste0('PF00579_10subsampleAA3Di_Musclealn_',i, '.fasta'))}

## sample 3Di partition
Musclealn_3Dionly <-  seqinr::read.fasta( 'PF00579_3Di_Musclealn.fasta', as.string=F)
PF00579aln_AAonly <- seqinr::read.fasta('PF00579protein210_3DiAAalnNoG000403645.fa', as.string=T)
PF00579aln_AAonly <- PF00579aln_AAonly[match(names(PF00579aln_3di), names(PF00579aln_AAonly))]

for (i in 3:100){
  subsample <- sample(1:755, 75, FALSE) # 453, 302 sites are 60% and 40% of aln
 # print(subsample)
  subsample3Di <- lapply(1:209,function(x){Musclealn_3Dionly[[x]][subsample]})
  subsample3Diesqs <-  sapply(1:209,function(x){paste(subsample3Di[[x]], collapse = "")})
  write.fasta(as.list(paste0(PF00579aln_AAonly, subsample3Diesqs)), names(PF00579aln_3diAA),
              paste0('PF00579_10AAsubsample3Di_Musclealn_',i, '.fasta'))}

nchar(as.list(paste0(PF00579aln_AAonly, subsample3Diesqs))[[1]])
writeLines(as.character(1:100), 'Bootstrapfiles')

#### Ancestral state reconstruction ####
#Trp LUCA node is 59 (0/2 sites are >50% W), HGT node is Node 2 (0/3 sites are >50% W), and bacterial tyr node 174 (4/6 sites 99% W)
state_file <- read.delim('PF00579_partNQ_asr.state', comment.char = "#")
state_file <- read.delim('GapPF00579protein210_3DiAAalnNoG000403645.fa.state', comment.char = "#")
root_state <- state_file[which(state_file$Node == 'Node1'),]
root_state[which(root_state$State == 'Y'),]
state_file[which(state_file$Node == 'Node59' & state_file$State == 'W'),]
state_file[which(state_file$Node == 'Node2' & state_file$State == 'Y'),]
state_file[which(state_file$Node == 'Node174' & state_file$State == 'Y'),]
state_file[which(state_file$Node == 'Node59' & state_file$State == 'W'),]
# AT THE ROOT THERE ARE 5 sites that are most likely W but with very low probabilities <50%
# 37,38,33,39,40. The second most like aa at these site are hydrophobic F, L, L, F,Y respectively
sort(root_state[which(root_state$State == 'W'),][1,4:23], decreasing=T)
state_file[which(state_file$Node == 'Node2' & state_file$Site == 134),]

anc_Seq <- state_file$State[which(state_file$Node == 'Node174')]
names(anc_Seq ) <- 1:length(anc_Seq )
anc_Seq[-which(anc_Seq == '-')][which(names(anc_Seq)[-which(anc_Seq == '-')] == '134')]
length(anc_Seq[-which(anc_Seq == '-')])

state_file[which(state_file$Node == 'Node174' & state_file$State == 'Y'),][
  which(state_file[which(state_file$Node == 'Node174' & state_file$State == 'Y'),]$p_Y > 0.7),]$Site



top4_probdf <- data.frame()
for (i in 1:length(which(state_file$Node == 'Node174'  & state_file$State == 'W'))){
ancsite_prob <- state_file[which(state_file$Node == 'Node174'  & state_file$State == 'W'),][i,4:23]
top4_prob <- ancsite_prob[which(ancsite_prob %in%  sort(as.numeric(ancsite_prob ), decreasing=T)[1:4] )]
top4_probvector <- as.numeric(top4_prob)
names(top4_probvector) <- gsub('p_','', colnames(top4_prob))
top4_probvector <- sort(top4_probvector, decreasing = T)
top4_probvector <- paste0(names(top4_probvector), " (" ,top4_probvector , ")")
top4_probdf <-  rbind(top4_probdf,top4_probvector )}
colnames(top4_probdf) <- c('1st most likely state','2nd most likely state','3rd most likely state','4th most likely state')

Node <- rep('Most recent common ancestor',length(which(state_file$Node == 'Node2'  & state_file$State == 'Y')))
Site <- state_file[which(state_file$Node == 'Node2'  & state_file$State == 'Y'),2]
top4_probdf <- cbind(Node,Site,top4_probdf )

Roottop4_probdf <- top4_probdf 
Y_top4_probdf <- rbind(Roottop4_probdf, top4_probdf )
plot(tableGrob(Y_top4_probdf, theme = ttheme_default(base_size = 20)))
write.csv(Y_top4_probdf, 'Tyrlikelysites.csv') #12X20

#Ancestral W sites in G000091085 YRS: W90, -,-, W205,L220, -,-,-
#Ancestral W sites in G000011125 WRS: R84, -,-,F197,L213,-,-,-
Site_G000091085YRS <- c('W90', '-','-', 'W205','L220', '-','-','-') #Chlamydophila pneumoniae strain AR39. Y
Site_G000011125WRS <- c('R84', '-','-','F197','L213','-','-','-') #Aeropyrum pernix strain K1 W

W_top4_probdfwithRefs <- cbind(W_top4_probdf,Site_G000091085YRS, Site_G000011125WRS)
W_top4_probdfwithRefs$Node[which(W_top4_probdfwithRefs$Node=="Most recent common ancestor")] <-
  'MRCA'

colnames(W_top4_probdfwithRefs )[7:8] <- c('In C. pneumoniae YRS', 'In A. pernix WRS')
plot(tableGrob(W_top4_probdfwithRefs , theme = ttheme_default(base_size = 20)))

#### YRS recombination analysis#####
all_aln_fasta <- seqinr::read.fasta('../Tryptophan paper/PF00579_alnCPRdropped.fasta', as.string = T)
aln_fasta <- seqinr::read.fasta('../Tryptophan paper/YRS_aln.fasta', as.string = T)
names(aln_fasta) <- gsub('/.*','', names(aln_fasta))
CPRdropped_arctree <- read.tree('../Tryptophan paper/PF00579_arc_CPRdropped.treefile')
plot.phylo(CPRdropped_arctree, show.tip.label = F)
BacYRScanon <- extract.clade(CPRdropped_arctree, interactive = T)
BacYRSnoncanon <- extract.clade(CPRdropped_arctree, interactive = T)
ArcYRS <- extract.clade(CPRdropped_arctree, interactive = T)
ArcWRS <- extract.clade(CPRdropped_arctree, interactive = T)
BacWRS <- extract.clade(CPRdropped_arctree, interactive = T)
AllBacWRS <- extract.clade(CPRdropped_arctree, interactive = T)
# HPLDLK region in weblogo 610-625 ... 413-430
# 55-75 INSERTIONS AT N TERMINAL 32-52
#write.tree(ArcYRS, 'ArcYRS.treefile')
PF00579_enforce <- read.tree('../PF00579_enforceroot.treefile')
PF00579_enforce <- drop.tip(PF00579_enforce, 
                   c(PF00579_enforce$tip.label[which(grepl('CPR',PF00579_enforce$tip.label))],
                    "G000403645_PF00579.Euryarchaeota.Archaea"))

                           
seqinr::write.fasta(as.list(all_aln_fasta[which(names(all_aln_fasta ) %in% ArcWRS$tip.label)] ), 
                    names(all_aln_fasta)[which(names(all_aln_fasta ) %in% ArcWRS$tip.label)],
                    '../Tryptophan paper/AllPF00579_ArcWRS.fasta')

seqinr::write.fasta(as.list(all_aln_fasta[which(names(all_aln_fasta ) %in% ArcYRS$tip.label)] ), 
                    names(all_aln_fasta)[which(names(all_aln_fasta ) %in% ArcYRS$tip.label)],
                    '../Tryptophan paper/AllPF00579_ArcYRS.fasta')

seqinr::write.fasta(as.list(all_aln_fasta[which(names(all_aln_fasta ) %in% BacWRS$tip.label)] ), 
                    names(all_aln_fasta)[which(names(all_aln_fasta ) %in% BacWRS$tip.label)],
                    '../Tryptophan paper/AllPF00579_BacWRS.fasta')

seqinr::write.fasta(as.list(all_aln_fasta[which(names(all_aln_fasta ) %in% BacYRSnoncanon$tip.label)] ), 
                    names(all_aln_fasta)[which(names(all_aln_fasta ) %in% BacYRSnoncanon$tip.label)],
                    '../Tryptophan paper/AllPF00579_BacYRSnoncanon.fasta')

seqinr::write.fasta(as.list(all_aln_fasta[which(names(all_aln_fasta ) %in% BacYRScanon$tip.label)] ), 
                    names(all_aln_fasta)[which(names(all_aln_fasta ) %in% BacYRScanon$tip.label)],
                    '../Tryptophan paper/AllPF00579_BacYRScanon.fasta')

# CANON 73, NON CANON 44, ARCHAEA 61
BacYRSnoncanon_aln <- seqinr::read.fasta('PF00579_BacYRSnoncanon.fasta', as.string = T)
BacYRScanon_aln <- seqinr::read.fasta('PF00579_BacYRScanon.fasta', as.string = T)
ArcYRS_aln <- seqinr::read.fasta('PF00579_ArcYRS.fasta', as.string = T)
names(BacYRScanon_aln) <- paste0(names(BacYRScanon_aln),'.canon')
names(BacYRSnoncanon_aln) <- paste0(names(BacYRSnoncanon_aln),'.noncanon')
seqinr::write.fasta(as.list(c(str_sub(BacYRSnoncanon_aln , 214, 419), str_sub(BacYRScanon_aln , 214, 419), str_sub(ArcYRS_aln , 214, 419))),
                    c(names(BacYRSnoncanon_aln), names(BacYRScanon_aln), names(ArcYRS_aln)),
                    'PF00579_YRSrecombination.fasta')

seqinr::write.fasta(as.list(c(str_sub(BacYRSnoncanon_aln , 214, 415), str_sub(BacYRScanon_aln , 214, 415), str_sub(ArcYRS_aln , 214, 415))),
                    c(names(BacYRSnoncanon_aln), names(BacYRScanon_aln), names(ArcYRS_aln)),
                    '../Tryptophan paper/PF00579_YRSrecombination214_415.fasta')


seqinr::write.fasta(as.list(c(lapply(str_sub_all(BacYRSnoncanon_aln , c(1,420), c(213,448)), toString),
                              lapply(str_sub_all(BacYRScanon_aln , c(1,420), c(213,448)), toString),
                              lapply(str_sub_all(ArcYRS_aln , c(1,420), c(213,448)), toString))),
                    c(names(BacYRSnoncanon_aln), names(BacYRScanon_aln), names(ArcYRS_aln)),
                    '../Tryptophan paper/PF00579_YRSnonrecombination.fasta')


AllPF00579_noCPR_labelled <- seqinr::read.fasta('../Tryptophan paper/AllPF00579_noCPR_labelled.fasta',as.string = T)
seqinr::write.fasta(as.list(str_sub(AllPF00579_noCPR_labelled , 286, 620)),
                    names(AllPF00579_noCPR_labelled),
                    '../Tryptophan paper/AllPF00579_noCPR_labelled_recombination.fasta')

aa_sequences <- readFASTA('../Tryptophan paper/NonCanon_BacYRS_recombination.fasta')
aa_sequences_abin <- toupper(aa_sequences )
names(aa_sequences_abin) <- names(aa_sequences)
RecombinationHMM <- derivePHMM(aa_sequences_abin , residues='AMINO')
writePHMM(RecombinationHMM , 'YRSrecombination.hmm')


sample1 <- sample(44, 10)
sample2 <- sample(73, 10)
sample3 <- sample(61, 10)
seqinr::write.fasta(as.list(c(BacYRSnoncanon_aln[sample1],BacYRScanon_aln[sample2],
     ArcYRS_aln[sample3] )), c(names(BacYRSnoncanon_aln)[sample1], 
     names(BacYRScanon_aln)[sample2],names(ArcYRS_aln)[sample3] ), 'YRS10sample2_GARD.fasta')

YRS_PD <- read.csv('YRS_PD.csv', header = F)

seqinr::write.fasta(as.list(aln_fasta[which(names(aln_fasta) %in% YRS_PD$V2)]), 
         names(aln_fasta)[which(names(aln_fasta) %in% YRS_PD$V2)], 'YRS_PDGARD15')

seqinr::write.fasta(as.list(c(lapply(str_sub_all(BacYRSnoncanon_aln , c(214,323,384), c(311,361,419)), toString),
                              lapply(str_sub_all(BacYRScanon_aln , c(214,323,384), c(311,361,419)), toString),
                              lapply(str_sub_all(ArcYRS_aln ,c(214,323,384), c(311,361,419)), toString))),
                    c(names(BacYRSnoncanon_aln), names(BacYRScanon_aln), names(ArcYRS_aln)),
                    '../Tryptophan paper/PF00579_YRSnoinsertion_recombination.fasta')

seqinr::write.fasta(as.list(all_aln_fasta[-which(names(all_aln_fasta ) %in% names(BacYRSnoncanon_aln))]),
                        names(all_aln_fasta)[-which(names(all_aln_fasta ) %in% names(BacYRSnoncanon_aln))] ,
                    '../Tryptophan paper/PF00579_NonCanonCPRdropped.fasta')

yrs_recon <- seqinr::read.fasta('../Tryptophan paper/PF00579_YRSrecombination214_415.fasta', as.string = T)
write.fasta(as.list(yrs_recon[-which(nchar(yrs_recon) < 80)]), 
            names(yrs_recon)[-which(nchar(yrs_recon) < 80)], '../Tryptophan paper/PF00579_YRSrecombination214_415.fasta')


write.fasta(as.list(c(BacYRSnoncanon_aln ,BacYRScanon_aln )), c(names(BacYRSnoncanon_aln ), names(BacYRScanon_aln )),
            'AllBacYRS_aln.fasta')
write.fasta(as.list(c(ArcYRS_aln  ,BacYRScanon_aln )), c(names(ArcYRS_aln  ), names(BacYRScanon_aln )),
            'ArcBacYRScanon_aln.fasta')
write.fasta(as.list(c(ArcYRS_aln  ,BacYRSnoncanon_aln )), c(names(ArcYRS_aln ), names(BacYRSnoncanon_aln )),
            'ArcBacYRSnoncanon_aln.fasta')

AllBacYRS_aln <- read.alignment('AllBacYRS_aln.fasta','fasta')
ArcBacYRScanon_aln <- read.alignment('ArcBacYRScanon_aln.fasta','fasta') 
ArcBacYRSnoncanon_aln <- read.alignment('ArcBacYRSnoncanon_aln.fasta','fasta') 
AllBacdist <- dist.alignment(AllBacYRS_aln ,"similarity")
AllBacdist <- as.matrix(AllBacdist)
AllBacdist <- AllBacdist[which(grepl('noncanon',AllBacYRS_aln$nam)),which(grepl('[.]canon',AllBacYRS_aln$nam)) ]

ArcBacYRScanondist <- dist.alignment(ArcBacYRScanon_aln ,"similarity")
ArcBacYRScanondist <- as.matrix(ArcBacYRScanondist)
ArcBacYRScanondist <- ArcBacYRScanondist[which(grepl('Archaea',ArcBacYRScanon_aln$nam)),which(grepl('[.]canon',ArcBacYRScanon_aln$nam)) ]

ArcBacYRSnoncanondist <- dist.alignment(ArcBacYRSnoncanon_aln ,"similarity")
ArcBacYRSnoncanondist <- as.matrix(ArcBacYRSnoncanondist)
ArcBacYRSnoncanondist <- ArcBacYRSnoncanondist[which(grepl('Archaea',ArcBacYRSnoncanon_aln$nam)),which(grepl('noncanon',ArcBacYRSnoncanon_aln$nam)) ]

1- 0.44^2
1 - (mean(AllBacdist))^2 # 66% similarity 30%
1 - (mean(ArcBacYRScanondist))^2 #58 similarity  20%
1 - (mean(ArcBacYRSnoncanondist))^2 #63% similarity  24%
(mean(AllBacdist)^2)
### Color bacterial tree with non canon vs canon YRS ####
BacYRStips <- c(names(BacYRScanon_aln),names(BacYRSnoncanon_aln))
bacterialYRSspeciestree <- keep.tip(Tol_tree, Tol_tree$tip.label[c(
  which(Tol_tree$tip.label %in% gsub('_.*','',names(BacYRScanon_aln)) ),
  which(Tol_tree$tip.label %in% gsub('_.*','',names(BacYRSnoncanon_aln)) ))])

bacterialYRSspeciestips <- BacYRStips[match(bacterialYRSspeciestree$tip.label,gsub('_.*','',BacYRStips))]
bacterialYRSspeciestips <- gsub('.*_PF00579.','',bacterialYRSspeciestips)
bacterialYRSspeciestips <- gsub('.Bacteria','',bacterialYRSspeciestips)

colortaxa <- rep('black', nrow(bacterialYRSspeciestree$edge))
colortaxa[which(bacterialYRSspeciestree$edge[,2] %in% which(bacterialYRSspeciestree$tip.label %in%  gsub('_.*','',names(BacYRScanon_aln)) ))] <- 'red'
colortaxa[which(bacterialYRSspeciestree$edge[,2] %in% which(bacterialYRSspeciestree$tip.label %in%  gsub('_.*','',names(BacYRSnoncanon_aln)) ))]  <- 'magenta'
bacterialYRSspeciestree$tip.label <- bacterialYRSspeciestips
plot.phylo(bacterialYRSspeciestree , show.tip.label = F,type= 'phylogram', 
           edge.color = colortaxa, edge.width = 2 , align.tip.label = T) +
  add.scale.bar(x = 0.5, y = -1,cex =2, font = 15, lwd =2)


bacterialYRSspeciestree$tip.label[which(bacterialYRSspeciestree$tip.label %in%gsub('_.*','',names(BacYRScanon_aln)) )] <- 
  paste0(bacterialYRSspeciestree$tip.label[which(bacterialYRSspeciestree$tip.label %in%gsub('_.*','',names(BacYRScanon_aln)) )], 'YRS_canon')
bacterialYRSspeciestree$tip.label[which(bacterialYRSspeciestree$tip.label %in%gsub('_.*','',names(BacYRSnoncanon_aln)) )] <- 
  paste0(bacterialYRSspeciestree$tip.label[which(bacterialYRSspeciestree$tip.label %in%gsub('_.*','',names(BacYRSnoncanon_aln)) )], 'YRS_noncanon')
write.tree(bacterialYRSspeciestree, 'BacterialYRSspeciestree.newick')

library(diversitree)
state <- rep(1, length(bacterialYRSspeciestree$tip.label))
state[which(bacterialYRSspeciestree$tip.label %in%  gsub('_.*','',names(BacYRScanon_aln)))] <- 0
phy <- force.ultrametric(bacterialYRSspeciestree)
names(state) <- phy$tip.label
#f <- make.bisse(phy, state)
#lik <- make.bisse.td(phy, state, n.epoch = 2)
lik <- make.bisse(phy, state)
pars <- c(.3, .3, .02, .02, 0.02, 0.02)
fit <- find.mle(lik, pars, method="minqa") 
logLik(fit)
st <- asr.marginal(lik, coef(fit))
nodelabels(thermo=t(st), piecol=5:6, cex=.5)


##### Plot roots on unrooted tree #####
YRSrecomb_conbac <- read.tree('../YRS_recombination132_reconciledarcbaccons_bac.treefile')
#YRSrecomb_conbac <- read.tree('../PF00579_YRSnonrecombination_arcbacreconciled_bac.treefile')
mp_YRSrecomb_conbac <- phangorn::midpoint(YRSrecomb_conbac )
unrootedtree <- unroot.phylo(mp_YRSrecomb_conbac )
MADroot <- mad(unrootedtree, 'full' )
MADrootree <-MADroot[[6]][[1]]


sum(mp_YRSrecomb_conbac$edge.length[which(mp_YRSrecomb_conbac$edge[,1] == getRoot(mp_YRSrecomb_conbac))])
2.1939350/2.83947 # recombinant, longer branch is towards canon bacteria 77%
1.4721509/2.421396 # nonrecombinant, longer branch is towards archaea 60%

sum(MADrootree$edge.length[which(MADrootree$edge[,1] == getRoot(MADrootree))])
1.569889/2.83947 # recombinant, longer branch is towards canon bacteria 55%
1.341029/2.421396 # nonrecombinant, longer branch is towards archaea 55%

colortaxa <- rep('black', nrow(unrootedtree$edge))
colortaxa[which(unrootedtree$edge[,2] %in% which(grepl('Archaea',unrootedtree$tip.label)))] <- 'blue'
colortaxa[which(unrootedtree$edge[,2] %in% which(grepl('canon',unrootedtree$tip.label)))] <- 'red'
colortaxa[which(unrootedtree$edge[,2] %in% which(grepl('noncanon',unrootedtree$tip.label)))] <- 'magenta'
plot.phylo(unrootedtree, show.tip.label = F,type= 'unrooted', edge.color = colortaxa ) 
edgelabels(pch=19,edge=which(unrootedtree$edge.length > 2.83947), adj = c(0.53, 0.4), col='green') 
edgelabels(pch=19, edge=which(unrootedtree$edge.length > 2.83947), adj = c(0.68, -0.05))

    
colortaxa <- rep('black', nrow(unrootedtree$edge))
colortaxa[which(unrootedtree$edge[,2] %in% which(grepl('Archaea',unrootedtree$tip.label)))] <- 'blue'
colortaxa[which(unrootedtree$edge[,2] %in% which(grepl('canon',unrootedtree$tip.label)))] <- 'red'
colortaxa[which(unrootedtree$edge[,2] %in% which(grepl('noncanon',unrootedtree$tip.label)))] <- 'magenta'
plot.phylo(unrootedtree, show.tip.label = F,type= 'unrooted', edge.color = colortaxa ) 
#edgelabels(pch=19,edge=205, adj = c(0.5, 0.5), col='red') 
edgelabels(pch=19,edge=which(unrootedtree$edge.length > 2.421396), adj = c(0.55, 0.4) , col ='green') 
edgelabels(pch=19, edge=which(unrootedtree$edge.length > 2.421396), adj = c(0.61, 0.3))


##### BENNU ANALYSIS #####
AA_properties <- read.csv('../PFAM Trees/AminoAcid_properties.csv', header = T)
Bennu_aa <- read.csv('../PFAM Trees/Bennu_AA.csv', header =T)
Bennu_presence <- read.csv('../BennuAApresence_Angel.csv', header = T)
Ancient_aa <- read.csv('../PFAM Trees/AncientAAusage.csv', header =T)
cor.test(AA_properties$Relative_aafreq_CR2meteorites[-which(is.na(AA_properties$Relative_aafreq_CR2meteorites))],
         Ancient_aa$LUCA.usage[-which(is.na(AA_properties$Relative_aafreq_CR2meteorites))])


cor.test(Ancient_aa$LUCA.usage ,Bennu_aa$Experimental.m.z , method= 'spearman' )
cor.test(Ancient_aa$Pre.LUCA.usage ,Bennu_aa$Experimental.m.z, method= 'spearman' )

summary(lm(Ancient_aa$LUCA.usage ~ Bennu_aa$Experimental.m.z, weights = (1/(Ancient_aa$LUCA.usage.SE)^2)))
summary(lm(Ancient_aa$Pre.LUCA.usage ~ Bennu_aa$Experimental.m.z, weights = (1/(Ancient_aa$Pre.LUCA.usage.SE)^2)))


summary(lm(Ancient_aa$LUCA.usage ~ Bennu_presence$OREX080010700...aggregate., weights = (1/(Ancient_aa$LUCA.usage.SE)^2)))
summary(lm(Ancient_aa$Pre.LUCA.usage ~ Bennu_presence$OREX080010700...aggregate., weights = (1/(Ancient_aa$Pre.LUCA.usage.SE)^2)))

summary(lm(Ancient_aa$LUCA.usage ~ Bennu_presence$OREX08000550113..angular., weights = (1/(Ancient_aa$LUCA.usage.SE)^2)))
summary(lm(Ancient_aa$Pre.LUCA.usage ~ Bennu_presence$OREX08000550113..angular., weights = (1/(Ancient_aa$Pre.LUCA.usage.SE)^2)))
# preluca correlates with angular and aggregate, LUCA with angular, none with hummocky and mottled

summary(lm(Ancient_aa$LUCA.usage ~ Bennu_presence$OREX080010700...aggregate. + AA_properties$MolecularWeightsDa,
           weights = (1/(Ancient_aa$LUCA.usage.SE)^2)))

rm(bennuvsusage_df)
bennuvsusage_df <- cbind(Ancient_aa , Bennu_aa , AA_properties, Bennu_presence)
bennuvsusage_df <- bennuvsusage_df[,-which(duplicated(colnames(bennuvsusage_df)))]

ggplot(bennuvsusage_df, aes(y = as.numeric(LUCA.usage), 
                            x = as.numeric(Experimental.m.z), label = Letter)) + 
  xlab('Bennu abundance') +
  ylab('LUCA clan usage') + 
  geom_text(aes(colour = factor(Moosmann_category)), size = 16) + 
  labs(color='Moosmann category') + 
  theme(axis.text=element_text(size=24),axis.title=element_text(size=32,face="bold"),
        legend.text=element_text(size=27), legend.title = element_text(size=16,face="bold")) + 
  scale_fill_brewer(palette="Dark2") +  theme(legend.position = 'none') +
  guides(color = guide_legend(override.aes = list(size = 16))) + 
  annotate(geom="text", y=.75, x=200, label=paste0("Weighted R2 = 0.25"), color="black", size = 14) +
  annotate(geom="text", y=.7, x=200, label="p = 0.01", color="black", size = 14) 


AngularvsLUCA <- ggplot(bennuvsusage_df, aes(y = as.numeric(LUCA.usage), 
                            x = as.numeric(Bennu_presence$OREX08000550113..angular.), label = Letter)) + 
  xlab('Detection in angular Bennu samples') +
  ylab('LUCA clan usage') + 
  geom_text(aes(colour = factor(Moosmann_category)), size = 16) + 
  labs(color='Moosmann category') + 
  theme(axis.text=element_text(size=24),axis.title=element_text(size=32,face="bold"),
        legend.text=element_text(size=27), legend.title = element_text(size=16,face="bold")) + 
  scale_fill_brewer(palette="Dark2") +  theme(legend.position = 'none') +
  guides(color = guide_legend(override.aes = list(size = 16))) + 
  annotate(geom="text", y=.75, x=0.5, label=paste0("Weighted R2 = 0.25"), color="black", size = 14) +
  annotate(geom="text", y=.7, x=0.5, label="p = 0.01", color="black", size = 14) 


AngularvsPreLUCA <- ggplot(bennuvsusage_df, aes(y = as.numeric(Pre.LUCA.usage), 
                                             x = as.numeric(Bennu_presence$OREX08000550113..angular.), label = Letter)) + 
  xlab('Detection in angular Bennu samples') +
  ylab('Pre-LUCA clan usage') + 
  geom_text(aes(colour = factor(Moosmann_category)), size = 16) + 
  labs(color='Moosmann category') + 
  theme(axis.text=element_text(size=24),axis.title=element_text(size=32,face="bold"),
        legend.text=element_text(size=27), legend.title = element_text(size=16,face="bold")) + 
  scale_fill_brewer(palette="Dark2") +  theme(legend.position = 'none') +
  guides(color = guide_legend(override.aes = list(size = 16))) + 
  annotate(geom="text", y=.922, x=0.5, label=paste0("Weighted R2 = 0.18"), color="black", size = 14) +
  annotate(geom="text", y=.885, x=0.5, label="p = 0.03", color="black", size = 14) 

AggregatevsPreLUCA <- ggplot(bennuvsusage_df, aes(y = as.numeric(Pre.LUCA.usage), 
  x = as.numeric(Bennu_presence$OREX080010700...aggregate.), label = Letter)) + 
  xlab('Detection in aggregate Bennu samples') +
  ylab('Pre-LUCA clan usage') + 
  geom_text(aes(colour = factor(Moosmann_category)), size = 16) + 
  labs(color='Moosmann category') + 
  theme(axis.text=element_text(size=24),axis.title=element_text(size=32,face="bold"),
        legend.text=element_text(size=27), legend.title = element_text(size=16,face="bold")) + 
  scale_fill_brewer(palette="Dark2") +  theme(legend.position = 'none') +
  guides(color = guide_legend(override.aes = list(size = 16))) + 
  annotate(geom="text", y=.922, x=0.5, label=paste0("Weighted R2 = 0.25"), color="black", size = 14) +
  annotate(geom="text", y=.885, x=0.5, label="p = 0.01", color="black", size = 14) 

AggregatevsLUCA <- ggplot(bennuvsusage_df, aes(y = as.numeric(LUCA.usage), 
                  x = as.numeric(Bennu_presence$OREX080010700...aggregate.), label = Letter)) + 
  xlab('Detection in aggregate Bennu samples') +
  ylab('LUCA clan usage') + 
  geom_text(aes(colour = factor(Moosmann_category)), size = 16) + 
  labs(color='Moosmann category') + 
  theme(axis.text=element_text(size=24),axis.title=element_text(size=32,face="bold"),
  legend.text=element_text(size=27), legend.title = element_text(size=16,face="bold")) + 
  scale_fill_brewer(palette="Dark2") +  theme(legend.position = 'none') +
  guides(color = guide_legend(override.aes = list(size = 16))) + 
  annotate(geom="text", y=.75, x=0.5, label=paste0("Weighted R2 = 0.04"), color="black", size = 14) +
  annotate(geom="text", y=.7, x=0.5, label="p = 0.18", color="black", size = 14) 

gridExtra::grid.arrange(AngularvsLUCA,AngularvsPreLUCA , AggregatevsLUCA,AggregatevsPreLUCA )
