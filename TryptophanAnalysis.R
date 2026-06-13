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
library(phangorn)
setwd("~/Desktop/Tryptophan paper")
Tol_tree <- read.tree('Fixed_ToL.newick')
Distance2aaRS <- read.csv('AARSprotoenzyme_LUCA.csv', header = T)
PF00579_function <- read.tree('../Tryptophan paper/PF00579_function_geneTree.newick')
PF00579_aln <- read.fasta('../Tryptophan paper/PF00579_aln400.fasta', as.string = TRUE)
PF00579_conbac <- read.tree('../Tryptophan paper/PF00579_conHGTdropped_bac_branches.treefile')


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


reconciled.treefiles <- read.tree('YRS_recombination_arcbac.newick')
dropped.reconciled.treefiles <- drop.tip(reconciled.treefiles, 'G001940725_PF00579.Asgard.Archaea/1-7')
write.tree(unroot(dropped.reconciled.treefiles), 'YRS_recombination_arcbac_dropped.newick')
write.tree(unroot(reconciled.treefiles), 'YRS_recombination_arcbac_unrooted.newick')
write.tree(dropped_treefiles, '../Tryptophan paper/PF00579_arcbac_reconciled_HGTdropped.newick')
write.tree(midpoint(reconciled.treefiles), '../Tryptophan paper/PF00579_arcbac_reconciled_midpoint.newick')
Pfam_ConAAC[which(Pfam_ConAAC$pfamIDs == 'PF20974'),]


CRS_tree <- read.tree('../PF01406_enforceroot.treefile')
CRS_seq <- read.csv('../PFAM trees/PF01406_rows.csv', header = T)
CRS_seq$tiplabels <- paste0(CRS_seq$SequenceNames, '.', CRS_seq$Domain)
write.fasta( as.list(CRS_seq$pfam_sequence[which( CRS_seq$tiplabels %in% CRS_tree$tip.label)]),
             CRS_seq$tiplabels[which( CRS_seq$tiplabels %in% CRS_tree$tip.label)],
             '../PFAM trees/PF01406_seq.fasta')
             
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
UniprotPDBidmapping <- UniprotPDBidmapping[-which(
  grepl('CPR', UniprotPDBidmapping$SequenceNames)),]
UniprotPDBidmapping$Protein.names <- gsub('-.*','',UniprotPDBidmapping$Protein.names)
UniprotPDBidmapping$SequenceFuncNames <- paste0(UniprotPDBidmapping$accesion_names,'_', UniprotPDBidmapping$Entry,".", 
               UniprotPDBidmapping$Protein.names, '.', UniprotPDBidmapping$SequenceNames,".",  UniprotPDBidmapping$Domain)
UniprotPDBidmapping$SequenceFuncNames <- gsub('..........._PF00579', '', UniprotPDBidmapping$SequenceFuncNames)

UniprotPDBnoncanonYRS <- UniprotPDBidmapping[which(UniprotPDBidmapping$accesion_names %in%
                              gsub('_PF00579.*','',names(BacYRSnoncanon_aln)))  ,]
UniprotPDBcanonYRS <- UniprotPDBidmapping[which(UniprotPDBidmapping$accesion_names %in%
                   gsub('_PF00579.*','',names(BacYRScanon_aln)))  ,]
#UniprotPDBnononcanonYRS <- UniprotPDBnononcanonYRS[-which(grepl('CPR',UniprotPDBnononcanonYRS$SequenceFuncNames)),]


## Using I3IM11  as a ref (length 333AA) PFAM MATCH 4:287 / 124:661 IN 210 protein 3di alignment 
# 123: 653 in 209 3di aln with no long archaea G000403645
# with PDB (full length protein) muscle aln 186-926
# with PDB muscle aln and cysRS 177-1040
Protein3Di_seq <- seqinr::read.fasta('210PDBfiles/foldmason/PDB_3Di_MusclealnCysRS.fasta',as.string = T)
Truncated_alignedsequences <- substr(as.list(Protein3Di_seq),177,1040)
names(Truncated_alignedsequences ) <- names(Protein3Di_seq )
write.fasta(as.list(Truncated_alignedsequences), names(Truncated_alignedsequences),
            '210PDBfiles/foldmason/TruncPF00579_3Di_MusclealnCysRS.fasta' )
names(ProteinAA_seq) == names(Truncated_alignedsequences ) 
ProteinAA_seq <- read.fasta('210PDBfiles/foldmason/TruncPF00579_AA_MusclealnCysRS.fasta',as.string = T)
nchar(gsub('-','',ProteinAA_seq[which(grepl('Cys', names(ProteinAA_seq)))]))

write.fasta(as.list(paste0(ProteinAA_seq,Truncated_alignedsequences )), names(Truncated_alignedsequences ),
            '210PDBfiles/foldmason/TruncPF00579_AA3Di_MusclealnCysRS.fasta' )

TruncProtein3Di_seq <- seqinr::read.fasta('210PDBfiles/foldmason/TruncPF00579_3Di_MusclealnCysRS.fasta',as.string = T)
TruncProteinAA_seq <- seqinr::read.fasta('210PDBfiles/foldmason/TruncPF00579_AA_MusclealnCysRS.fasta',as.string = T)
Trim_Truncated_3Dialigned <- str_sub_all(as.list(TruncProtein3Di_seq),c(1,170,183,269,322,790,848),c(152,174,246,318,778,799,864))
Trim_Truncated_3Dialigned <- gsub(', ', "",lapply(Trim_Truncated_3Dialigned, toString))
Trim_Truncated_AAaligned <- str_sub_all(as.list(TruncProteinAA_seq),c(1,170,183,269,322,790,848),c(152,174,246,318,778,799,864)) 
Trim_Truncated_AAaligned <- gsub(', ', "",lapply(Trim_Truncated_AAaligned, toString))
write.fasta(as.list(Trim_Truncated_AAaligned) , names(TruncProteinAA_seq ), 
            '210PDBfiles/foldmason/Trimmed_TruncPF00579_AA_MusclealnCysRS.fasta')
write.fasta(as.list(paste0(Trim_Truncated_AAaligned,Trim_Truncated_3Dialigned)) , names(TruncProteinAA_seq ), 
            '210PDBfiles/foldmason/Trimmed_TruncPF00579_AA3Di_MusclealnCysRS.fasta')

nchar(gsub('-','',Trim_Truncated_AAaligned[which(grepl('Cys', names(TruncProteinAA_seq )))]))

#157-206,265-533, 558-589,633-659,674-709
TruncProtein3Di_seq <- seqinr::read.fasta('210PDBfiles/foldmason/TruncPF00579_3Di_Musclealn.fasta',as.string = T)
Trim_Truncated_AAaligned3 <- str_sub_all(as.list(TruncProtein3Di_seq),c(1,207,534,590,660,710),c(156,264,557,632,673,741)) 
Trim_Truncated_AAaligned3 <- gsub(', ', "",lapply(Trim_Truncated_AAaligned, toString))
write.fasta(as.list(Trim_Truncated_AAaligned) , names(TruncProtein3Di_seq), 
            '210PDBfiles/foldmason/Trimmed_TruncPF00579_AA_Musclealn.fasta')
write.fasta(as.list(paste0(Trim_Truncated_AAaligned,Trim_Truncated_AAaligned3)) , names(TruncProtein3Di_seq), 
            '210PDBfiles/foldmason/Trimmed_TruncPF00579_AA3Di_Musclealn.fasta')

Fig3C <- read.tree('../TruncPF00579_reconciledcon_bac.treefile')
SuppFig2 <- read.tree('../PF00579_alnCPRdropped_NQbac.contree')
0.33/sum(Fig3C$edge.length) # 0.3% of tree length
0.15/sum(SuppFig2$edge.length) # 0.08% of tree length

CPRdropped_reconbac <- read.tree('../PF00579_CPRdropped_reconbac.treefile')
CPRdropped_reconbac <- multi2di(CPRdropped_reconbac )
MADCPRdropped_reconbac <- mad(CPRdropped_reconbac, 'full')
MADCPRdropped_reconbac <- MADCPRdropped_reconbac[[6]][1][[1]]
MADCPRdropped_reconbac_252seq <- keep.tip(MADCPRdropped_reconbac, MADCPRdropped_reconbac$tip.label[which(gsub('_.*', "",MADCPRdropped_reconbac$tip.label ) %in% gsub('_.*', "",names(TruncProtein3Di_seq)) |
   gsub('_.*', "",MADCPRdropped_reconbac$tip.label ) %in% gsub('_.*', "",names(NonCanonYRS_recomb )))])


plot.phylo(MADCPRdropped_reconbac_252seq, show.tip.label = F)
write.tree(MADCPRdropped_reconbac_252seq, '../MADCPRdropped_reconbac_252seq.treefile')

trees <- read.tree('../PF00579_partition.rec_uml')
write.tree(midpoint.root(trees),'../MPTruncPF00579_AA_bac.treefile')
45*20

AA209_bacTree <- read.tree('../TruncPF00579_AA_bac.contree.treefile')
AA209_bacTree <- multi2di(AA209_bacTree)
MADAA209_bacTree <- mad(AA209_bacTree, 'full')
MADAA209_bacTree <- MADAA209_bacTree[[6]][1][[1]]
write.tree(MADAA209_bacTree, '../MADAA209_bacTree.treefile')
write.tree(midpoint.root(AA209_bacTree ), '../MPAA209_bacTree.treefile')

Rootstrap_str <- str_extract_all( toString( AA209_ArcTree[4]) ,  "\\[(.*?)\\]")
Rootstrap_str <- str_split(Rootstrap_str, ',')
Rootstrap_str <- as.numeric(str_extract(Rootstrap_str[[1]][grepl('rootstrap=', Rootstrap_str[[1]])], "\\d+"))
sum(Rootstrap_str) - 96- 96
51584.40 - 18026.00
FoldmasonTrunc_AAseq <- read.fasta('PF00579protein210_3DiAAalnNoG000403645.fa',as.string = T)
FoldmasonTrunc_AAseq <- FoldmasonTrunc_AAseq[match( names(Truncated_alignedsequences), names(FoldmasonTrunc_AAseq))]
nchar(gsub("-","", FoldmasonTrunc_AAseq )) - nchar(gsub("-","", Truncated_alignedsequences))



PF00579aln_3di <- seqinr::read.fasta('PF00579protein210_3DialnNoG000403645.fa', as.string=T)
PF00579aln_3diAA <- seqinr::read.fasta('PF00579protein210_3DiAAalnNoG000403645.fa', as.string=T)
PF00579aln_3diAA <- PF00579aln_3diAA[match(names(PF00579aln_3di), names(PF00579aln_3diAA))]
Musclealn3DI <- read.fasta('PF00579_3Di_Musclealn.fasta', as.string=T)
nchar(gsub('-','',Musclealn3DI [1]))


yrs_nonrecombtree <- read.tree('../YRSnonrecomb_reconciledcon_bac.treefile')
yrs_nonrecombtree <- multi2di(yrs_nonrecombtree )
MADyrs_nonrecombtree <- mad(yrs_nonrecombtree, 'full')
write.tree(MADyrs_nonrecombtree[[6]][[1]], '../YRSnonrecomb_reconciledcon_bac_MAD.treefile')

# 1:755 AA aln, 756:1286 3Di aln
write.fasta(as.list(paste0(PF00579aln_3diAA,PF00579aln_3di)), names(PF00579aln_3diAA),
            'PF00579_AA3Dialn.fasta')

arcyrs <- read.fasta('PF00579_ArcYRS.fasta')
EuryYRS <- gsub('G','GCF_',gsub('_.*','',names(arcyrs)[which(grepl('Euryarchaeota',names(arcyrs)))]) )
halobacteriales <- read.delim('Halobacteriales.tsv')
halobacteriales$Assembly.Accession <- gsub('[.].*','',halobacteriales$Assembly.Accession)
HaloYRS <- gsub('GCF_','G',EuryYRS[which(EuryYRS %in% halobacteriales$Assembly.Accession)])

write.fasta(arcyrs[which(gsub('_.*','',names(arcyrs)) %in% HaloYRS)],
            names(arcyrs)[which(gsub('_.*','',names(arcyrs)) %in% HaloYRS)],
            'HaloYRS.fasta')
length(arcyrs[[1]])
nchar(BacYRSnoncanon_aln[[1]])
aatree <- read.tree('PF00579_Musclepartition_conAA.newick')
aatree <- drop.tip(aatree, aatree$tip.label[which(gsub('_.*','',aatree$tip.label) %in% HaloYRS)])
plot.phylo(midpoint_root(aatree), show.tip.label = F)


## replace aa sites in muscle alignment with 3Di
PF00579aln_AA <- read.fasta('210PDBfiles/foldmason/Renamedfoldmason_full209aaseqCysRSaln.fa', as.string = F)
names(PF00579aln_AA)

PF00579aln_3di <- read.fasta('210PDBfiles/foldmason/foldmason_209ssseqCysRS.fa', as.string = F)
PF00579aln_3di <- PF00579aln_3di[match(gsub('[.].*', "", gsub('.*_', "" ,names(PF00579aln_AA)) ), names(PF00579aln_3di))]

#PF00579aln_3di <- PF00579aln_3diCys[match( names(PF00579aln_3diAACys),names(PF00579aln_3diCys))]
Musclealn_3Di <- list()
for ( seq in 1:211){
template_chars <-  PF00579aln_AA[[seq]]
#fill_chars <- PF00579aln_3di[[seq]][-which(PF00579aln_3di[[seq]] == '-')]
fill_chars <- PF00579aln_3di[[seq]]
fill_index <- 1
for (i in seq_along(template_chars)) {
  if (template_chars[i] != "-") {
    template_chars[i] <- fill_chars[fill_index]
    fill_index <- fill_index + 1 }}

Musclealn_3Di <- append (Musclealn_3Di, paste(template_chars, collapse = "")) }

nchar(PF00579aln_3diAAseqsCys)

PF00579aln_3diAAseqs <- sapply(1:211,function(x){paste(PF00579aln_AA[[x]], collapse = "")})
#PF00579aln_3diAAseqsCys <- sapply(1:211,function(x){paste(PF00579aln_3diAACys[[x]], collapse = "")})
write.fasta(as.list(PF00579aln_3diAAseqs), names(PF00579aln_AA),
            '210PDBfiles/foldmason/PDB_AA_MusclealnCysRS.fasta')
write.fasta(as.list(paste0(PF00579aln_3diAAseqs,Musclealn_3Di)), names(PF00579aln_AA),
            'PDB_AA3Di_MusclealnCysRS.fasta')
write.fasta(as.list(Musclealn_3Di), names(PF00579aln_AA),
            '210PDBfiles/foldmason/PDB_3Di_MusclealnCysRS.fasta')

YRS_UniprotID <- read.fasta('PF00579_AA3Di_MusclealnBacYRS.fasta')
YRS_UniprotIDnames <-  gsub('.*_', "",names(YRS_UniprotID))
gsub('.Tyrosine...*', "",YRS_UniprotIDnames)
accesion_line <- gsub(',', "",toString(accesion_Df[,2]))
paste('./dataformat download gene accesion', accesion_line )
      
JD_WAG <- read.tree('../iqtree-3.0.1-macOS/JD_WAG.contree')
JDrealign_WAG <- read.tree('../iqtree-3.0.1-macOS/TrpTyr_JDalign_noRecombaln.fasta.contree')
JD_WAGultra <- force.ultrametric(multi2di(JD_WAG))
JDrealign_WAGultra <- force.ultrametric(multi2di(JDrealign_WAG))
plot.phylo(midpoint.root(JD_WAGultra))
plot.phylo(midpoint.root(JDrealign_WAGultra))

bacstates <- readLines('../arcstates')
writeLines(Pfam_ConAAC$pfamIDs[-which(Pfam_ConAAC$pfamIDs %in% bacstates  )], '../MissingArc')
Pfam_ConAAC$ancestor[-which(Pfam_ConAAC$pfamIDs %in% bacstates  )]

plot.phylo(PF00579_enforce, show.tip.label = T)
PNAStree_noncanon <- extract.clade(PF00579_enforce, interactive = T)
PF00579_enforcenoncanondropped <- drop.tip(PF00579_enforce, PNAStree_noncanon$tip.label)
write.tree(PF00579_enforcenoncanondropped,'PF00579_enforcenoncanondropped.treefile')


GF_WAGG8 <- read.tree('../239_2015_9672_MOESM4_ESM.fasta.contree')
GF_WAGG8 <- midpoint.root(GF_WAGG8 )
plot.phylo(GF_WAGG8 , show.tip.label = F)
GF_WAGG8recombYRS <- extract.clade(GF_WAGG8 , interactive = T)
GF_WAGG8recombdropped <- drop.tip(GF_WAGG8 , GF_WAGG8recombYRS$tip.label)
plot.phylo(midpoint(GF_WAGG8recombdropped), show.tip.label = T)

GF_droppedbac <- read.tree('../JD47realign_NQbac.contree')
plot.phylo(midpoint.root(GF_droppedbac), show.tip.label = T)
MADGF_droppedbac <- mad(GF_droppedbac, 'full'  )
plot.phylo(MADGF_droppedbac[[6]][[1]], show.tip.label = T)
write.tree(MADGF_droppedbac[[6]][[1]],'../MAD_JD47realign_NQbac.contree')

GF_aln <- read.fasta('../239_2015_9672_MOESM4_ESM.fasta', as.string=T)
write.fasta(as.list(GF_aln[-which(names(GF_aln) %in% GF_WAGG8recombYRS$tip.label )]),
        names(GF_aln)[-which(names(GF_aln) %in% GF_WAGG8recombYRS$tip.label )],
        '../239_2015_9672_MOESM4_YRSrecombdropped.fasta' )
GF_alndropped <- read.fasta('../239_2015_9672_MOESM4_YRSrecombdropped.fasta',as.string = T)
write.fasta(as.list(str_remove_all(GF_alndropped,'-')), names(GF_alndropped),
            '../239_2015_9672_MOESM4_YRSrecombdropped_seq.fasta')

CysRS3Ditree <- read.tree('../TruncPF00579_reconciledconCysRS_3Di.treefile')
RootedCysRS3Ditree <- root.phylo(CysRS3Ditree, outgroup = 
                      CysRS3Ditree$tip.label[grepl('Cys',CysRS3Ditree$tip.label)], resolve.root = T)
RootedCysRS3Ditree$root.edge <- 1 
plot.phylo(RootedCysRS3Ditree, show.tip.label = F, use.edge.length = T , root.edge = 22)


# Assuming you have a list of trees named 'trees'
library(ggtree)
#parttrees <- read.tree('../PF00579_Musclepartition_con.parttrees')
pf00579_AAtree <- read.tree('../TruncPF00579_reconciledconCysRS_bac.treefile')
pf00579_AAtree <- unroot(pf00579_AAtree)
pf00579_AAtree$tip.label[which(grepl('Cys',pf00579_AAtree$tip.label))] <- 'CysRS_outgroup'
plot.phylo(pf00579_AAtree , show.tip.label = F)
pf00579_AAtree <- root.phylo(pf00579_AAtree, interactive = T)
pf00579_3ditree <- read.tree('../TruncPF00579_reconciledcon_3Di.treefile')
pf00579_3ditree <- unroot(pf00579_3ditree)

madtree  <- mad(multi2di(pf00579_AAtree), 'full')
madtree <- madtree[[6]][[1]]
write.tree(madtree,'../MAD_TruncPF00579_reconciledcon_arc.treefile' )

clsAA <-list(bacteria=grep('Bacteria',pf00579_AAtree$tip.label ),
          Archaea=grep('Archaea',pf00579_AAtree$tip.label ),
          cysteine=grep('Cys',pf00579_AAtree$tip.label ))
pf00579_AAtree_new <- groupOTU(pf00579_AAtree,.node=clsAA)


getMRCA(pf00579_3ditree, pf00579_3ditree$tip.label[which(
  grepl('Tyrosine', pf00579_3ditree$tip.label))])
#pf00579_3ditree_new$node <- seq(1:pf00579_3ditree_new$Nnode)
pf00579_3ditree_new <- rotateNodes(pf00579_3ditree, 280)
plot.phylo(pf00579_3ditree, show.tip.label = F, 'unrooted') 
plot.phylo(pf00579_3ditree_new, show.tip.label = F, 'unrooted') 

cls3Di <-list(bacteria=grep('Bacteria',pf00579_3ditree_new$tip.label ),
              Archaea=grep('Archaea',pf00579_3ditree_new$tip.label ))
pf00579_3ditree_new <- groupOTU(pf00579_3ditree_new ,.node=cls3Di)


ggtree( c(pf00579_3ditree_new, pf00579_AAtree_new) , root.position = 0) + 
  geom_tree(aes(color=group), size =0.8) +  
  facet_wrap(~.id, scales = "fixed") + 
  geom_treescale(width = 1.0, x=-1.5, , y= 4, offset = 0.2, linesize=0.8, fontsize = 8) +
  #geom_tree(aes(color=group), size =0.8) +  
  geom_rootedge(1, size= 0.8) +
  scale_color_manual(values=c("blue","red", 'green', 'black')) + 
  theme(legend.position = 'none', strip.background = element_blank(), 
   strip.text = element_blank())

1.4244+1.8675
6.1 = 3.2919
  = 1.8675 #NQ
 = 2.7268 # midpoint
 = 2.4239 # MAD

(2.4239*6.1)/3.2919 
max_x <- max(ggplot_build(tree3di_plot )$layout$panel_scales_x[[1]]$range$range)
p1_scaled <- tree3di_plot  + xlim(-2, max_x * 20) # Extend a bit for spacing
p2_scaled <- treeAA_plot + xlim(-2, max_x * 20)

tree3di_plot <- ggtree( pf00579_3ditree_new , root.position = 2, , layout="equal_angle") + 
 # geom_treescale(width = 0.5, x=-0.5, y= -0.5, offset = 0.05, linesize=0.8) +
  geom_tree(aes(color=group), layout="equal_angle", size =0.8) +  
  scale_color_manual(values=c("blue","red")) + 
  theme(legend.position = 'none', strip.background = element_blank(), 
        strip.text = element_blank())

treeAA_plot <- ggtree(  pf00579_AAtree_new , root.position = 2) + 
 # geom_treescale(width = 0.5, x=2, offset = 2, linesize=0.8) +
  geom_tree(aes(color=group), size =0.8) +  geom_rootedge(2) +
  scale_color_manual(values=c("blue","red")) + 
  theme(legend.position = 'none', strip.background = element_blank(), 
        strip.text = element_blank())

ggtree(c(p1_scaled, p2_scaled)) 




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
# Archaeal YRS node is 3
# New tree: NODE 2 is archaeal YRS + WRS, Node 3 is archaeal YRS, 
#Node 61 is WRS, Node 175 is bacterial YRS
state_tree <- read.tree('../iqtree-3.0.1-macOS/TruncPF00579_AA_Musclealn.fasta.treefile')
plot.phylo(state_tree, show.node.label = T, show.tip.label = F)
#state_file <- read.delim('../iqtree-3.0.1-macOS/bin/TruncPF00579_ASR_Bac.state', comment.char = "#")
state_file <- read.delim('../iqtree-3.0.1-macOS/TruncPF00579_AA_Musclealn.fasta.state', comment.char = "#")#nq.pfam
root_state <- state_file[which(state_file$Node == 'Node1'),]
root_state[which(root_state$State == 'W'),] #3
state_file[which(state_file$Node == 'Node61' & state_file$State == 'W'),]
state_file[which(state_file$Node == 'Node2' & state_file$State == 'W'),] #6
state_file[which(state_file$Node == 'Node3' & state_file$State == 'Y'),] #9, 12
state_file[which(state_file$Node == 'Node175' & state_file$State == 'W'),] #10,9
root_state[which(root_state$Site == 55),]
# nq.bac , nq.arc and nq.pfam all have 4 W sites at 174 >99%
# nq.bac shows 4 W sites at root with ~50%

state_file[which(state_file$Node == 'Node175' & state_file$Site ==606),] #F(0.59), R(0.99), F(0.84), Y(0.99)

# AT THE ROOT THERE ARE 5 sites that are most likely W but with very low probabilities <50%
# 37,38,33,39,40. The second most like aa at these site are hydrophobic F, L, L, F,Y respectively
sort(root_state[which(root_state$State == 'W'),][1,4:23], decreasing=T)
state_file[which(state_file$Node == 'Node2' & state_file$Site == 134),]

anc_Seq <- state_file$State[which(state_file$Node == 'Node175')]
names(anc_Seq ) <- 1:length(anc_Seq )
which(names(anc_Seq)[-which(anc_Seq == '-')] == '606')
anc_Seq[-which(anc_Seq == '-')][which(names(anc_Seq)[-which(anc_Seq == '-')] == '127')]
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

W_sitesproblist <- readRDS("W_sitesproblist.rds")
W_sitesproblist[1:5]
W_LACA <- unlist(replace(W_sitesproblist[which(names(W_sitesproblist) == 'LACA')],
                         lengths(W_sitesproblist[which(names(W_sitesproblist) == 'LACA')]) == 0, 0))
length(which(W_LACA > 0.9 ))/length(W_LACA ) # 49%
which(W_LACA == 1 )
W_LBCA <- unlist(replace(W_sitesproblist[which(names(W_sitesproblist) == 'LBCA')],
                         lengths(W_sitesproblist[which(names(W_sitesproblist) == 'LBCA')]) == 0, 0))
length(which(W_LBCA < 0.9 ))/length(W_LBCA ) # 67%

W_LUCA <- unlist(replace(W_sitesproblist[which(names(W_sitesproblist) == 'LUCA')],
                         lengths(W_sitesproblist[which(names(W_sitesproblist) == 'LUCA')]) == 0, 0))
length(which(W_LUCA == 1))/length(which(names(W_sitesproblist) == 'LUCA')) #39%

W_preLUCA <- unlist(replace(W_sitesproblist[which(names(W_sitesproblist) == 'preLUCA')],
                            lengths(W_sitesproblist[which(names(W_sitesproblist) == 'preLUCA')]) == 0, 0))
length(which(W_preLUCA == 1))/length(which(names(W_sitesproblist) == 'preLUCA')) #32%

W_sitesproblist <- readRDS("W_sitesproblistPfamIDs0.4_arc.rds")
W_preLACA <- unlist(replace(W_sitesproblist[which(names(W_sitesproblist) %in% AncientPostLUCA$PFAM_IDs[which(AncientPostLUCA$classified_ancestor == 'preLACA')])],
                         lengths(W_sitesproblist[which(names(W_sitesproblist) %in% AncientPostLUCA$PFAM_IDs[which(AncientPostLUCA$classified_ancestor == 'preLACA')])]) == 0, 0))
names(W_preLACA ) <- str_extract(names(W_preLACA ),'.......')
length(which(W_preLACA < 0.7))/ length(W_preLACA) # 46%
length(which(W_preLACA >= 0.7))
unique(names(W_preLACA )[which(W_preLACA > 0.9)])
# 31 total, 11 prelaca clans with confident W: 2 Clans and 9 pfams
diversifiedLACAclans[which(diversifiedLACAclans %in% 
 unique(Pfam_ConAAC$clans[which(Pfam_ConAAC$pfamIDs  %in% names(W_preLACA)[which(W_preLACA > 0.9)]  )] )  ) ]

Pfam_ConAAC[which(Pfam_ConAAC$clans == 'CL0734'),]
# CL0539 has 3 pfams 1 unclassifiable, one prelaca PF11469 and laca; CL0734 1 prelaca PF01900

W_sitesproblist[[which(names(W_sitesproblist) == 'PF01984')]]
#PreLACA not NA interacting CL0229, CL0200, CL0639?,CL0402,CL0214, CL0668
W_preLBCA <- unlist(replace(W_sitesproblist[which(names(W_sitesproblist) %in% AncientPostLUCA$PFAM_IDs[which(AncientPostLUCA$classified_ancestor == 'preLBCA')])],
                         lengths(W_sitesproblist[which(names(W_sitesproblist) %in% AncientPostLUCA$PFAM_IDs[which(AncientPostLUCA$classified_ancestor == 'preLBCA')])]) == 0, 0))
length(which(W_preLBCA < 0.7))/ length(W_preLBCA) #70%
diversifiedLBCAclans[which(diversifiedLBCAclans %in% 
     unique(Pfam_ConAAC$clans[which(Pfam_ConAAC$pfamIDs  %in% names(W_preLBCA)[which(W_preLBCA > 0.9)]  )] )  ) ]
Clan_ConAAfreq_Df$Clan_ancestor[which(Clan_ConAAfreq_Df$Clans %in% diversifiedLBCAclans[which(diversifiedLBCAclans %in% 
                  unique(Pfam_ConAAC$clans[which(Pfam_ConAAC$pfamIDs  %in% names(W_preLBCA)[which(W_preLBCA > 0.9)]  )] )  ) ])]
PreLACAecdf_fun <- ecdf(W_preLACA )
PreLACA_contr <- PreLACAecdf_fun(W_preLACA )*W_preLACA 
PreLBCAecdf_fun <- ecdf(W_preLBCA )
PreLBCA_contr <- PreLBCAecdf_fun(W_preLBCA )*W_preLBCA 
W_ecdf <- bind_rows(
  data.frame(val = W_preLACA, contribution= PreLACA_contr,group = "preLACA"),
  data.frame(val = W_preLBCA, contribution= PreLBCA_contr, group = "preLBCA"))
# PF04032 double check W, PF04034 (one R), PF05693,PF11469,PF05559? CONSERVED W

conserved_W <- read.csv('Conserved_W.csv', header = T)

Y_sitesproblist <- readRDS("Y_sitesproblistPfamIDs0.4_arc.rds")
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
# PF04032 double check W, PF04034 (one R), PF05693,PF11469,PF05559? CONSERVED W

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


library(AMR)
contigency_table <- matrix(c(32,1, 745, 155), nrow = 2)
g.test(contigency_table )
contigency_table <- matrix(c(3,0, 79, 74), nrow = 2)
g.test(contigency_table )
contigency_table <- matrix(c(35,1, 79+745, 74+155), nrow = 2)
fisher.test(contigency_table )


#### YRS recombination analysis#####
all_aln_fasta <- seqinr::read.fasta('../Tryptophan paper/PF00579_alnCPRdropped.fasta', as.string = T)
aln_fasta <- seqinr::read.fasta('../Tryptophan paper/YRS_aln.fasta', as.string = T)
names(aln_fasta) <- gsub('/.*','', names(aln_fasta))
CPRdropped_arctree <- read.tree('../Tryptophan paper/PF00579_arc_CPRdropped.treefile')
plot.phylo(CPRdropped_arctree, show.tip.label = T)
BacYRScanon <- extract.clade(CPRdropped_arctree, interactive = T)
BacYRSnoncanon <- extract.clade(CPRdropped_arctree, interactive = T)
ArcYRS <- extract.clade(CPRdropped_arctree, interactive = T)
ArcWRS <- extract.clade(CPRdropped_arctree, interactive = T)
BacWRS <- extract.clade(CPRdropped_arctree, interactive = T)
AllBacYRS <- extract.clade(CPRdropped_arctree, interactive = T)
# HPLDLK region in weblogo 610-625 ... 413-430
# 55-75 INSERTIONS AT N TERMINAL 32-52
#write.tree(ArcYRS, 'ArcYRS.treefile')
PF00579_enforce <- read.tree('../Tryptophan paper/PF00579_enforceroot.treefile')
PF00579_enforce <- drop.tip(PF00579_enforce, 
                   c(PF00579_enforce$tip.label[which(grepl('CPR',PF00579_enforce$tip.label))],
                    "G000403645_PF00579.Euryarchaeota.Archaea"))
write.tree(PF00579_enforce, '../Tryptophan paper/PF00579_enforcerootdropped.treefile')
                           
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
# after removing all - from PF00579_YRSrecombination.fasta
yrs_recon <- seqinr::read.fasta('../Tryptophan paper/PF00579_YRSrecombination.fasta', as.string = T)
write.fasta(as.list(yrs_recon[-which(nchar(yrs_recon) < 90)]), 
            names(yrs_recon)[-which(nchar(yrs_recon) < 90)], '../Tryptophan paper/PF00579_YRSrecombination162.fasta')

YRSrecomb_aln <- read.fasta( 'PF00579_YRSrecombination.fasta', as.string = T)
nonrecomb <- read.fasta('PF00579_YRSnonrecombination.fasta', as.string = T)
length(which(nchar(YRSrecomb_aln) >= 90))
median(nchar(nonrecomb))
median(nchar(YRSrecomb_aln))
which(grepl('G000147335_PF00579.Terrabacteria.Bacteria.noncanon',names(YRSrecomb_aln)))
which(grepl('G000143845_PF00579.Terrabacteria.Bacteria.noncanon',names(YRSrecomb_aln)))
nchar(YRSrecomb_aln)[1]
#44=23, 73=53, 61=56
length(which(nchar(NonCanonYRS_recomb) > 90))
NonCanonYRS_recomb <- read.fasta('NonCanon_BacYRS_recombination.fasta',as.string = T) #hmm


YRS_CRSaln <- read.fasta('YRSprofile_alnPF01406.fasta', as.string = T)
YRS_recomb <-  read.fasta('../Tryptophan paper/YRSprofile_alnCRSrecombination.fasta', as.string = T)
seqinr::write.fasta(as.list(str_sub(YRS_CRSaln, 268, 476)), names(YRS_CRSaln), 
                    '../Tryptophan paper/YRSprofile_alnPF01406recombination.fasta')

seqinr::write.fasta(as.list(str_sub(YRS_CRSaln, 448, 678)), names(YRS_CRSaln), 
                    '../Tryptophan paper/PF00579_YRS_PF01406_PF00750recombination.fasta')
seqinr::write.fasta(as.list(lapply(str_sub_all(YRS_CRSaln , c(30,679), c(447,690)), toString)),
                    names(YRS_CRSaln), '../Tryptophan paper/PF00579_YRS_PF01406_PF00750nonrecombination.fasta')

seqinr::write.fasta(as.list(c(str_sub(BacYRSnoncanon_aln , 214, 415), str_sub(BacYRScanon_aln , 214, 415), str_sub(ArcYRS_aln , 214, 415))),
                    c(names(BacYRSnoncanon_aln), names(BacYRScanon_aln), names(ArcYRS_aln)),
                    '../Tryptophan paper/PF00579_YRSrecombinationCRS.fasta')

seqinr::write.fasta(as.list(c(lapply(str_sub_all(BacYRSnoncanon_aln , c(1,420), c(213,448)), toString),
                              lapply(str_sub_all(BacYRScanon_aln , c(1,420), c(213,448)), toString),
                              lapply(str_sub_all(ArcYRS_aln , c(1,420), c(213,448)), toString))),
                    c(names(BacYRSnoncanon_aln), names(BacYRScanon_aln), names(ArcYRS_aln)),
                    '../Tryptophan paper/PF00579_YRSnonrecombination.fasta')

seqinr::write.fasta(as.list(c(lapply(str_sub_all(BacYRSnoncanon_aln , c(1,420), c(213,430)), toString),
                lapply(str_sub_all(BacYRScanon_aln , c(1,420), c(213,430)), toString),
                 lapply(str_sub_all(ArcYRS_aln , c(1,420), c(213,430)), toString))),
                    c(names(BacYRSnoncanon_aln), names(BacYRScanon_aln), names(ArcYRS_aln)),
                    '../Tryptophan paper/PF00579_YRSnonrecombination1_430.fasta')

seqinr::write.fasta(as.list(c(lapply(str_sub_all(BacYRSnoncanon_aln , 1, 213), toString),
                              lapply(str_sub_all(BacYRScanon_aln ,1, 213 ), toString),
                              lapply(str_sub_all(ArcYRS_aln , 1, 213), toString))),
                    c(names(BacYRSnoncanon_aln), names(BacYRScanon_aln), names(ArcYRS_aln)),
                    '../Tryptophan paper/PF00579_YRSnonrecombination1_213.fasta')

AllPF00579_noCPR_labelled <- seqinr::read.fasta('../Tryptophan paper/AllPF00579_noCPR_labelled.fasta',as.string = T)
seqinr::write.fasta(as.list(str_sub(AllPF00579_noCPR_labelled , 286, 620)),
                    names(AllPF00579_noCPR_labelled),
                    '../Tryptophan paper/AllPF00579_noCPR_labelled_recombination.fasta')
WithCysRS <- read.fasta('foldmason (4)/foldmason_ss.fa', as.string = T) 
# 25-396 in foldmason aln, 25-388 in muscle aln 
seqinr::write.fasta(as.list(str_sub(WithCysRS, 25, 396)),
                    names(WithCysRS),
                    'foldmason (4)/foldmason_sstruncated.fasta')

aa_sequences <- readFASTA('../Tryptophan paper/NonCanon_BacYRS_recombination.fasta')
aa_sequences_abin <- toupper(aa_sequences )
names(aa_sequences_abin) <- names(aa_sequences)
RecombinationHMM <- derivePHMM(aa_sequences_abin , residues='AMINO')
writePHMM(RecombinationHMM , 'YRSrecombination.hmm')
UniprotPDBidmapping$Organism[grepl('pyrococcus', UniprotPDBidmapping$Organism)]
pf00579_matching$V5[grepl('Bacillus subtilis', pf00579_matching$V5)]
bacillus_subtilis <- read.delim('Bacillus_subtilis.tsv', header = T)
bacillus_subtilis_Acc <- gsub('GCA_','G',bacillus_subtilis$Assembly.Accession)
bacillus_subtilis_Acc <- gsub('GCF_','G',bacillus_subtilis_Acc )
bacillus_subtilis_Acc <- gsub('[.].*','',bacillus_subtilis_Acc )
Tol_tree$tip.label[which(Tol_tree$tip.label %in% bacillus_subtilis_Acc)]

sample1 <- sample(44, 10)
sample2 <- sample(73, 10)
sample3 <- sample(61, 10)
seqinr::write.fasta(as.list(c(BacYRSnoncanon_aln[sample1],BacYRScanon_aln[sample2],
     ArcYRS_aln[sample3] )), c(names(BacYRSnoncanon_aln)[sample1], 
     names(BacYRScanon_aln)[sample2],names(ArcYRS_aln)[sample3] ), 'YRS10sample2_GARD.fasta')

YRS_PD <- read.csv('YRS_PD.csv', header = F)
GARD_YRS10 <- read.fasta('YRS10sample_GARD.fasta', as.string = T)
seqinr::write.fasta(as.list(str_sub(GARD_YRS10 , 1, 430)),
                   names(GARD_YRS10),  '../Tryptophan paper/YRS10sample430_GARD.fasta')

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
1 - (mean(AllBacdist))^2 # 66% similarity 30%


#2.82 = 3.0827 aa sub UNROOTED
#2.5935 = 2.8351 aa sub MP
#1.414894 = 1.5467 AA sub NQ
#2.21240 = 2.4185 AA sub MAD



# 0.88 = 1.087 3DI sub unrooted
# 0.8309402 = 1.0264 3DI sub MP
#  0.6937994 = 0.857 3DI SUB MAD


### Color bacterial tree with non canon vs canon YRS ####
BacYRStips <- c(names(BacYRScanon_aln),names(BacYRSnoncanon_aln))
bacterialYRSspeciestree <- keep.tip(Tol_tree, Tol_tree$tip.label[c(
  which(Tol_tree$tip.label %in% gsub('_.*','',names(BacYRScanon_aln)) ),
  which(Tol_tree$tip.label %in% gsub('_.*','',names(BacYRSnoncanon_aln)) ))])

UniprotPDBnoncanonYRS[which(UniprotPDBnoncanonYRS$accesion_names =="G000409775"),]
bacterialYRSspeciestips <- BacYRStips[match(bacterialYRSspeciestree$tip.label,gsub('_.*','',BacYRStips))]
bacterialYRSspeciestips <- gsub('.*_PF00579.','',bacterialYRSspeciestips)
bacterialYRSspeciestips <- gsub('.Bacteria','',bacterialYRSspeciestips)

colortaxa <- rep('black', nrow(bacterialYRSspeciestree$edge))
colortaxa[which(bacterialYRSspeciestree$edge[,2] %in% which(bacterialYRSspeciestree$tip.label %in%  gsub('_.*','',names(BacYRScanon_aln)) ))] <- 'red'
colortaxa[which(bacterialYRSspeciestree$edge[,2] %in% which(bacterialYRSspeciestree$tip.label %in%  gsub('_.*','',names(BacYRSnoncanon_aln)) ))]  <- 'magenta'
bacterialYRSspeciestree$tip.label <- bacterialYRSspeciestips
plot.phylo(bacterialYRSspeciestree , show.tip.label = T,type= 'phylogram', 
           edge.color = colortaxa, edge.width = 2 , align.tip.label = T) +
  add.scale.bar(x = 0.5, y = -1,cex =2, font = 15, lwd =2)


# info on canon vs non canon and number of pfams in protein
AllBacYRS_S4tree <- keep.tip(AllBacYRS, AllBacYRS$tip.label[which(
  AllBacYRS$tip.label %in% c(paste0(UniprotPDBcanonYRS$SequenceNames,'.', UniprotPDBcanonYRS$Domain),
                             paste0(UniprotPDBnoncanonYRS$SequenceNames,'.', UniprotPDBnoncanonYRS$Domain)) )])
BacYRSnbofPfams <- nchar(c(UniprotPDBcanonYRS$Pfam, UniprotPDBnoncanonYRS$Pfam))
library(RRphylo)
#AllBacYRS_S4tree <-  multi2di(AllBacYRS_S4tree)
AllBacYRS_S4tree <- force.ultrametric(AllBacYRS_S4tree)
AllBacYRS_S4tree <- fix.poly(AllBacYRS_S4tree, type = "resolve")
tipshape <- rep('NA', length(AllBacYRS_S4tree$tip.label ))
tipshape[which(AllBacYRS_S4tree$tip.label %in% 
       c(paste0(UniprotPDBcanonYRS$SequenceNames,'.', UniprotPDBcanonYRS$Domain),
       paste0(UniprotPDBnoncanonYRS$SequenceNames,'.', UniprotPDBnoncanonYRS$Domain))[
       which( BacYRSnbofPfams > 10)] )] <- 16
AllBacYRS_S4tree$tip.label [39]
S4tree_color <- rep('black', nrow(AllBacYRS_S4tree$edge))
S4tree_color[1:63] <- 'red'
S4tree_color[64:110] <- 'magenta'
AllBacYRS_S4tree$root.edge <-  1

153-169
175-182
247-268 
779-789
800-847
16+7+21+10+47
103/6
plot.phylo(AllBacYRS_S4tree, show.tip.label = T,type= 'phylogram', 
     edge.width = 2 , align.tip.label = T, edge.color = S4tree_color, root.edge = T ) +
  tiplabels(pch = as.numeric(tipshape) , adj = 0.5, cex = 1.2) +
  add.scale.bar(x = 0.5, y = -1,cex =2, font = 15, lwd =2) 

which.min(UniprotPDBidmapping$Length)
UniprotPDBidmapping[128,]
yrsnonrecomb430 <- read.tree('../YRS_nonrecomb430rec_bac.treefile')
MAD_yrsnonrecomb430 <- mad(yrsnonrecomb430 , 'full')
write.tree(MAD_yrsnonrecomb430[[6]][[1]], '../MAD_YRS_nonrecomb430rec_bac.treefile')
4 =1.93
2.17= 1.05 #MAD
2.32= 1.12 #MP

bacterialYRSspeciestree$tip.label[which(bacterialYRSspeciestree$tip.label %in%gsub('_.*','',names(BacYRScanon_aln)) )] <- 
  paste0(bacterialYRSspeciestree$tip.label[which(bacterialYRSspeciestree$tip.label %in%gsub('_.*','',names(BacYRScanon_aln)) )], 'YRS_canon')
bacterialYRSspeciestree$tip.label[which(bacterialYRSspeciestree$tip.label %in%gsub('_.*','',names(BacYRSnoncanon_aln)) )] <- 
  paste0(bacterialYRSspeciestree$tip.label[which(bacterialYRSspeciestree$tip.label %in%gsub('_.*','',names(BacYRSnoncanon_aln)) )], 'YRS_noncanon')
write.tree(bacterialYRSspeciestree, 'BacterialYRSspeciestree.newick')

PF00579_209AAarc <- read.tree('../PF00579_209AA_rec_arcrootstrap.rootstrap.nex')
PF00579_AANQ <- read.tree('../NQbranches_part_recon_contree.treefile')
PF00579_209AAarc <- force.ultrametric(PF00579_209AAarc$`treetree_1=`)
PF00579_AANQ <- force.ultrametric(PF00579_AANQ)
plot.phylo(midpoint(PF00579_AANQ), show.tip.label = F)



##### Plot roots on unrooted tree #####
YRSrecomb_conbac <- read.tree('../YRS_recombination132_reconciledarcbaccons_bac.treefile')
#YRSrecomb_conbac <- read.tree('../PF00579_YRSnonrecombination_arcbacreconciled_bac.treefile')
mp_YRSrecomb_conbac <- phangorn::midpoint(YRSrecomb_conbac )
unrootedtree <- unroot.phylo(mp_YRSrecomb_conbac )
MADroot <- mad(unrootedtree, 'full' )
MADrootree <-MADroot[[6]][[1]]

sum(mp_YRSrecomb_conbac$edge.length[which(mp_YRSrecomb_conbac$edge[,1] == getRoot(mp_YRSrecomb_conbac))])
#2.1939350/2.83947 # recombinant, longer branch is towards canon bacteria 77%
#1.4721509/2.421396 # nonrecombinant, longer branch is towards archaea 60%

sum(MADrootree$edge.length[which(MADrootree$edge[,1] == getRoot(MADrootree))])
#1.569889/2.83947 # recombinant, longer branch is towards canon bacteria 55%
#1.341029/2.421396 # nonrecombinant, longer branch is towards archaea 55%

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
