library(dada2)
packageVersion("dada2")

path1 <- "ANT17_project" 
list.files(path1) # Verify the file list

fnFs1 <- sort(list.files(path1, pattern="_R1_001.fastq", full.names = TRUE))
fnRs1 <- sort(list.files(path1, pattern="_R2_001.fastq", full.names = TRUE))

sample.names1 <- sapply(strsplit(basename(fnFs1), "_"), `[`, 1)
sample.names1

filtFs1 <- file.path(path1, "filtered", paste0(sample.names1, "_F_filt.fastq.gz"))
filtRs1 <- file.path(path1, "filtered", paste0(sample.names1, "_R_filt.fastq.gz"))
	

png("ANT17_project/plots/QplotF.png", width=1600, height=1200)
plotQualityProfile(fnFs1[3:8]) # Forward sequences
dev.off()

png("ANT17_project/plots/QplotR.png", width=1600, height=1200)
plotQualityProfile(fnRs1[3:8]) # reverse sequences
dev.off()

out1 <- filterAndTrim(fnFs1, filtFs1, fnRs1, filtRs1, truncLen=c(290,220), maxN=0, maxEE=c(6,7), truncQ=2, rm.phix=TRUE, trimLeft=c(20,21), compress=TRUE, multithread=TRUE)
head(out1)

errF1 <- learnErrors(filtFs1, multithread=TRUE, randomize=TRUE, verbose=TRUE)

errR1 <- learnErrors(filtRs1, multithread=TRUE, randomize=TRUE, verbose=TRUE)

png("ANT17_project/plots/errF.png", width=1600, height=1200)
plotErrors(errF1, nominalQ=TRUE)
dev.off()

png("ANT17_project/plots/errR.png", width=1600, height=1200)
plotErrors(errR1, nominalQ=TRUE)
dev.off()

# Dereplication step
derepFs1 <- derepFastq(filtFs1, verbose=TRUE)

derepRs1 <- derepFastq(filtRs1, verbose=TRUE)

names(derepFs1) <- sample.names1
names(derepRs1) <- sample.names1

dadaFs1 <- dada(derepFs1, err=errF1, pool=TRUE, multithread=TRUE)
dadaRs1 <- dada(derepRs1, err=errR1, pool=TRUE, multithread=TRUE)

bim <-isBimeraDenovo(seqtab2, minFoldParentOverAbundance = 8, allowOneOff = TRUE, verbose = FALSE)

# Mate Pairing the reads
mergers1 <- mergePairs(dadaFs1, derepFs1, dadaRs1, derepRs1, verbose=TRUE)

save.image(file="ANT17_project/RData/mergers1.RData")

seqtab1 <- makeSequenceTable(mergers1)
dim(seqtab1) # Gives you info on the number of Samples and ASVs identified in run_1

table(nchar(getSequences(seqtab1)))

#seqtab2 <- seqtab1[,nchar(colnames(seqtab1)) %in% seq(363,375)]

table(nchar(getSequences(seqtab2)))
dim(seqtab2) # Gives you info on the number of Samples and ASVs after tail sequence dropping

getN <- function(x) sum(getUniques(x))
track1 <- cbind(out1, sapply(dadaFs1, getN), sapply(dadaRs1, getN), sapply(mergers1, getN))
colnames(track1) <- c("input", "filtered", "denoisedF", "denoisedR", "merged")
rownames(track1) <- sample.names1
head(track1) 
write.csv(track1, "ANT17_project/csv/run1_asv_stats.csv")

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab2, minFoldParentOverAbundance = 8, method="consensus", multithread=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab2) # Gives you the percentage of sequences recovered after chimera removal
## This is a ggod check points. Even if a lot of ASVs have been removed, the majority of reads sould
## be retained. Usually >0.80 (aka 80%) are retained

## Assign Taxonomy. Point to where the silva database actually is
taxa <- assignTaxonomy(seqtab.nochim, "silva_db/silva_nr_v138_train_set.fa.gz", multithread=TRUE)
#taxa <- addSpecies(taxa, "silva_db/silva_species_assignment_v138.fa.gz")
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
write.csv(taxa.print, "ANT17_project/csv/taxa_print.csv") # For inspection in bash or excel
save.image(file="ANT17_project/RData/taxa_assigned.RData")

############ Making Tree
seqs <- getSequences(seqtab.nochim) ## Get the sequences
names(seqs) <- seqs
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA) # Generating the multiple sequence alignment
writeXStringSet(alignment, "alignment.fasta", format="fasta") # Exporting the alignment to an external file

## Build the tree using FastTree in bash using the GTR model. See more detail about
## FastTree at http://www.microbesonline.org/fasttree/

system('fasttree -gtr -nt alignment.fasta > alignment.tree', intern=TRUE)
tree <- read.tree(file = "alignment.tree") # Reading back the tree into the R environment

## Sanity check. The two numbers should match!
dim(seqtab.nochim)
tree$Nnode+2

save.image(file="ANT17_project/RData/tree.RData")

saveRDS(taxa,"ANT17_taxa.rds")
saveRDS(seqtab.nochim,"ANT17_seqtab_nochim.rds")
