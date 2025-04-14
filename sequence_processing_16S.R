#16S metabarcoding sequence processing 
#R.K.Salis 
BiocManager::install(version = '3.14')
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2", version = "3.14")
library(dada2); packageVersion("dada2")
library(phyloseq)
library(seqinr)
library(DESeq2)

path <- "all" # directory containing the fastq files
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#check quality
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

#Filter and trim
#Assign the filenames and place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(290,250),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) 
head(out)
#Learn the Error Rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
#Sample Inference
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaFs[[1]]

#Merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

#Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

#Remove chimeras
seqtab.nochim.16S <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim.16S)
sum(seqtab.nochim.16S)/sum(seqtab)
saveRDS(seqtab.nochim.16S, "seqtab.nochim.16S.rds")

#Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track16S <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track16S) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track16S) <- sample.names
head(track16S)
saveRDS(track16S, "track16S.rds")
write.csv(track16S, "track16S.csv")

#Assign taxonomy
taxa.16S <- assignTaxonomy(seqtab.nochim.16S, "silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE, tryRC=TRUE)
saveRDS(taxa.16S, "taxa.16S.rds")
taxa_sp.16S <- addSpecies(taxa.16S, "silva_species_assignment_v138.1.fa.gz")
saveRDS(taxa_sp.16S, "taxa_sp.16S.rds")

taxa.print <- taxa.16S 
rownames(taxa.16S) <- NULL
head(taxa.16S)
taxa.print <- taxa_sp.16S 
rownames(taxa_sp.16S) <- NULL
head(taxa_sp.16S)


#create fasta file
seqnum <- paste0("ASV", seq(ncol(seqtab.nochim)))
uniqueSeqs <- as.list(colnames(seqtab.nochim))
write.fasta(uniqueSeqs, seqnum, "16Sfulldataset.fasta")

#create phyloseq object
eDNA1.16S = phyloseq(tax_table(taxa.sp.16S), otu_table(seqtab.nochim.16S, taxa_are_rows = FALSE))
eDNA1.16S
rownames(otu_table(eDNA1.16S))
rownames(tax_table(eDNA1.16S))
## Upload metadata ##
map.eDNA1.16S1 = read.csv("mapping_eDNA1_16S.csv")
head(map.eDNA1.16S1)
###Rename samples for metadata
map.eDNA1.16S <- map.eDNA1.16S1[, -1]
row.names(map.eDNA1.16S)<- map.eDNA1.16S1$Sample #sample data row names must align with dada2 rowname outputs
map.eDNA1.16S = as.data.frame(map.eDNA1.16S)
head(map.eDNA1.16S)
eDNA1.16S = phyloseq(tax_table(taxa.sp.16S), otu_table(seqtab.nochim.16S, taxa_are_rows = FALSE), sample_data(map.eDNA1.16S))
eDNA1.16S
###Rename sequence variants
a.vec = as.vector(1:22143)  #total ASVs
a.nam = cbind("ASV", a.vec)
a.nam = as.data.frame(a.nam)
asv.names = paste0(a.nam$V1, a.nam$a.vec)
asv.names = as.data.frame(asv.names)
head(asv.names)
# apply ASV names to sequence table
taxa_names(eDNA1.16S) = asv.names$asv.names
taxa_names(eDNA1.16S)
eDNA1.16S

# write files for full dataset, no taxonomic filtering or sample filtering
write.csv(otu_table(eDNA1.16S), "asv_table_eDNA1_16S.csv")
write.csv(tax_table(eDNA1.16S), "tax_table_eDNA1_16S.csv")


###taxonomic filtering
rank_names(eDNA1.16S)
table(tax_table(eDNA1.16S)[, "Kingdom"], exclude = NULL)
eDNA1.16S
eDNA1.16SB <- subset_taxa(eDNA1.16S, Kingdom %in% "Bacteria")
eDNA1.16SB
table(tax_table(eDNA1.16SB)[, "Phylum"], exclude = NULL)
#remove Chloroplasts (at Order level) and Mitochondria (Family)
eDNA1.16S_nc <- subset_taxa(eDNA1.16SB, (Order!="Chloroplast") | is.na(Order))
eDNA1.16S_nc
eDNA1.16S_ncm <- subset_taxa(eDNA1.16S_nc, (Family!="Mitochondria") | is.na(Family))
eDNA1.16S_ncm
write.csv(otu_table(eDNA1.16S_ncm), file = "asv_table_eDNA1.16S_ncm.csv")
write.csv(tax_table(eDNA1.16S_ncm), file = "tax_table_eDNA1.16S_ncm.csv")
write.csv(sample_data(eDNA1.16S_ncm), file = "sample_data_eDNA1.16S_ncm.csv")
saveRDS(eDNA1.16S_ncm, "eDNA1.16S_ncm.rds")

#mesocosm samples
eDNA1.16S_noINV  = subset_samples(eDNA1.16S_ncm, treatment != "INV")
eDNA1.16S_noINV <- filter_taxa(eDNA1.16S_noINV, function (x) {sum(x > 0) > 0.5}, prune=TRUE)
eDNA1.16S_noINV
eDNA1.16S_noNEG  = subset_samples(eDNA1.16S_noINV, treatment != "NEG")
eDNA1.16S_noNEG <- filter_taxa(eDNA1.16S_noNEG, function (x) {sum(x > 0) > 0.5}, prune=TRUE)
eDNA1.16S_noNEG
write.csv(otu_table(eDNA1.16S_noNEG), file = "asv_table_eDNA1.16S_noNEG.csv")
write.csv(tax_table(eDNA1.16S_noNEG), file = "tax_table_eDNA1.16S_noNEG.csv")
write.csv(as.matrix(sample_data(eDNA1.16S_noNEG)), file = "sample_data_eDNA1.16S_noNEG.csv")
saveRDS(eDNA1.16S_noNEG, "eDNA1_16S_noNEG.rds")


#subset invasion control samples
eDNA1.16S_ncm.INV  = subset_samples(eDNA1.16S_ncm, treatment == "INV")
eDNA1.16S_ncm.INV 
eDNA1.16S_ncm.INV  <- filter_taxa(eDNA1.16S_ncm.INV, function (x) {sum(x > 0) > 0.5}, prune=TRUE)
eDNA1.16S_ncm.INV

#subset field controls (NEG)-
eDNA1.16S_ncm.NEG  = subset_samples(eDNA1.16S_ncm, treatment == "NEG")
eDNA1.16S_ncm.NEG 
eDNA1.16S_ncm.NEG  <- filter_taxa(eDNA1.16S_ncm.NEG, function (x) {sum(x > 0) > 0.5}, prune=TRUE)
eDNA1.16S_ncm.NEG
eDNA1.16S_ncm.NEG  = subset_samples(eDNA1.16S_ncm.NEG, sampling.week != "0")
eDNA1.16S_ncm.NEG <- filter_taxa(eDNA1.16S_ncm.NEG, function (x) {sum(x > 0) > 0.5}, prune=TRUE)
eDNA1.16S_ncm.NEG
write.csv(otu_table(eDNA1.16S_ncm.NEG), "asv_table_eDNA1_16S_ncm_NEG.csv")
write.csv(tax_table(eDNA1.16S_ncm.NEG), "tax_table_eDNA1_16S_ncm_NEG.csv")
write.csv(as.matrix(sample_data(eDNA1.16S_ncm.NEG)), "sample_data_eDNA1_16S_ncm_NEG.csv")
saveRDS(eDNA1.16S_ncm.NEG, "eDNA1_16S_ncm_NEG.rds")

#investigate negative controls
NEG_asvs <- read.csv("tax_table_eDNA1_16S_ncm_NEG.csv")
Mesocosm_asvs <- read.csv("tax_table_eDNA1_16S_noNEG.csv")
NEG_asvs_inMesos <- merge(NEG_asvs, Mesocosm_asvs, by = c("X","Kingdom","Phylum","Class","Order","Family","Genus","Species"))
# list of ASVs
asvs_inNegMesos <- NEG_asvs_inMesos$X 

# Prune the phyloseq object to remove these ASVs
eDNA1_16S_ncm_NEG_pruned <- prune_taxa(taxa_names(eDNA1.16S_ncm.NEG) %in% asvs_inNegMesos, eDNA1.16S_ncm.NEG)
eDNA1_16S_ncm_NEG_pruned
asvs_neg <- as.data.frame(otu_table(eDNA1_16S_ncm_NEG_pruned))
asvs_neg_all <- as.data.frame(otu_table(eDNA1.16S_ncm.NEG))
eDNA1_16S_ncm_NEG_pruned2 <- prune_taxa(taxa_names(eDNA1.16S_noNEG) %in% asvs_inNegMesos, eDNA1.16S_noNEG)
eDNA1_16S_ncm_NEG_pruned2
asvs_neg2 <- as.data.frame(otu_table(eDNA1_16S_ncm_NEG_pruned2))
asvs_all <- as.data.frame(otu_table(eDNA1.16S_noNEG))
asvs_all2 <- as.data.frame(otu_table(eDNA1.16S_noNEG_vs))

# calculate normalised read counts with DESeq2 - varianceStabilizingTransformation - fulldataset
dds = phyloseq_to_deseq2(eDNA1.16S_ncm, ~ 1)
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(dds), 1, gm_mean)
dds = estimateSizeFactors(dds, geoMeans = geoMeans)
sizeFactors(dds)
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vs_counts <- assay(vsd)
vs_counts[vs_counts<0]=0
taxa.eDNA1.16S_ncm1 <- read.csv("tax_table_eDNA1.16S_ncm.csv")
taxa.eDNA1.16S_ncm <- taxa.eDNA1.16S_ncm1[, -1]
row.names(taxa.eDNA1.16S_ncm) <- taxa.eDNA1.16S_ncm1$X #sample data row names must align with dada2 rowname outputs
taxa.eDNA1.16S_ncm = as.matrix(taxa.eDNA1.16S_ncm)
head(taxa.eDNA1.16S_ncm)
eDNA1.16S_ncm_vs = phyloseq(tax_table(taxa.eDNA1.16S_ncm), otu_table(vs_counts, taxa_are_rows = TRUE))
eDNA1.16S_ncm_vs
map.eDNA1.16S_ncm_vs = read.csv("sample_data_eDNA1.16S_ncm.csv")
head(map.eDNA1.16S_ncm_vs)
row.names(map.eDNA1.16S_ncm_vs)<- map.eDNA1.16S_ncm_vs$X #sample data row names must align with dada2 rowname outputs
head(map.eDNA1.16S_ncm_vs)
eDNA1.16S_ncm_vs = phyloseq(tax_table(taxa.eDNA1.16S_ncm), otu_table(vs_counts, taxa_are_rows = TRUE), 
                              sample_data(map.eDNA1.16S_ncm_vs))
eDNA1.16S_ncm_vs

#mesocosm samples
eDNA1.16S_noINV_vs  = subset_samples(eDNA1.16S_ncm_vs, treatment != "INV")
eDNA1.16S_noINV_vs <- filter_taxa(eDNA1.16S_noINV_vs, function (x) {sum(x > 0) > 0.5}, prune=TRUE)
eDNA1.16S_noINV_vs
eDNA1.16S_noNEG_vs  = subset_samples(eDNA1.16S_noINV_vs, treatment != "NEG")
eDNA1.16S_noNEG_vs <- filter_taxa(eDNA1.16S_noNEG_vs, function (x) {sum(x > 0) > 0.5}, prune=TRUE)
eDNA1.16S_noNEG_vs
write.csv(otu_table(eDNA1.16S_noNEG_vs), file = "asv_table_eDNA1_16S_noNEG_vs.csv")
write.csv(tax_table(eDNA1.16S_noNEG_vs), file = "tax_table_eDNA1_16S_noNEG_vs.csv")
write.csv(sample_data(eDNA1.16S_noNEG_vs), file = "sample_data_eDNA1_16S_noNEG_vs.csv")
saveRDS(eDNA1.16S_noNEG_vs, "eDNA1_16S_noNEG_vs.rds")