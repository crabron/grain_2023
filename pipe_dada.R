require(dada2)
require(Biostrings)
require(DECIPHER)
require(phyloseq)
require(seqinr)
require(data.table)
require(metagMisc)
require(tibble)
library(phyloseq)
library(tidyverse)
library(vegan)

mult = 10
mlt = 10
path_trein_set = "~/storage/somebases/silva_nr99_v138.1_train_set.fa.gz"
path_trein_set_species = "~/storage/somebases/silva_species_assignment_v138.1.fa.gz"
name_Brief = "nameBrief.txt"
path <- "/home/gladkov/storage/rnf_2023/raw"
active_dir <- "/home/gladkov/storage/june_2023/"
fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
on.exit()
filtFs <- file.path(active_dir, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(active_dir, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(200,180), trimLeft=c(19,20), maxN=0, truncQ=2, maxEE=c(2,5), rm.phix=TRUE, compress=TRUE, multithread=mult)
errF <- learnErrors(filtFs, multithread=mult)
errR <- learnErrors(filtRs, multithread=mult)
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names
dadaFs <- dada(derepFs, err=errF, multithread=mult)
dadaRs <- dada(derepRs, err=errR, multithread=mult)
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))
getN <- function(x) sum(getUniques(x))
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=mult, verbose=TRUE)

briefToSeq <- colnames(seqtab.nochim)
names(briefToSeq) <- paste0("Seq", seq(ncol(seqtab.nochim))) 
st.brief <- seqtab.nochim
colnames(st.brief) <- names(briefToSeq)

write.table(briefToSeq, file = name_Brief, sep = "\t")


taxa.dada2 <- assignTaxonomy(briefToSeq, path_trein_set , multithread=mult)
taxa.dada2.species <- assignSpecies(briefToSeq, path_trein_set_species) # maybe use not seqs but brief 
rownames(taxa.dada2.species) <- rownames(briefToSeq)
briefToSeq.df <- data.frame(briefToSeq)
rownames(taxa.dada2.species) <- rownames(briefToSeq.df)
rownames(taxa.dada2) <- rownames(taxa.dada2.species)
taxa <- cbind2(taxa.dada2, taxa.dada2.species[,2])
colnames(taxa)[7] <- "Species"

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
write.table(track, file = "track.tsv", sep= "\t", col.names = NA, quote=FALSE)

st.brief.t.df <- data.frame(t(st.brief))
write.table(st.brief.t.df, file = "otu_table.txt", sep= "\t", col.names = NA, quote=FALSE)

briefToSeq.ls <- as.list(briefToSeq.df[,c("briefToSeq")])
briefToSeq.names <- as.list(rownames(briefToSeq.df))
seqinr::write.fasta(briefToSeq.ls, briefToSeq.names , "rep_seq.fasta", as.string = FALSE)

write.table(taxa, file = "taxa.txt", sep= "\t", col.names = NA, quote=FALSE)
