#####################################################################
########## Processing Set −up #######################################
#####################################################################
#Rpath <- "~/R/x86_64-redhat-linux-gnu-library/3.6/"
########## Call libraries for use ##########
# Packages from CRAN
library("ggplot2") # version 3.2.1
#library("reshape2") # version 1.4.3
library("data.table") # version 1.12.6
#library("RColorBrewer") # design new color palette for data # version 1.1-2
#library("scales") # for scientific notation in plots # version 1.0.0
#library("cowplot") # used to prepare final figuresMd # version 1.0.0
#library("tidyr") # version 1.0.0
library("dplyr") # version 0.8.3
#library("vegan") # version 2.5-6

# Packages from Bioconductor
library("phyloseq") # version 1.28.0
library("decontam") # identify contaminant ASVs # version 1.4.0
#library("DESeq2") # The DESeq2 model internally corrects for library size, so transformed or normalized values such as counts scaled by library size should not be used as input. # version 1.24.0
library("dada2") # version 1.12.1
#library("Biostrings") # 2.52.0
#library("microbiome") # version 1.6.0

##################################################################################
########## If running from scrum (or similar), use the following block ###########
##################################################################################
#input <- commandArgs(trailingOnly = TRUE) # Read arguments when beginning script in scrum
#path <- as.character(input[1])
#setwd(path) # set working directory
#dir.create("output") # create directory for output files to go

################################################################################
########## If running by line (i.e. RStudio) use the following block ###########
################################################################################
path <- "~/Desktop/Wannigan_CC/" # change this to directory where you will be working
setwd(path) # set working directory
dir.create("output") # create directory for output files to go
dir.create("plots") # create directory for plots to go

#####################################################################
##### ASV assignment with DADA2 for V3-4 samples #####################################
#####################################################################
fnFs <- sort(list.files(paste0(path,"raw/"), pattern="R1_001", full.names = TRUE))
fnRs <- sort(list.files(paste0(path,"raw/"), pattern="R2_001", full.names = TRUE))

if(length(fnFs) != length(fnRs)) stop("At least one sample is unpaired. Please check forward and reverse reads are present for all samples")
if(length(fnFs)==0) stop("No forward reads. Do you need to modify directory location?")
if(length(fnRs)==0) stop("No reverse reads. Do you need to modify directory location?")

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
#sample.namesR <- sapply(strsplit(basename(fnRs), "_"), `[`, 1) # to check if reverse reads have same sample names as forward
#write.csv(sample.names,"output/sampleNames.csv") # print list of sample names to check if they match map file
table_map <- data.frame(read.csv("raw/table_map.csv", header = TRUE, row.names = 1, check.names=FALSE)) # reads csv file into data.frame with row names in column 1

save(fnFs,file=("output/output_fnFs.RData")) # Save in a .RData file 
save(fnRs,file=("output/output_fnRs.RData")) # Save in a .RData file 

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
save(filtFs,file=("output/output_filtFs.RData")) # Save in a .RData file
save(filtRs,file=("output/output_filtRs.RData")) # Save in a .RData file

filtered <- filterAndTrim(fnFs, filtFs, 
                          fnRs, filtRs,
                          trimLeft=13, trimRight=0,
                          maxLen=500, minLen = 200,
                          maxN=1, maxEE=c(2,5), truncQ=2, rm.phix=TRUE,
                          compress=TRUE, verbose=TRUE, multithread=TRUE)
save(filtered,file=("output/output_filtered.RData")) # Save .RData file

### Subset filtFs and filtRs to include files with > 0 reads (some files may have been emptied during error analysis) 
filtFs <- filtFs[file.exists(filtFs)]
filtRs <- filtRs[file.exists(filtRs)]
save(filtFs,file=("output/output_filtFs.RData")) # Save .RData file
save(filtRs,file=("output/output_filtRs.RData")) # Save .RData file 

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
save(errF,file=("output/output_errF.RData")) # Save .RData file
save(errR,file=("output/output_errR.RData")) # Save .RData file
#plotErrors(errF, nominalQ=TRUE)


dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
save(mergers,file=("output/output_mergers.RData")) # Save in a .RData file

seqtab <- makeSequenceTable(mergers)
# table(nchar(getSequences(seqtab))) # use to check distribution of sequence lengths
seqtab.all <- seqtab[,nchar(colnames(seqtab)) %in% 410:440]
seqtab.noBim <- removeBimeraDenovo(seqtab.all, method="consensus", multithread=TRUE)

#####################################################################
########## Transfer data into Phyloseq ##############################
##### This preliminary Phyloseq object will be used to trim #########
##### singletons from the dataset before assigning taxonomy #########
#####################################################################


OTUraw = otu_table(seqtab.noBim, taxa_are_rows=FALSE) # assigns ASV table from DADA output
MAP = sample_data(table_map) # reads csv file into data.frame with row names in column 1
phyW = phyloseq(OTUraw, MAP) # prepare phyloseq object
save(phyW, file=("output/phyW.RData"))

### Dataset too small (too few samples), so no trimming of ASVs in only 1 sample

taxRaw <- assignTaxonomy(OTUraw, "~/Masters/silva_nr_v132_train_set_RossMod.fa", multithread=TRUE) # assign taxonomy to trimmed dataset
save(taxRaw,file=("output/output_taxRaw.RData")) # Save the phyloseq data object in a .RData file 

#####################################################################
########## Transfer data into Phyloseq ##############################
#####################################################################
#load("output/output_taxRaw.RData")
#table_map <- data.frame(read.csv("raw/table_map.csv", header = TRUE, row.names = 1, check.names=FALSE)) # reads csv file into data.frame with row names in column 1

OTU = otu_table(OTUraw, taxa_are_rows=FALSE) # assigns ASV table from trimmed DADA output
TAX = tax_table(taxRaw) # assigns taxonomy table
MAP = sample_data(table_map) # reads csv file into data.frame with row names in column 1
physeq = phyloseq(OTU, TAX, MAP) # prepare phyloseq object

### create and assign ASV numbers to consensus sequences for easier reference ###
ASV <- paste0("ASV", seq(ntaxa(physeq))) 
Sequence <- row.names(tax_table(physeq)) 
bind.asv <- cbind(tax_table(physeq),ASV)
bind.seq <- cbind(bind.asv,Sequence) 
TAX2 = tax_table(as.matrix(bind.seq)) # define new taxonomy table
colnames(TAX2) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species","Strain","ASV","Sequence")
physeqR = phyloseq(OTU, TAX2, MAP) # prepare phyloseq object modified with ASV numbers

df.map <- as.data.frame(sample_data(physeqR))
df.map$LibrarySize <- sample_sums(physeqR)
df.map$AnimalReads <- sample_sums(subset_taxa(physeqR, Kingdom=="Animalia"))
df.map$BacteriaReads <- sample_sums(subset_taxa(physeqR, Kingdom=="Bacteria"))
df.map$RatioReads <- df.map$BacteriaReads / df.map$AnimalReads
df.map$BacteriaPercent <- df.map$BacteriaReads / df.map$LibrarySize

physeqR = phyloseq(OTU,TAX2,sample_data(df.map))
save(physeqR,file=("output/physeq_initial.RData")) # Save the phyloseq data object in a .RData file 
write.csv(tax_table(physeqR),"output/table_tax.csv") # Save taxonomy table as .csv
write.csv(otu_table(physeqR),"output/table_otu.csv") # Save ASV table as .csv
#load("output/physeq_initial.RData")

###################################################################
########## You have now successfully imported all #################
########## the raw data into R and performed preliminary ##########
########## clean-up so that you are ready to proceed ##############
########## with further analysis ##################################
###################################################################