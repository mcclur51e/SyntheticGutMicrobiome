########## Call libraries for use ##########
# Packages from CRAN
library("ggplot2") # version 3.2.1

########## Call functions for use ##########
#source("~/Masters/taxa_summary.R",local=TRUE) # load fast_melt function
########## Define color palettes ##########
pal.pairBiome<-c("#1c541c","#2f8f2f","#49c349","#bfeabf",
                 "#B80000","#F00000","#FF7777","#ffcccc",
                 "#000080","#0000cd","#8282ff","#cfcfff",
                 "#623800","#c47000","#ff9914","#ffddb1",
                 "#430059","#7d00a7","#cc32ff","#eebbff",
                 "#626200","#c4c400","#ffff14","#ffffb1",
                 "#005f6c","#00a4bb","#1ee3ff","#bbf7ff",
                 "#750063","#c400a5","#ff13da","#ffb0f3",
                 "#1a1a1a","#808080","#d9d9d9","#ffffff","#808080")
pal.pairMini<-c(
  "#A6CEE3","#1F78B4",
  "#B2DF8A","#33A02C",
  "#FB9A99","#E31A1C",
  "#FDBF6F","#FF7F00",
  "#CAB2D6","#6A3D9A",
  "#FFFF99","#eded09",
  "#9aebff","#01cdff",
  "#e6e6e6","#404040"
)
pal.CB<-c("#e69f00","#56b4e9","#009e73","#f0e442","#0072b2","#D55e00","#cc79a7","#999999") # colorblind-friendly color palette
pal.CBmod<-c("#e69f00","#56b4e9","#009e73","#f0e442","#000000","#0072b2","#D55e00","#ac59ff","#cc79a7","#999999","#ffffff") # colorblind-friendly color palette
pal.SGM<-c("#5087E3","#8F0238","#F049CB","#FFFA00","#008136","#A7A7A7","#000000","#F2A122",
           "#7548B4","#55C41F","#5087E3","#CFC6F1","#e90909","#A16910","#00CFFF"
) # palette used for synthetic gut microbiome by Desai et al, 2016

load("output/physeq_initial.RData") # file produced from code_AFMprelim.R

########## Modify sample data ##########
sample_data(physeqR)$DateCollected <- as.Date(sample_data(physeqR)$DateCollected , format = "%m/%d/%y") # assign format to DateCollected to column
sample_data(physeqR)$genotype <- with(sample_data(physeqR),ifelse(Cage=="1","WT",
                                                              ifelse(Cage=="2","delta",
                                                                     "other"))) # define genotypes by Cage number
########## Modify taxonomy data ##########
df.tax <- as.data.frame(tax_table(physeqR)) # create data frame from taxonomy tale
# assign some taxonomy that database didn't do well on and create G.epithet notation for the rest
df.tax$Spp <- with(df.tax, ifelse(Sequence=="GCAGTGAGGAATATTGGTCAATGGGCGCAGGCCTGAACCAGCCAAGTAGCGTGAAGGATGACTGCCCTATGGGTTGTAAACTTCTTTTATATGGGAATAAAGTTTTCCACGTGTGGAATTTTGTATGTACCATATGAATAAGGATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATCCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGAGCGTAGGTGGACAGTTAAGTCAGTTGTGAAAGTTTGCGGCTCAACCGTAAAATTGCAGTTGATACTGGCTGTCTTGAGTACAGTAGAGGTGGGCGGAATTCGTGGTGTAGCGGTGAAATGCTTAGATATCACGAAGAACTCCGATTGCGAAGGCAGCTCACTGGACTGCAACTGACACTGATGCTCGAAAGTGTGGGTATCAAACAGGATTAGA","B.thetaiotaomicron",
                                  ifelse(Sequence=="GCAGTGAGGAATATTGGTCAATGGACGCGAGTCTGAACCAGCCAAGTAGCGTGAAGGATGACTGCCCTATGGGTTGTAAACTTCTTTTATATGGGAATAAAGTGGTCCACGTGTGGACTTTTGTATGTACCATATGAATAAGGATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATCCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGAGCGTAGGCGGATTGTTAAGTCAGTTGTGAAAGTTTGCGGCTCAACCGTAAAATTGCAGTTGATACTGGCAGTCTTGAGTGCAGTAGAGGTGGGCGGAATTCGTGGTGTAGCGGTGAAATGCTTAGATATCACGAAGAACTCCGATTGCGAAGGCAGCTCACTGGAGTGTAACTGACGCTGATGCTCGAAAGTGTGGGTATCAAACAGGATTAGA","B.caccae",
                                         ifelse(Sequence=="GCAGTGAGGAATATTGGTCAATGGACGGGAGTCTGAACCAGCCAAGTAGCGTGAAGGATGACTGCCCTATGGGTTGTAAACTTCTTTTATATGGGAATAAAGTGATCCACGTGTGGATTTTTGTATGTACCATATGAATAAGGATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATCCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGAGCGTAGGTGGACAGTTAAGTCAGTTGTGAAAGTTTGCGGCTCAACCGTAAAATTGCAGTTGATACTGGCTGTCTTGAGTACAGTAGAGGTGGGCGGAATTCGTGGTGTAGCGGTGAAATGCTTAGATATCACGAAGAACTCCGATTGCGAAGGCAGCTCACTGGACTGCAACTGACACTGATGCTCGAAAGTGTGGGTATCAAACAGGATTAGA","B.faecis",
                                                ifelse(is.na(Genus),Family,
                                                       ifelse(is.na(Species),Genus,paste0(substr(df.tax[,"Genus"],1,1),".",Species)))))))
df.tax$Spp2 <- with(df.tax, ifelse(Genus=="Bacteroides",Spp,Genus)) # create new column for species designations that will be used for plots
physeqR <- phyloseq(otu_table(physeqR),tax_table(as.matrix(df.tax)),sample_data(physeqR)) # return modified taxonomy table to phyloseq object

### Subset to only animal reads 
phyR.animal <- subset_taxa(physeqR, Kingdom=="Animalia") # subset data to only animal reads
write.csv(tax_table(phyR.animal),"output/table_taxaAnimal.csv") # write taxonomy table of only animal reads


### Continue with only bacteria reads 
phyR.bacteria <- subset_taxa(physeqR,Kingdom=="Bacteria" & taxa_sums(physeqR) > 1) # subset data to only bacteria reads
write.csv(tax_table(phyR.bacteria),"output/table_taxaBacteria.csv") # write taxonomy table of only bacteria reads

phyR.control <- subset_samples(phyR.bacteria, Control%in%c("negExtrac","posExtract")) # subset to only controls

### Sample pre-processing/trimming
phyR.samples <- subset_samples(phyR.bacteria, Control%in%c("sample") & sample_sums(phyR.bacteria)>1) # subset to only samples 
phyS.samples <- subset_taxa(phyR.samples, taxa_sums(phyR.samples)>1) # remove empty taxa
phyN.samples <- subset_samples(phyS.samples, sample_sums(phyS.samples)>1) # remove empty samples
phyT.samples <- transform_sample_counts(phyN.samples, function(x){x / sum(x)}) # convert raw read counts to abundance percentages (phyloseq Transformed . controls)

#########################################
########## Rearrange data ##########
#########################################
### Pre-process data
ls.SGM <- c("Akkermansia","Collinsella","Desulfovibrio","Escherichia/Shigella","Faecalibacterium",
            "Roseburia","Marvinbryantia","Cenarchaeum","Eubacterium","Agathobacter",
            "Barnesiella","Bacteroides","Clostridium","Lachnoclostridium") # list of taxa to plot from synthetic gut microbiome by Desai et al, 2016
phyT.bacteroidesSpp <- subset_taxa(phyT.samples, Genus%in%ls.SGM) # subset to taxa in synthetic gut microbiome by Desai et al, 2016
phyT.bacteroidesSpp <- subset_taxa(phyT.bacteroidesSpp, !Spp%in%c("B.faecis","B.fragilis","Bacteroides","Eubacterium")) # additional subset of taxa 
phyG.bacteroidesSpp <- tax_glom(phyT.bacteroidesSpp,"Spp2", NArm=F) # combine taxa by Spp2
TAXsm <- as.data.frame(tax_table(phyG.bacteroidesSpp))[,c("Kingdom","Phylum","Class","Order","Family","Genus","Spp2")] # remove columns from taxonomy table

### Add column to taxonomy table for hex code Desai colors
TAXsm$Desai <- with(TAXsm, ifelse(Spp2=="Akkermansia","#8F0238",
                                  ifelse(Spp2=="Collinsella","#F2A122",
                                         ifelse(Spp2=="Desulfovibrio","#7548B4",
                                                ifelse(Spp2=="Escherichia/Shigella","#55C41F",
                                                       ifelse(Spp2=="Faecalibacterium","#CFC6F1",
                                                              ifelse(Spp2=="Roseburia","#00CFFF",
                                                                     ifelse(Spp2=="Marvinbryantia","#A16910",
                                                                            ifelse(Spp2=="Cenarchaeum","#e90909",
                                                                                   ifelse(Spp2=="Lachnoclostridium","#e90908",
                                                                                          ifelse(Spp2=="Clostridium","#e90909",
                                                                                                 ifelse(Spp2=="Eubacterium","#5087E3",
                                                                                                        ifelse(Spp2=="Barnesiella","#000000",
                                                                                                               ifelse(Spp2=="B.caccae","#F049CB",
                                                                                                                      ifelse(Spp2=="B.uniformis","#A7A7A7",
                                                                                                                             ifelse(Spp2=="B.ovatus","#FFFA00",
                                                                                                                                    ifelse(Spp2=="B.thetaiotaomicron","#008136",
                                                                                                                                           "#FFFFFF")))))))))))))))))
phyG.bacteroidesSpp <- phyloseq(otu_table(phyG.bacteroidesSpp),tax_table(as.matrix(TAXsm)),sample_data(phyG.bacteroidesSpp)) # return modified taxonomy table to phyloseq object 

#########################################
########## Barplot of SGM taxa ##########
#########################################

### Create barplot with Desai et al, 2016 color scheme
phy.plot <- phyG.bacteroidesSpp # assign data to plot
phyT.plot <- transform_sample_counts(phy.plot, function(x){x / sum(x)}) # convert raw read counts to abundance percentages (phyloseq Transformed . controls)
ps.plot <- psmelt(phyT.plot)
pBar.bacteroidesSpp <- ggplot(ps.plot, aes(x=Mouse,y=Abundance,fill=Spp2)) +
  geom_bar(position="stack", stat="identity",width=1) +
  ylim(0,1) +
  facet_grid(genotype+ExperimentReplicate~DateCollected, scales="free_x", space="free_x") +
  theme_bw() +
  theme(text=element_text(size=12),axis.text.x = element_text(angle=90),panel.spacing = unit(0, "lines")) +
  scale_fill_manual(values = unique(ps.plot$Desai[order(ps.plot$Spp2)]))
pBar.bacteroidesSpp
ggsave(pBar.bacteroidesSpp,filename="plots/bar_bacteroidesSpp.pdf", dpi="retina",width=20,height=10,units="in") # Save figure to .pdf file

#########################################
########## NMDS ##########
#########################################
phyT.SGM <- subset_taxa(phyG.bacteroidesSpp, Genus%in%ls.SGM) # subset to taxa in synthetic gut microbiome by Desai et al, 2016

ord.samples <- ordinate(phyT.SGM, "NMDS", "bray") # calculate ordination
### Plot using sample data
pNMDS.samples <- plot_ordination(phyT.SGM, ord.samples, type="sample", color="genotype") +
  stat_ellipse() +
  geom_point(aes(color=genotype,shape=ExperimentReplicate,size=2)) +
  theme_bw() +
  theme(text=element_text(size=12),axis.text.x = element_text(angle=90),panel.spacing = unit(0, "lines")) +
  scale_color_manual(values = pal.CB)
pNMDS.samples
ggsave(pNMDS.samples,filename="plots/NMDS_genotype.pdf", dpi="retina",width=8,height=8,units="in") # Save figure to .pdf file

#plot using taxonomy data
pal.desai <- as.data.frame(unique(tax_table(phyT.SGM)[,"Desai"][order(tax_table(phyT.SGM)[,"Spp2"])]))$Desai # palette to color taxonomy like Desai et al, 2016
pNMDS.samplesTax <- plot_ordination(phyT.SGM, ord.samples, type="taxa", color="Spp2") +
  geom_point(aes(color=Spp2,size=2)) +
  facet_wrap(~Phylum,3) +
  theme_bw() +
  theme(text=element_text(size=12),axis.text.x = element_text(angle=90),panel.spacing = unit(0, "lines")) +
  scale_color_manual(values=pal.desai)
pNMDS.samplesTax
ggsave(pNMDS.samplesTax,filename="plots/NMDS_taxa.pdf", dpi="retina",width=10,height=10,units="in") # Save figure to .pdf file












#################################
########## Stream plot ##########
#################################

########## THIS CODE DOES NOT YET WORK ##########
library("ggstream")
phy.stream <- phyG.bacteroidesSpp
sample_data(phy.stream) <- sample_data(phy.stream)[,c("genotype","ExperimentReplicate","DateCollected")] # reduce sample data to only necessary columns
sample_data(phy.stream)$stream <- paste0(sample_data(phy.stream)$genotype,sample_data(phy.stream)$ExperimentReplicate,sample_data(phy.stream)$DateCollected) # create groups for averaging data
phyMS.stream <- merge_samples(phy.stream, "stream",fun=mean) # merge samples by group

### Recreate columns that were affected by merge ###
# recreate genotype column
sample_data(phyMS.stream)$genotype <- with(sample_data(phyMS.stream), ifelse(substr(rownames(sample_data(phyMS.stream)),1,2)=="WT",
                                                                           "WT",substr(rownames(sample_data(phyMS.stream)),1,5))) 
# recreate ExperimentalReplicate column
sample_data(phyMS.stream)$ExperimentReplicate <- with(sample_data(phyMS.stream), ifelse(substr(rownames(sample_data(phyMS.stream)),1,2)=="WT",
                                                                                      substr(rownames(sample_data(phyMS.stream)),3,6),
                                                                                      substr(rownames(sample_data(phyMS.stream)),6,9)))
sample_data(phyMS.stream)$DateCollected <- as.Date(sample_data(phyMS.stream)$DateCollected,origin="1970-01-01") # return dates from numerical data


phyG.stream <- tax_glom(phyMS.stream,"Desai",NArm=F)
phyG.stream <- subset_taxa(phyG.stream,taxa_sums(phyG.stream)>0.001)

#phyMS.stream2 <- subset_taxa(phyMS.stream,taxa_sums(phyMS.stream)>0.001)
phyT.stream <- transform_sample_counts(phyG.stream, function(x){x / sum(x)}) # convert raw read counts to abundance percentages (phyloseq Transformed . controls)

ps.stream <- psmelt(subset_samples(phyT.stream,genotype=="WT"))

ps.stream <- ps.stream %>% mutate(Spp2=factor(Spp2, levels=c("Akkermansia", "Collinsella", "Desulfovibrio", 
                                                             "Escherichia/Shigella", "Faecalibacterium", "Roseoburia", 
                                                             "Marvinbryantia", "Lachnoclostridium","Agathobacter",
                                                             "Barnesiella","B.caccae","B.uniformis","B.ovatus",
                                                             "B.thetaiotaomicron")))

ps.stream$DateCollected <- as.character(ps.stream$DateCollected)

ps.temp <- subset(ps.stream, DateCollected < "2023-04-03")

ps.a <- ps.stream[1:50,]
pStream <- ggplot(ps.a, aes(x=DateCollected,y=Abundance,fill=Spp2)) +
  geom_stream(type="proportional") +
  #scale_x_date(date_breaks="1 week") +
  theme_bw() +
  theme(text=element_text(size=12),axis.text.x = element_text(angle=90),panel.spacing = unit(0, "lines")) +
  scale_fill_manual(values=pal.SGM2) 
pStream
ggsave(pStream,filename="plots/stream_grant.pdf", dpi="retina",width=10,height=10,units="in") # Save figure to .pdf file



#################################
########## END ##################
#################################


























