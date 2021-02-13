##############------------------------------------------------------------------------------------
#PACKAGES
##############------------------------------------------------------------------------------------
library(ggplot2)
library(velocyto.R)
library(loomR)
library(biomaRt)

##############------------------------------------------------------------------------------------
#PREPARATIONS
##############------------------------------------------------------------------------------------

#Set output directory
dir.create("/media/sf_aphid-probing-feeding_ms/analysis/009_velocyto/")
parDir <- ("/media/sf_aphid-probing-feeding_ms/analysis/009_velocyto/")
setwd(parDir)

#Load loom file
ldat.raw <- read.loom.matrices("/media/sf_aphid-probing-feeding_ms/data/EPG-prob_sq4_aphid_looms.loom")

#Set DE gene directory
DE.dir <- ("/media/sf_aphid-probing-feeding_ms/analysis/002_DESeq2/")

#Define time-points
t1 <- "0hr"
t2 <- "05hr"
t3 <- "1hr"
t4 <- "2hr"
t5 <- "3hr"
t6 <- "5hr"
t7 <- "7hr"
t8 <- "24hr"

times <- c(t1,t2,t3,t4,t5,t6,t7,t8)

#Prepare Biomart
mart = useMart(biomart="plants_mart", dataset="athaliana_eg_gene", host="plants.ensembl.org")

#Biomart function
BM_call <- function(genes) {
  data.frame(getBM(attributes = c("ensembl_gene_id",
                                  "tair_symbol",
                                  "external_gene_name",
                                  "name_1006"),
                   filters="ensembl_gene_id", values=as.character(unique(genes)), mart= mart))
}

#Velocyto plotting function for gene subsets
velocyto_geneset <- function(gene.subset,fit.quantile,filename){
  # exonic read (spliced) expression matrix
  emat <- ldat$spliced
  # intronic read (unspliced) expression matrix
  nmat <- ldat$unspliced
  # spanning read (intron+exon) expression matrix
  smat <- ldat$spanning
  # filter expression matrices based on some minimum max-cluster averages
  emat <- filter.genes.by.cluster.expression(emat,cell.colors,min.max.cluster.average = 5)
  nmat <- filter.genes.by.cluster.expression(nmat,cell.colors,min.max.cluster.average = 1)
  smat <- filter.genes.by.cluster.expression(smat,cell.colors,min.max.cluster.average = 0.5)
  #Create gene subsets
  emat <- as.matrix(subset(as.matrix(emat), as.character(tolower(row.names(emat))) %in% as.character(gene.subset)))
  nmat <- as.matrix(subset(as.matrix(nmat), as.character(tolower(row.names(nmat))) %in% as.character(gene.subset)))
  smat <- as.matrix(subset(as.matrix(smat), as.character(tolower(row.names(smat))) %in% as.character(gene.subset)))
  
  #Print length of gene set
  print(length(intersect(rownames(emat),rownames(nmat))))
  
  #Velocyto plotting
  rvel.qf <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells = 5,fit.quantile = fit.quantile)
  pdf(filename, height = 5, width = 5)
  pca.velocity.plot(rvel.qf,
                    nPcs=2,
                    plot.cols=1,
                    cell.colors=ac(cell.colors,alpha=1),
                    cex=3,pcount=0.1,pc.multipliers=c(1,-1),
                    arrow.scale = 1)
  dev.off()
}



#Suffix for filenames
suffix="aphid"

##############------------------------------------------------------------------------------------
#LOAD GENES
##############------------------------------------------------------------------------------------

#DESeq2 DEgenes
dfDE <- data.frame()
for(file in grep("_DESEq2_padj005.csv",list.files(DE.dir),value=T)){
  df <- read.csv(file.path(DE.dir,file))
  #Change column names
  colnames(df) <- c("gene",colnames(df)[2:ncol(df)])
  #Assign timepoint
  df$time <- gsub("05","0.5",gsub("hr","hap",unlist(strsplit(file, '_'))[1]))
  dfDE <- rbind(dfDE, df)
}

#Load all genes
dfAll <- data.frame()
for(file in grep("_DESEq2_no-padj-filter.csv",list.files(DE.dir),value=T)){
  df <- read.csv(file.path(DE.dir,file))
  #Change column names
  colnames(df) <- c("gene",colnames(df)[2:ncol(df)])
  #Assign timepoint
  df$time <- gsub("05","0.5",gsub("hr","hap",unlist(strsplit(file, '_'))[1]))
  dfAll <- rbind(dfAll, df)
}

#Load all Biomart genes for subsetting
BioMartGenes <- getBM(attributes = c("ensembl_gene_id"), mart= mart)$ensembl_gene_id

#Retrieve biomart information of all genes
df.bm.all <- BM_call(dfAll$gene)

#Retrieve biomart information of DE genes
df.bm.DE <- BM_call(dfDE$gene)

#Extract circadian genes from all detected genes
df.circadian.all <- df.bm.all[grep("circadian|clock", df.bm.all$name_1006),]
write.csv(df.circadian.all, "all.genes.circadian.clock.csv")

#Extract circadian genes from DE detected genes
df.circadian.DE <- df.bm.DE[grep("circadian|clock", df.bm.DE$name_1006),]
write.csv(df.circadian.DE, "DE.genes.circadian.clock.csv")

#Define subsets of genes to be processed for velocyto analysis

#All genes - ensembl gene id, external gene name and tair symbol for catching all genes from looms
gene.names.all <- c(tolower(df.bm.all$ensembl_gene_id),
                    tolower(df.bm.all$external_gene_name),
                    tolower(df.bm.all$tair_symbol))
#DE genes - ensembl gene id, external gene name and tair symbol for catching all genes from looms
gene.names.DE <- c(tolower(df.bm.DE$ensembl_gene_id),
                   tolower(df.bm.DE$external_gene_name),
                   tolower(df.bm.DE$tair_symbol))

#All genes circadian subset - ensembl gene id, external gene name and tair symbol for catching all genes from looms
gene.names.circadian.all <- c(tolower(df.circadian.all$ensembl_gene_id),
                              tolower(df.circadian.all$external_gene_name),
                              tolower(df.circadian.all$tair_symbol))

#DE genes circadian subset - ensembl gene id, external gene name and tair symbol for catching all genes from looms
gene.names.circadian.DE <- c(tolower(df.circadian.DE$ensembl_gene_id),
                             tolower(df.circadian.DE$external_gene_name),
                             tolower(df.circadian.DE$tair_symbol))

##############------------------------------------------------------------------------------------
#VELOCYTO PREPARATION
##############------------------------------------------------------------------------------------

#Modify loom filenames
str(ldat.raw)
ldat <- lapply(ldat.raw,function(x){
  colnames(x) <-  gsub(".aphid","",gsub("_",".",gsub(":Aligned.out.sorted.bam","",gsub("EPG-prob.","p.",colnames(x)))))
  x
})
str(ldat)

#Create named vector for colouring the PCA plot
samples <- c()
samples <- colnames(ldat$spliced)
cell.colors <- c(rep("#F0F8FF",4),
                 rep("#9f4646",4),
                 rep("#550c5c",4),
                 rep("#3d6838",4),
                 rep("#7c7c7c",4),
                 rep("#1be246",3),
                 rep("#e2df1b",4),
                 rep("#f94e04", 4)
names(cell.colors) <- samples

###################
#######VELOCYTO ALL GENES
###################

# exonic read (spliced) expression matrix
emat <- ldat$spliced
# intronic read (unspliced) expression matrix
nmat <- ldat$unspliced
# spanning read (intron+exon) expression matrix
smat <- ldat$spanning
# filter expression matrices based on some minimum max-cluster averages
emat <- filter.genes.by.cluster.expression(emat,cell.colors,min.max.cluster.average = 5)
nmat <- filter.genes.by.cluster.expression(nmat,cell.colors,min.max.cluster.average = 1)
smat <- filter.genes.by.cluster.expression(smat,cell.colors,min.max.cluster.average = 0.5)
# look at the resulting gene set
length(intersect(rownames(emat),rownames(nmat)))
fit.quantile <- 0.05
rvel.qf <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells = 5,fit.quantile = fit.quantile)

#Number of reads per gene
hist(log10(rowSums(ldat$spliced)+1),col='wheat',xlab='log10[ number of reads + 1]',main='number of reads per gene')

pdf(paste0("velocyto_PCA_all-genes_",suffix,".pdf"), height = 5, width = 5)
pca.velocity.plot(rvel.qf,
                  nPcs=2,
                  plot.cols=1,
                  cell.colors=ac(cell.colors,alpha=1),
                  cex=3,pcount=0.1,pc.multipliers=c(1,-1),
                  arrow.scale = 1)
dev.off()

###################
####### DE GENES
###################
velocyto_geneset(gene.subset=gene.names.DE,fit.quantile=0.05,filename=paste0("velocyto_PCA_",suffix,"_DE-genes-only.pdf"))

###################
####### ALL CIRCADIAN GENES
###################
velocyto_geneset(gene.subset=gene.names.circadian.all,fit.quantile=0.05,filename=paste0("velocyto_PCA_",suffix,"_all-circadian-genes.pdf"))

###################
####### ALL DE WITHOUT CIRCADIAN
###################
velocyto_geneset(gene.subset=as.character(setdiff(gene.names.DE, gene.names.circadian.all)),fit.quantile=0.05,filename=paste0("velocyto_PCA_",suffix,"_DE-genes-without-all-circadian-genes.pdf"))

###################
####### ALL GENES WITHOUT DE
###################
velocyto_geneset(gene.subset=setdiff(gene.names.all, gene.names.DE),fit.quantile=0.05,filename=paste0("velocyto_PCA_",suffix,"_all-genes-without-DE-genes.pdf"))

###################
####### ALL GENES WITHOUT CIRCADIAN
###################
velocyto_geneset(gene.subset=setdiff(gene.names.all, gene.names.circadian.all),fit.quantile=0.05,filename=paste0("velocyto_PCA_",suffix,"_all-genes-without-circadian-genes.pdf"))
