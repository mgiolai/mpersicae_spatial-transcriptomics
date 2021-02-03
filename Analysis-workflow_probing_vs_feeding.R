######################
## LOAD PACKAGES
######################
library(DESeq2)
library(tidyverse)
library(ggsci)
library(reshape2)
library(biomaRt)
library(clusterProfiler)
library(GOSemSim)
library(org.At.tair.db)
library(apcluster)
library(treemap)
library(factoextra)
library(reticulate)
library(rvest)
library(stringr)
library(STRINGdb)
library(igraph)
library(ggrepel)
library(dichromat)
library(RColorBrewer)

##############------------------------------------------------------------------------------------
#DIRECTORIES
##############------------------------------------------------------------------------------------

#Working directory
parDir <- ("C:/Users/michael/Desktop/_SCRIPTS_AND_PROJECTS/aphid-probing-feeding_ms/analysis")
dir.create(file.path(parDir))

#HTSeq-count files
htseqDirProbing <- ("C:/Users/michael/Desktop/_SCRIPTS_AND_PROJECTS/aphid-probing-feeding_ms/data/probing/HTSeq-count/")
htseqDirFeeding <- ("C:/Users/michael/Desktop/_SCRIPTS_AND_PROJECTS/aphid-probing-feeding_ms/data/feeding/HTSeq-count/")

##############------------------------------------------------------------------------------------
#PREPARE PLANT BIOMART
##############------------------------------------------------------------------------------------

mart = useMart(biomart="plants_mart", dataset="athaliana_eg_gene", host="plants.ensembl.org")

##############------------------------------------------------------------------------------------
#DEFINITIONS AND FUNCTIONS
##############------------------------------------------------------------------------------------

#Timepoints and times
t1 <- "0hr"
t2 <- "05hr"
t3 <- "1hr"
t4 <- "2hr"
t5 <- "3hr"
t6 <- "5hr"
t7 <- "7hr"
t8 <- "24hr"
times <- c(t1,t2,t3,t4,t5,t6,t7,t8)

#Squares
sq1 <- "sq1"
sq2 <- "sq2"
sq3 <- "sq3"
sq4 <- "sq4"
sq5 <- "sq5"
sq6 <- "sq6"
sq7 <- "sq7"
sq8 <- "sq8"
squares <- c(sq1,sq2,sq3,sq4,sq5,sq6,sq7,sq8)

#GO term enrichment function
GO_enrich <- function(filename, genes, type){
  dfego <- enrichGO(gene = unique(genes),OrgDb="org.At.tair.db",ont=type,keyType="TAIR",pAdjustMethod="BH",readable=FALSE,pvalueCutoff=0.05,qvalueCutoff=0.05)
  if(nrow(dfego@result) > 0){
    dfego <- dfego@result[which(dfego@result$qvalue < 0.05),]
    if(nrow(dfego) > 0){
      dfego$GeneRatio <- gsub("/","|",dfego$GeneRatio, fixed=TRUE)
      dfego$BgRatio <- gsub("/","|",dfego$BgRatio, fixed=TRUE)
      dfego$geneID <- gsub("/","|",dfego$geneID, fixed=TRUE)
      write.csv(dfego, file = paste0(filename,".csv"))
    }
  }
}

#Save PCA plots function
plotPCA <- function(pca, filename){
  pdf(paste0(filename,".pdf"),height=5,width=6)
  for(dimA in seq(1,10)){
    for(dimB in seq(1,10)){
      df <- pca$x %>%
        as.data.frame() %>%
        dplyr::select(paste0("PC",dimA), paste0("PC",dimB)) %>%
        mutate(index=row.names(.)) %>%
        mutate(treatment=ifelse(grepl("hap",index),"probing","feeding")) %>%
        mutate(replicate=sapply(strsplit(index, '.',fixed=T), '[[', 2)) %>%
        mutate(Timepoint=gsub("05","0.5",sapply(strsplit(index, '.',fixed=T), '[[', 1))) %>%
        mutate(index=gsub("05","0.5",index))
      
      p <- ggplot(df, aes(x=pull(df, paste0("PC",dimA)), y=pull(df, paste0("PC",dimB)))) +
        geom_point(aes(color=Timepoint),size=4.5) +
        geom_text_repel(aes(color=Timepoint, label=index), size=6) +
        scale_color_nejm() +
        theme_bw(base_size = 14) +
        theme(plot.title = element_text(hjust = 0.5, size=20),
              axis.text=element_text(size=20),
              axis.title=element_text(size=20),
              axis.text.x=element_text(angle=60,hjust = 1.0),
              legend.position = "right",
              legend.text=element_text(size=20),
              legend.title=element_text(size=20)) +
        xlab(paste0("PC",dimA," ",round(get_eig(res.pca)[dimA,2],2),"%")) +
        ylab(paste0("PC",dimB," ",round(get_eig(res.pca)[dimB,2],2),"%"))
      plot(p)
    }
  }
  dev.off()
}

##############------------------------------------------------------------------------------------
#COUNT READS
##############------------------------------------------------------------------------------------

dir.create(file.path(parDir,"001_read-counts"))
setwd(file.path(parDir,"001_read-counts"))

#Count number of mapped reads, exclude files with too few mapped reads 
dfMapF <- data.frame()
for(file in list.files(htseqDirFeeding)){
  dfMapF <- rbind(dfMapF,
                 data.frame("file"=file,
                            "reads"=sum(head(read.table(file.path(htseqDirFeeding,file), header = F), 33977)$V2)))
}
dfMapF[which(dfMapF$reads < 100000),]$file
write.csv(dfMapF,"htseq.readcounts.feeding.csv")

#Count number of mapped reads, exclude files with too few mapped reads 
dfMapP <- data.frame()
for(file in list.files(htseqDirProbing)){
  dfMapP <- rbind(dfMapP,
                 data.frame("file"=file,
                            "reads"=sum(head(read.table(file.path(htseqDirProbing,file), header = F), 33977)$V2)))
}
dfMapP[which(dfMapP$reads < 100000),]$file
write.csv(dfMapP,"htseq.readcounts.probing.csv")

#Determine samples with less than 'n' reads
samplesExclude <- dfMapP[which(dfMapP$reads < 100000),]$file

#Load files and exclude samples with less than n reads
probingFiles <- setdiff(list.files(htseqDirProbing), samplesExclude)

#Number of detected genes with >= 1 read
dfngenes <- data.frame()
for(file in list.files(htseqDirFeeding)){
  df <- head(read.table(file.path(htseqDirFeeding,file), header = F), 33977)
  dfngenes <- rbind(dfngenes, data.frame(library=file,"genes"=length(df[df$V2 >= 1,]$V1)))
}
for(file in list.files(htseqDirProbing)){
  df <- head(read.table(file.path(htseqDirProbing,file), header = F), 33977)
  dfngenes <- rbind(dfngenes, data.frame(library=file,"genes"=length(df[df$V2 >= 1,]$V1)))
}
write.csv(dfngenes,"found_genes_per_library.csv")
##############------------------------------------------------------------------------------------
#BUILD DESEQ TABLE FROM READ FILTERED FILENAMES
##############------------------------------------------------------------------------------------

dir.create(file.path(parDir,"002_DESeq2"))
setwd(file.path(parDir,"002_DESeq2"))

#Build DEseq2 table
sampleTable <- data.frame("replicateName"=gsub("htseq-count_","",gsub(".txt","",probingFiles)),
                          "filename"=probingFiles,
                          "replicate"=sapply(strsplit(probingFiles, '_'), '[[', 5),
                          "treatment"=gsub(".txt","",sapply(strsplit(probingFiles, '_'), '[[', 6)))
sampleTable$treatment <- as.factor(sampleTable$treatment)
write.csv(sampleTable,"sampleTable.csv")

for(time in times){
  #Build DESeq2Table
  st <- sampleTable[grep(paste0("_",time,"_"),sampleTable$filename),]
  st <- st[grep("sq4",st$filename),]
  DESeq2Table<-DESeqDataSetFromHTSeqCount(sampleTable=st, directory=htseqDirProbing, design=~treatment)
  #Subset for counts
  DESeq2Table <- DESeq2Table[rowSums(counts(DESeq2Table))>10,]
  #Estimate size factors
  DESeq2Table <- estimateSizeFactors(DESeq2Table)
  #Perform PCA analysis
  res.pca <- prcomp(as.data.frame(t(counts(DESeq2Table, normalized=TRUE)))+1,  scale=T)

  p <- res.pca$x %>%
    as.data.frame() %>%
    mutate(index=row.names(.)) %>%
    mutate(treatment=sapply(strsplit(index, '_'), '[[', 5)) %>%
    mutate(replicate=sapply(strsplit(index, '_'), '[[', 4)) %>% 
    ggplot(aes(x=PC1, y=PC2, color=treatment, label=replicate)) +
    geom_point(size=2) +
    geom_text() +
    scale_color_nejm() +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, size=20),
          axis.text=element_text(size=20),
          axis.title=element_text(size=20),
          axis.text.x=element_text(angle=60,hjust = 1.0),
          legend.position = "right") +
    ggtitle(paste0(time," hap")) +
    xlab(paste0("PC1 ",round(get_eig(res.pca)[1,2],2),"%")) +
    ylab(paste0("PC2 ",round(get_eig(res.pca)[2,2],2),"%"))
  plot(p)
  
  #DE gene calling
  results <- as.data.frame(results(DESeq(DESeq2Table), alpha=0.05, contrast=c("treatment", "aphid", "ctrl")))
  write.csv(results,paste0(time,"_DESEq2_no-padj-filter.csv"))
  results <- results[which(results$padj < 0.05),]
  #Export DE gene tables
  if(nrow(results) > 0){
    write.csv(results,paste0(time,"_DESEq2_padj005.csv"))
  }
}

##############------------------------------------------------------------------------------------
#EXPORT A TAIR ANNOTATED LIST OF ALL DE GENES WITH TIME OF DIFFERENTIAL EXPRESSION
##############------------------------------------------------------------------------------------

dfDE <- data.frame()
for(file in grep("_DESEq2_padj005.csv",list.files(),value=T)){
  df <- read.csv(file)
  #Change column names
  colnames(df) <- c("gene",colnames(df)[2:ncol(df)])
  #Assign timepoint
  df$time <- gsub("05","0.5",gsub("hr","hap",unlist(strsplit(file, '_'))[1]))
  dfDE <- rbind(dfDE, df)
}

#Merge with TAIR name
df <- getBM(attributes = c("ensembl_gene_id","tair_symbol"), filters="ensembl_gene_id", values=unique(dfDE$gene), mart= mart)
dfDE$TAIR <- unlist(lapply(dfDE$gene, function(x) unique(ifelse(df[df$ensembl_gene_id==x,]$tair_symbol=="",x,df[df$ensembl_gene_id==x,]$tair_symbol[1]))))

write.csv(dfDE %>% arrange(factor(time, levels=gsub("05","0.5",gsub("hr","hap",times)))),"DEgenes.all.table.csv")
   
##############------------------------------------------------------------------------------------
#COUNT AND VISUALISE THE NUMBER OF DE GENES OVER TIME
##############------------------------------------------------------------------------------------

#Count number of DE genes for each timepoint 
dfGcount <- data.frame()
for(file in grep("_DESEq2_padj005.csv",list.files(),value=T)){
  dfGcount <- rbind(dfGcount,
                 data.frame("time"=gsub("05","0.5",gsub("hr","hap",unlist(strsplit(file, '_'))[1])),
                            "square"=unlist(strsplit(file, '_'))[2],
                            "Higher"=length(which(read.csv(file, header = T)$log2FoldChange>0)),
                            "Lower"=length(which(read.csv(file, header = T)$log2FoldChange<0))))
}
write.csv(dfGcount %>% arrange(factor(time, levels=gsub("05","0.5",gsub("hr","hap",times)))),"DEgenes.counts.csv")

#Plot number of DE genes
dfGcount$time <- factor(dfGcount$time, levels=gsub("05","0.5",gsub("hr","hap",times)))
dfGcount %>%
  mutate(Lower=Lower*-1) %>%
  melt() %>%
  ggplot(aes(x=time, y=value, fill=variable)) + 
  geom_bar(stat="identity", position="identity") +
  theme_bw(base_size = 14) +
  scale_fill_nejm() +
  theme(plot.title = element_text(hjust = 0.5, size=16),
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        axis.text.x=element_text(angle=60, hjust = 1.0),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16)) +
  labs(fill = "Expression") +
  xlab("") +
  ylab ("Number of DE genes") + #This is the axis label and title
  scale_y_continuous(breaks = scales::pretty_breaks(n = 8), limits=c(-2500,2500)) +
  guides(color=guide_legend(title="")) -> p
plot(p)

#Save
for(i in c(".pdf",".png")){ggsave(file=paste0("DEnumbers-over-time",i), plot=p, width=15, height=12, units = c("cm"), dpi = 300)}

##############------------------------------------------------------------------------------------
#Analyse GO-terms of DE genes for each timepoint
##############------------------------------------------------------------------------------------

dir.create(file.path(parDir,"003_DEgene_GOterms_timepoints"))
setwd(file.path(parDir,"003_DEgene_GOterms_timepoints"))

#Analyse GO terms of each timepoint
for(time in times){
  file <- grep(paste0("^",time,"_DESEq2_padj005.csv"),list.files(file.path(parDir,"002_DESeq2")),value=T)
  GO_enrich(filename=paste0("GO-BP_",time,"_DESEq2_padj005"),genes=read.csv(file.path(parDir,"002_DESeq2",file))$X,type="BP")
}

#Create a list with the GOterms of each timepoint
dfGcount <- data.frame()
for(file in grep("_DESEq2_padj005.csv",list.files(),value=T)){
  dfGO <- read.csv(file, header = T, row.names=1)
  time <- unlist(strsplit(file, '_'))[2]
  #Assign direction of DE to each gene
  dfDE <- read.csv(file.path(parDir,"002_DESeq2",paste0(time,"_DESEq2_padj005.csv")))
  for(GOID in dfGO$ID){
    geneID <- paste(dfDE[which(dfDE$X %in% unlist(strsplit(dfGO[grep(GOID,dfGO$ID),]$geneID,"|",fixed=T))),]$X,collapse="|")
    l2FC <- paste(dfDE[which(dfDE$X %in% unlist(strsplit(dfGO[grep(GOID,dfGO$ID),]$geneID,"|",fixed=T))),]$log2FC,collapse="|")
    df <- dfGO[grep(GOID,dfGO$ID),][c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","qvalue","Count")]
    df$geneID <- geneID
    df$l2FC <- l2FC
    df$time <- gsub("05","0.5",gsub("hr","hap",time))
    dfGcount <- rbind(dfGcount, df)
  }
}
write.csv(dfGcount %>% arrange(factor(time, levels=gsub("05","0.5",gsub("hr","hap",times)))), "GO-BP.list.csv")

#Count number of GO terms for each timepoint 
dfGcount <- data.frame()
for(file in grep("_DESEq2_padj005.csv",list.files(),value=T)){
  dfGcount <- rbind(dfGcount,
                    data.frame("time"=gsub("05","0.5",gsub("hr","hap",unlist(strsplit(file, '_'))[2])),
                               "GO-terms"=nrow(read.csv(file, header = T))))
}
write.csv(dfGcount %>% arrange(factor(time, levels=gsub("05","0.5",gsub("hr","hap",times)))),"GO-BP.counts.csv")
dfGcount$GO.terms <- log(dfGcount$GO.terms)
dfGcount <- rbind(dfGcount, data.frame(time="2hap",GO.terms=0))
dfGcount$time <- factor(dfGcount$time, levels=gsub("05","0.5",gsub("hr","hap",times)))
dfGcount %>%
  ggplot(aes(x=time, y=GO.terms, fill=time)) + 
  geom_bar(stat="identity", position="identity") +
  theme_bw(base_size = 14) +
  scale_fill_nejm() +
  theme(plot.title = element_text(hjust = 0.5, size=16),
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        axis.text.x=element_text(angle=60, hjust = 1.0),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16)) +
  labs(fill = "Timepoint") +
  xlab("") +
  ylab ("log(Number of BP GO-terms)") + #This is the axis label and title
  scale_y_continuous(breaks = scales::pretty_breaks(n = 8)) +
  guides(color=guide_legend(title="")) -> p
plot(p)

#Save
for(i in c(".pdf",".png")){ggsave(file=paste0("GO-BP-numbers-over-time",i), plot=p, width=15, height=12, units = c("cm"), dpi = 300)}

##############------------------------------------------------------------------------------------
#Condensation of GO-terms and plotting of Gantt Chart
##############------------------------------------------------------------------------------------

dir.create(file.path(parDir,"004_DEgene_GOterms_simplified_timepoints"))
setwd(file.path(parDir,"004_DEgene_GOterms_simplified_timepoints"))

#Create a table with all GO terms over time for grouping and plotting
dfGOall <- data.frame()
for(file in grep("_DESEq2_padj005.csv",list.files(file.path(parDir,"003_DEgene_GOterms_timepoints")),value=T)){
  dfGO <- read.csv(file.path(parDir,"003_DEgene_GOterms_timepoints",file), header = T, row.names=1)
  time <- unlist(strsplit(file, '_'))[2]
  #Assign direction of DE to each gene
  dfDE <- read.csv(file.path(parDir,"002_DESeq2",paste0(time,"_DESEq2_padj005.csv")))
  for(GOID in dfGO$ID){
    DE_l2FC <- paste(paste(dfDE[which(dfDE$X %in% unlist(strsplit(dfGO[grep(GOID,dfGO$ID),]$geneID,"|",fixed=T))),]$X,dfDE[which(dfDE$X %in% unlist(strsplit(dfGO[grep(GOID,dfGO$ID),]$geneID,"|",fixed=T))),]$log2FoldChange,sep="_"),collapse="|")
    df <- dfGO[grep(GOID,dfGO$ID),][c("ID","Description","geneID","GeneRatio","BgRatio","pvalue","p.adjust","qvalue","Count")]
    df$DE_l2FC <- DE_l2FC
    df$time <- gsub("05","0.5",gsub("hr","hap",time))
    dfGOall <- rbind(dfGOall, df)
  }
}

#Determine GO-term occurence groups based on identical enrichment times
dfGOall <- dfGOall %>%
  arrange(factor(time, levels=gsub("05","0.5",gsub("hr","hap",times)))) %>%
  group_by(ID) %>%
  mutate(times=paste(time, collapse=",")) %>% 
  mutate(group=paste(time, collapse=",")) %>% 
  ungroup() %>%
  mutate(group = as.numeric(factor(group)))

#Count the DE genes per group
dfGOall <- dfGOall %>%
  group_by(group) %>%
  #Condensed description headers per group
  mutate(Descriptions=paste(unique(Description), collapse=", ")) %>%
  #Condensed IDs per group
  mutate(IDs=paste(unique(ID), collapse=", ")) %>%
  #Calculate genes per group
  mutate(nGenesGroup=length(unique(unlist(strsplit(geneID,"|",fixed=T))))) %>%
  ungroup()

#Count the DE genes per timepoint per group
dfGOall <- dfGOall %>%
  group_by(group, time) %>%
  #Calculate genes for each timepoint per group
  mutate(nGenesGroupTime=length(unique(unlist(strsplit(geneID,"|",fixed=T))))) %>%
  mutate(nGenesGroupTimeHigher=length(which(as.numeric(sapply(unique(strsplit(unlist(strsplit(DE_l2FC,"|",fixed=T)),"_",fixed=T)), '[[', 2)) > 0))) %>%
  mutate(nGenesGroupTimeLower=length(which(as.numeric(sapply(unique(strsplit(unlist(strsplit(DE_l2FC,"|",fixed=T)),"_",fixed=T)), '[[', 2)) < 0))) %>%
  ungroup()

#Save
write.csv(dfGOall, "GO-BP.grouped.csv")

#Scrape REViGO
dfGOallRevigo <- data.frame()
for(identifier in unique(dfGOall$times)){
  goterms <- dfGOall %>% filter(times == identifier) %>% pull(IDs) %>% unique() %>% gsub(", ","\n",.)
  #Define homepage and pick submission form
  session <- html_session('http://revigo.irb.hr/')
  submitform <- html_form(session)[[1]]
  #Set values and submit GO IDs
  submitform <- set_values(submitform, goList = goterms)
  #Convert result to condensed GO table
  result <- submit_form(session, submitform)
  tables <- html_nodes(result, "table") 
  #convert the first table into a dataframe
  table <- html_table(tables[1])
  table <- as.data.frame(table)
  colnames(table) <- table[2,]
  table <- table[3:nrow(table),]
  #Add times information
  for(timepoint in unlist(strsplit(identifier,","))){
    group <- unique(dfGOall[dfGOall$times == identifier,]$group)
    table$group <- group
    table$time <- timepoint
    table$nGenesGroup <- dfGOall %>% filter(times == identifier) %>% filter(time == timepoint) %>% dplyr::select(nGenesGroup) %>% unique() %>% pull()
    table$nGenesGroupTime <- dfGOall %>% filter(times == identifier) %>% filter(time == timepoint) %>% dplyr::select(nGenesGroupTime) %>% unique() %>% pull()
    table$nGenesGroupTimeHigher <- length(which(as.numeric(unlist(lapply(unique(strsplit(unlist(strsplit(dfGOall[dfGOall$time==timepoint & dfGOall$group==group,]$DE_l2FC,"|",fixed=T)),"_")),'[[',2))) > 0))
    table$nGenesGroupTimeLower <- length(which(as.numeric(unlist(lapply(unique(strsplit(unlist(strsplit(dfGOall[dfGOall$time==timepoint & dfGOall$group==group,]$DE_l2FC,"|",fixed=T)),"_")),'[[',2))) < 0))
    dfGOallRevigo <- rbind(dfGOallRevigo, table)
  }
}

#Save
write.csv(dfGOallRevigo %>% arrange(factor(time, levels=gsub("05","0.5",gsub("hr","hap",times)))),"GO-BP.groups.csv")

#Filter for plotting
dfGOallRevigo <- dfGOallRevigo %>%
  dplyr::select(description, frequency, dispensability, uniqueness, group, time, nGenesGroupTimeHigher, nGenesGroupTimeLower) %>%
  unique() %>%
  dplyr::filter(dispensability <= 0.7) %>%
  group_by(group) %>%
  slice_max(order_by = frequency, n = 5) %>%
  mutate(descriptions=paste(unique(description),collapse=", ")) %>%
  ungroup() %>%
  dplyr::select(descriptions, group, time, nGenesGroupTimeHigher, nGenesGroupTimeLower) %>%
  unique()

#Make some terms shorter for better reading on y-axis
dfGOallRevigo$descriptions <- gsub("metabolic process", "m.p.",dfGOallRevigo$descriptions)
dfGOallRevigo$descriptions <- gsub("biosynthetic process", "b.p.",dfGOallRevigo$descriptions)
dfGOallRevigo$descriptions <- gsub("process", "p.",dfGOallRevigo$descriptions)
dfGOallRevigo$descriptions <- gsub("jasmonic acid", "JA",dfGOallRevigo$descriptions)
dfGOallRevigo$descriptions <- gsub("salicylic acid", "SA",dfGOallRevigo$descriptions)

#Create custom order for y-axis
o<-unique(dfGOallRevigo$descriptions)
order <- c(o[7],o[8],o[3],o[2],o[5],o[6],o[1],o[4],o[12],o[13],o[11],o[9],o[10],o[14],o[15],o[17],o[18],o[19],o[20],o[21],o[22],o[23],o[24],o[25],o[26],o[27],o[16])
#Add numerical to each group for sorting and plotting
dfGOallRevigo$descriptions <- paste0(sapply(dfGOallRevigo$descriptions, function(x) which(order == x)),") ",dfGOallRevigo$descriptions)
#Custom order on y axis
dfGOallRevigo$descriptions <- factor(dfGOallRevigo$descriptions, levels=paste0(seq(1,length(order)),") ",order))

dfGOallRevigo$time <- gsub("05","0.5",gsub("hr","hap",dfGOallRevigo$time))
p <- ggplot(dfGOallRevigo, aes(x=factor(time, levels=gsub("05","0.5",gsub("hr","hap",times))), y=descriptions, fill=log2((1+nGenesGroupTimeHigher)/(1+nGenesGroupTimeLower)))) + 
  geom_tile() +
  #coord_equal() +
  scale_fill_viridis_c() +
  #scale_fill_manual(values=colorRampPalette(pal_nejm("default")(5))(length(unique(dfGOallRevigo$group)))) +
  geom_text(aes(label=paste0(nGenesGroupTimeHigher,"/-",nGenesGroupTimeLower))) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, size=17),
        axis.text=element_text(size=15),
        axis.title=element_text(size=16),
        axis.text.x=element_text(angle=60, hjust = 1.0, size=25),
        axis.text.y=element_text(size=18),
        legend.text=element_text(size=20),
        legend.title=element_text(size=16),
        legend.position="bottom") +
  labs(fill = "log2(higher/lower DE)") +
  scale_y_discrete(labels=function(x) str_wrap(x, width=85)) +
  xlab("") +
  ylab ("") + #This is the axis label and title
  guides(color=guide_legend(title=""))
plot(p)

#Save
for(i in c(".pdf",".png")){ggsave(file=paste0("GO-BP-processes-over-time",i), plot=p, width=40, height=33, units = c("cm"), dpi = 300)}

##############------------------------------------------------------------------------------------
#Analyse GO-terms of DE genes for both waves
##############------------------------------------------------------------------------------------

dir.create(file.path(parDir,"005_DEgene_GOterms_waves"))
setwd(file.path(parDir,"005_DEgene_GOterms_waves"))

#First wave
genes <- c(read.csv(file.path(parDir,"002_DESeq2","0hr_DESEq2_padj005.csv"))$X,
           read.csv(file.path(parDir,"002_DESeq2","05hr_DESEq2_padj005.csv"))$X,
           read.csv(file.path(parDir,"002_DESeq2","1hr_DESEq2_padj005.csv"))$X)

GO_enrich(filename=paste0("GO-BP_first-wave_DESEq2_padj005"),genes=unique(genes),type="BP")

#Second wave
genes <- c(read.csv(file.path(parDir,"002_DESeq2","2hr_DESEq2_padj005.csv"))$X,
           read.csv(file.path(parDir,"002_DESeq2","3hr_DESEq2_padj005.csv"))$X,
           read.csv(file.path(parDir,"002_DESeq2","5hr_DESEq2_padj005.csv"))$X,
           read.csv(file.path(parDir,"002_DESeq2","7hr_DESEq2_padj005.csv"))$X,
           read.csv(file.path(parDir,"002_DESeq2","24hr_DESEq2_padj005.csv"))$X)

GO_enrich(filename=paste0("GO-BP_second-wave_DESEq2_padj005"),genes=genes,type="BP")

#Create a list with the GOterms of each wave
dfGcount <- data.frame()
for(file in grep("_DESEq2_padj005.csv",list.files(),value=T)){
  df <- read.csv(file, header = T, row.names=1)
  df$time <- unlist(strsplit(file, '_'))[2]
  dfGcount <- rbind(dfGcount, df)
}
write.csv(dfGcount, "GO-BP.list.csv")

##############------------------------------------------------------------------------------------
#Study differences between probing and feeding by PCA
##############------------------------------------------------------------------------------------

dir.create(file.path(parDir,"006_PCA-compare_probing-feeding"))
setwd(file.path(parDir,"006_PCA-compare_probing-feeding"))

###################
#LOAD PROBING DATA
###################

dfProbingFeeding <- data.frame()
for(time in c(t2,t3)){
  for(rep in c("r1","r2","r3","r4")){
    if(sum(head(read.table(file.path(htseqDirProbing,paste0("htseq-count_EPG-prob_",time,"_sq4_",rep,"_aphid.txt")), header = F), 33977)$V2) >= 100000){
      if(sum(head(read.table(file.path(htseqDirProbing,paste0("htseq-count_EPG-prob_",time,"_sq4_",rep,"_ctrl.txt")), header = F), 33977)$V2) >= 100000){
        #Read and normalise treatment and control tables
        dfTreat <- read.csv(file.path(htseqDirProbing,paste0("htseq-count_EPG-prob_",time,"_sq4_",rep,"_aphid.txt")),sep="\t",header=F)
        dfTreat <- head(dfTreat, 33977)
        dfTreat$V2 <- dfTreat$V2+1e-100
        dfTreat$V2 <- dfTreat$V2/(sum(dfTreat$V2)/1e6)
        
        dfCtrl <- read.csv(file.path(htseqDirProbing,paste0("htseq-count_EPG-prob_",time,"_sq4_",rep,"_ctrl.txt")),sep="\t",header=F)
        dfCtrl <- head(dfCtrl, 33977)
        dfCtrl$V2 <- dfCtrl$V2+1e-100
        dfCtrl$V2 <- dfCtrl$V2/(sum(dfCtrl$V2)/1e6)
        
        #Calculate l2FC (treatment vs control)
        dfTreat$V2 <- log2(dfTreat$V2/dfCtrl$V2)
    
        #Reformat to a single table
        row.names(dfTreat) <- dfTreat$V1
        dfTreat <- dfTreat[c('V2')]
        colnames(dfTreat) <- paste0(gsub("hr","hap",time),".",rep)
        dfProbingFeeding <- rbind(dfProbingFeeding, t(dfTreat))
      }
    }
  }
}

###################
#LOAD FEEDING DATA
###################

#Build table for feeding data
feedingFiles <- c("htseq-count_E2-feed-aphid_sq4_r1.txt",
                  "htseq-count_E2-feed-aphid_sq4_r2.txt",
                  "htseq-count_E2-feed-aphid_sq4_r3.txt")

for(file in feedingFiles){
  if(sum(head(read.csv(file.path(htseqDirFeeding,file),sep="\t",header=F), 33977)$V2) >= 100000){
      if(sum(head(read.csv(file.path(htseqDirFeeding,gsub("aphid","ctrl",file)),sep="\t",header=F), 33977)$V2) >= 100000){
        #Read and normalise treatment and control tables
        dfTreat <- read.csv(file.path(htseqDirFeeding,file),sep="\t",header=F)
        dfTreat <- head(dfTreat, 33977)
        dfTreat$V2 <- dfTreat$V2+1e-100
        dfTreat$V2 <- dfTreat$V2/(sum(dfTreat$V2)/1e6)
        
        dfCtrl <- read.csv(file.path(htseqDirFeeding,gsub("aphid","ctrl",file)),sep="\t",header=F)
        dfCtrl <- head(dfCtrl, 33977)
        dfCtrl$V2 <- dfCtrl$V2+1e-100
        dfCtrl$V2 <- dfCtrl$V2/(sum(dfCtrl$V2)/1e6)
        
        #Calculate l2FC (treatment vs control)
        dfTreat$V2 <- log2(dfTreat$V2/dfCtrl$V2)
        
        #Reformat to a single table
        row.names(dfTreat) <- dfTreat$V1
        dfTreat <- dfTreat[c('V2')]
        colnames(dfTreat) <- paste0("feeding.",gsub(".txt","",unlist(strsplit(file,"_"))[4]))
        dfProbingFeeding <- rbind(dfProbingFeeding, t(dfTreat))
    }
  }
}

#NA to 0
dfProbingFeeding[is.na(dfProbingFeeding)] <- 0

###################
#Perform PCA analysis to compare the samples using t2, t3 DE genes
###################
genes <- c()
for(time in c(t2,t3)){
  file <- grep(paste0("^",time,"_DESEq2_padj005.csv"),list.files(file.path(parDir,"002_DESeq2")),value=T)
  genes <- c(genes,read.csv(file.path(parDir,"002_DESeq2",file))$X)
}
res.pca <- prcomp(dfProbingFeeding[which(colnames(dfProbingFeeding) %in% unique(genes))],  scale=T)
plotPCA(pca=res.pca,filename="PCA_0.5hap_1hap_DE-genes")

###################
#Perform PCA analysis with subsets of DE genes (of the size as the number of t2,t3 DE genes)
###################
random.genes <- c()
for(time in c(t1,t2,t3,t4,t5,t6,t7,t8)){
  file <- grep(paste0("^",time,"_DESEq2_padj005.csv"),list.files(file.path(parDir,"002_DESeq2")),value=T)
  random.genes <- c(random.genes,read.csv(file.path(parDir,"002_DESeq2",file))$X)
}
random.genes <- unique(setdiff(random.genes,genes))

dfRandomGenesVarPercent <- data.frame()
for(i in seq(1,1000)){
  dfRandomGenesVarPercent <- rbind(dfRandomGenesVarPercent, t(data.frame(get_eig(prcomp(dfProbingFeeding[which(colnames(dfProbingFeeding) %in% sample(random.genes, length(genes)))],  scale=T))[,2])))
}
colnames(dfRandomGenesVarPercent) <- paste0("PC",seq(1,ncol(dfRandomGenesVarPercent)))
dfRandomGenesVarPercent <- melt(dfRandomGenesVarPercent[c(paste0("PC",seq(1,10)))])

#Mean contribution to PC1
mean(dfRandomGenesVarPercent[dfRandomGenesVarPercent$variable=="PC1",]$value)
sd(dfRandomGenesVarPercent[dfRandomGenesVarPercent$variable=="PC1",]$value)

p <- ggplot(dfRandomGenesVarPercent, aes(x=variable, y=value)) +
  geom_boxplot(aes(fill=variable)) +
  scale_fill_brewer(palette="Set3") +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, size=16),
        axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        axis.text.x=element_text(angle=60, hjust = 1.0),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.position="") +
  xlab("") +
  ylab("Variance Percent") +
  ylim(c(-1,30))
plot(p)
#Save
for(i in c(".pdf",".png")){ggsave(file=paste0("PC_variance-percent_random-DE-subsets",i), plot=p, width=20, height=12, units = c("cm"), dpi = 300)}

###################
#Find genes which separate probing from feeding; hereafter called feeding.genes
###################
feeding.genes <- c()
dfContrib <- as.data.frame(get_pca_var(res.pca)$contrib)
for(gene in row.names(dfContrib)){
  df <- dfContrib[row.names(dfContrib) %in% gene,]
  if(max(df) == df$Dim.1){
    feeding.genes <- c(feeding.genes, gene)
  }
}
GO_enrich(filename=paste0("GO-BP_feeding-genes"),genes=unique(feeding.genes),type="BP")

#Control contributions of selected feeding.genes to the PC1
res.pca <- prcomp(dfProbingFeeding[which(colnames(dfProbingFeeding) %in% unique(feeding.genes))], scale=T)
plotPCA(pca=res.pca,filename="PCA_feeding-genes")

#Visualise the log2FCs of PCA selected DE genes in probing and feeding samples
df <- dfProbingFeeding[(colnames(dfProbingFeeding) %in% feeding.genes)]
df$type <- c(rep("0.5hap",4),rep("1hap",4),rep("feeding",3))

df <- df %>%
  melt() %>%
  group_by(type, variable) %>%
  mutate(l2FC=mean(value)) %>%
  dplyr::select(type, variable, l2FC) %>%
  unique()

df <- merge(df, getBM(attributes = c("ensembl_gene_id","tair_symbol"), filters="ensembl_gene_id", values=df$variable, mart= mart), by.x="variable", by.y="ensembl_gene_id")
df$Direction <- ifelse(df$variable %in% df[df$type=="1hap",][which(df[df$type=="1hap",]$l2FC > 0),]$variable, "Higher DE", "Lower DE")
df$type <- factor(df$type, levels=c("feeding", "0.5hap", "1hap"))

p <- ggplot(df, aes(x=type, y=l2FC)) +
  geom_boxplot(aes(fill=type),alpha=0.5) +
  geom_line(aes(group=variable),alpha=0.25,linetype="dashed") +
  geom_point(aes(color=variable)) +
  scale_fill_nejm() +
  scale_color_manual(values=colorRampPalette(pal_nejm("default")(5))(length(unique(df$variable)))) +
  facet_wrap(~Direction) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, size=20),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        axis.text.x=element_text(angle=60,hjust = 1.0),
        strip.text = element_text(size = 20),
        legend.position = "") +
  xlab("") +
  ylab("log2FC vs control at aphid site") +
  ylim(c(-4,4))
plot(p)
#Save
for(i in c(".pdf",".png")){ggsave(file=paste0("feeding_genes-changes_feeding-0.5hap-1hap",i), plot=p, width=25, height=15, units = c("cm"), dpi = 300)}

#Export table with selected DE genes
df <- dfProbingFeeding[(colnames(dfProbingFeeding) %in% feeding.genes)]
df$type <- c(rep("0.5hap",4),rep("1hap",4),rep("feeding",3))

df <- df %>%
  melt() %>%
  group_by(type, variable) %>%
  mutate(l2FC=mean(value)) %>%
  dplyr::select(type, variable, l2FC) %>%
  unique()

df <- dcast(df, variable ~ type)
df <- merge(df, getBM(attributes = c("ensembl_gene_id","tair_symbol"), filters="ensembl_gene_id", values=df$variable, mart= mart), by.x="variable", by.y="ensembl_gene_id")
df$Direction <- ifelse(df$`1hap` > 0, "Higher DE", "Lower DE")
write.csv(df, "PCA.feeding.genes.l2FCmeans.csv")

##############------------------------------------------------------------------------------------
#Network analysis
##############------------------------------------------------------------------------------------

dir.create(file.path(parDir,"007_feeding-network"))
setwd(file.path(parDir,"007_feeding-network"))

write.csv(data.frame(feeding.genes),"feeding-genes_forString.csv")
#########
#Perform web analysis of feeding.genes listed in the "feeding-genes_forString.csv" file on https://string-db.org
#1. Select: Multiple Proteins, Organism: A.thaliana and upload gene list without changing advanced settings
#2. Create network
#3. Settings: Select 'Full network', set 1st and 2nd shell to 'none'
#4. Analysis: Download GO-terms (Biological process table) and other tables of interest
#5. Clusters: Perform MCL clustering using standard settings (inflation parameter 3)
#6. Clusters: Download TSV file of the clustering to the 007_feeding-network folder: 'clusters string_MCL_clusters.tsv'
#7. Exports: Download TSV files of the network to the 007_feeding-network folder: 'string_interactions.tsv'
#8. Continue with the code below
#########

#Read the STRING web exported tsv tables
stringNw <- read.delim("string_interactions.tsv",header=T,sep="\t",fill=T)
stringCl <- read.delim("string_MCL_clusters.tsv",header=T,sep="\t",fill=T)

#Add TAIR loci to the STRING network
stringNw$from <- as.character(gsub("[.].*$","",gsub("3702.","",stringNw$node1_string_id)))
stringNw$to <- as.character(gsub("[.].*$","",gsub("3702.","",stringNw$node2_string_id)))

stringNw <- stringNw %>%
  dplyr::select(X.node1, node2, from, to, coexpression, experimentally_determined_interaction, database_annotated, automated_textmining, combined_score)

#Reformat stringNw dataframe for network calculations below
df <- stringNw
colnames(df) <- c("node2", "X.node1", "to", "from", "coexpression", "experimentally_determined_interaction", "database_annotated", "automated_textmining", "combined_score")
stringNw <- unique(rbind(stringNw, df))

#Add cluster information to the STRING network, select genes which belong to clusters
stringCl$protein.identifier <- gsub("[.].*$","",gsub("3702.","",stringCl$protein.identifier))
stringNw$cluster.number <- unlist(lapply(stringNw$from, function(x) stringCl[stringCl$protein.identifier==x,]$cluster.number))
stringNw <- stringNw %>% dplyr::select("from", "to", "cluster.number") %>% unique()

#Merge with genes that don't have any interactions and give them the cluster number 0
stringNw <- rbind(stringNw, data.frame("from"=setdiff(feeding.genes,c(stringNw$from, stringNw$to)),
                                       "to"=setdiff(feeding.genes,c(stringNw$from, stringNw$to)),
                                       "cluster.number"=0))
  
##Analyse each cluster
stringNwClusterScores <- data.frame()
for(cluster in unique(stringNw$cluster.number)){
  df <- stringNw[stringNw$cluster.number == cluster,][c(1,2)]
  g <- simplify(graph.data.frame(df, directed=FALSE))
  df$degreecentrality <- unlist(lapply(df$from, function(x) as.numeric(degree(g)[x])))
  df$betweennesscentrality <- unlist(lapply(df$from, function(x) as.numeric(betweenness(g)[x])))
  df$eigencentrality <- unlist(lapply(df$from, function(x) as.numeric(eigen_centrality(g)$vector[x])))
  df <- df %>% dplyr::select(from, degreecentrality, betweennesscentrality, eigencentrality) %>% unique()
  stringNwClusterScores <- rbind(stringNwClusterScores, df)
}
stringNw <- merge(stringNw, stringNwClusterScores, by="from")

#Add TAIR names
df <- getBM(attributes = c("ensembl_gene_id","tair_symbol"), filters="ensembl_gene_id", values=unique(c(stringNw$from,stringNw$to)), mart= mart)
stringNw$A <- unlist(lapply(stringNw$from, function(x) unique(ifelse(df[df$ensembl_gene_id==x,]$tair_symbol=="",x,df[df$ensembl_gene_id==x,]$tair_symbol[1]))))
stringNw$B <- unlist(lapply(stringNw$to, function(x) unique(ifelse(df[df$ensembl_gene_id==x,]$tair_symbol=="",x,df[df$ensembl_gene_id==x,]$tair_symbol[1]))))

#Export for cytoscape
write.table(stringNw,"for-cytoscape-interactions_table_STRING_MCL.txt",row.names=FALSE,sep="\t", quote = FALSE)

#GOenrich clusters
for(cluster in unique(stringNw$cluster.number)){
  GO_enrich(filename=paste0("GO-BP_MCL-cluster",cluster),genes=unique(c(stringNw[stringNw$cluster.number==cluster,]$from,stringNw[stringNw$cluster.number==cluster,]$to)),type="BP")
  GO_enrich(filename=paste0("GO-MF_MCL-cluster",cluster),genes=unique(c(stringNw[stringNw$cluster.number==cluster,]$from,stringNw[stringNw$cluster.number==cluster,]$to)),type="MF")
}

##############------------------------------------------------------------------------------------
#Analyse spatial expression patterns of selected feeding genes: Network clusters and affinity propagation clustering
##############------------------------------------------------------------------------------------

dir.create(file.path(parDir,"008_spatial-plotting"))
setwd(file.path(parDir,"008_spatial-plotting"))

#Probing; calculate RPM
dfCts <- data.frame()
for(file in grep("aphid",list.files(htseqDirProbing),value=T)){
  df <- read.table(file.path(htseqDirProbing,file), header = F)
  #Select rows corresponding to genes
  df <- head(df, 33977)
  #Normalise
  df$V2 <- df$V2+1e-100
  df$V2 <- df$V2/(sum(df$V2)/1e6)
  #Select string genes
  df <- df[df$V1 %in% feeding.genes,]
  #Add additional information
  df$type <- "probing"
  df$Timepoint <- unlist(strsplit(gsub("htseq-count_","",gsub(".txt","",file)),'_'))[2]
  df$square <- unlist(strsplit(gsub("htseq-count_","",gsub(".txt","",file)),'_'))[3]
  df$replicate <- unlist(strsplit(gsub("htseq-count_","",gsub(".txt","",file)),'_'))[4]
  dfCts <- rbind(dfCts, df)
}
dfCts <- dfCts %>% dplyr::filter(Timepoint %in% c(t2,t3))

#Feeding; calculate RPM
for(file in grep("aphid",list.files(htseqDirFeeding),value=T)){
  df <- read.table(file.path(htseqDirFeeding,file), header = F)
  #Select rows corresponding to genes
  df <- head(df, 33977)
  #Normalise
  df$V2 <- df$V2+1e-100
  df$V2 <- df$V2/(sum(df$V2)/1e6)
  #Select string genes
  df <- df[df$V1 %in% feeding.genes,]
  #Add additional information
  df$type <- "feeding"
  df$Timepoint <- "feeding"
  df$square <- unlist(strsplit(gsub("htseq-count_","",gsub(".txt","",file)),'_'))[2]
  df$replicate <- unlist(strsplit(gsub("htseq-count_","",gsub(".txt","",file)),'_'))[3]
  dfCts <- rbind(dfCts, df)
}

#Calculate normalised counts
dfCts <- dfCts %>%
  group_by(V1, Timepoint, square) %>% #Mean RPM of replicates per square and timepoint (e.g. 0.5hap, 1hap, feeding)
  mutate(meanCts=mean(V2)) %>%
  ungroup() %>%
  unique() %>%
  group_by(V1) %>%
  mutate(l2vsMean=log2(meanCts/mean(V2))) %>% #Mean RPM per square and timepoint normalised with mean RPM of all Timepoints
  ungroup()

#Merge with MCL cluster
dfCts <- merge(dfCts,stringNw %>% dplyr::select("from","cluster.number") %>% unique(),by.x="V1",by.y="from")

#Merge with TAIR name
df <- getBM(attributes = c("ensembl_gene_id","tair_symbol"), filters="ensembl_gene_id", values=unique(dfCts$V1), mart= mart)
dfCts$name <- unlist(lapply(dfCts$V1, function(x) unique(ifelse(df[df$ensembl_gene_id==x,]$tair_symbol=="",x,df[df$ensembl_gene_id==x,]$tair_symbol[1]))))

#Modify columns for plotting
dfCts$Timepoint <- gsub("hr","hap",dfCts$Timepoint)

#Calculate mean for each cluster
dfCts <- dfCts %>% group_by(cluster.number, square, Timepoint) %>% mutate(avgl2vsMean=mean(l2vsMean)) %>% ungroup()

#Plot
p <- ggplot(dfCts) +
  geom_line(aes(x=square, y=l2vsMean, group=interaction(name,Timepoint), color=Timepoint),size=0.6,alpha=0.2) +
  geom_line(aes(x=square, y=avgl2vsMean, group=interaction(name,Timepoint), color=Timepoint),size=1, linetype="dashed") +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, size=20),
        axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        axis.text.x=element_text(angle=60,hjust = 1.0,size=10),
        strip.text = element_text(size = 16),
        legend.position = "right",
        legend.text=element_text(size=16),
        legend.title=element_text(size=16)) +
  xlab("") +
  ylab("log2FC vs.mean RPM\nof timepoints") +
  facet_wrap(~cluster.number,scales="free_y", nrow = 2) +
  scale_color_nejm() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot(p)
#Save
for(i in c(".pdf",".png")){ggsave(file=paste0("DE-expression_patterns_probing-vs-feeding_MCL-cluster",i), plot=p, width=28, height=10, units = c("cm"), dpi = 300)}

#Prepare dataframe for expression pattern clustering
df <- dfCts
df$var <- paste0(df$Timepoint,df$square)
df <- dcast(data = df[c("V1","l2vsMean","var")],formula = V1~var,value.var ="l2vsMean",fun.aggregate=mean)
row.names(df) <- df$V1
df <- df[seq(2,ncol(df))]
df <- apcluster(s=corSimMat(df, sel=NA, r=1, signed=TRUE, method="pearson"))

dfAP <- data.frame()
for(cluster in seq(1,length(df@clusters))){
  dfAP <- rbind(dfAP, data.frame("V1"=names(df[[cluster]]), "ap.cluster"=cluster))
}
dfCts <- merge(dfCts, dfAP, by="V1")

#Calculate mean for each affinity propagation cluster cluster
dfCts <- dfCts %>% group_by(ap.cluster, square, Timepoint) %>% mutate(avgl2vsMean=mean(l2vsMean)) %>% ungroup()

#Plot
p <- ggplot(dfCts) +
  geom_line(aes(x=square, y=l2vsMean, group=interaction(name,Timepoint), color=Timepoint),size=0.6,alpha=0.2) +
  geom_line(aes(x=square, y=avgl2vsMean, group=interaction(name,Timepoint), color=Timepoint),size=1, linetype="dashed") +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, size=20),
        axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        axis.text.x=element_text(angle=60,hjust = 1.0,size=10),
        strip.text = element_text(size = 16),
        legend.position = "right",
        legend.text=element_text(size=16),
        legend.title=element_text(size=16)) +
  xlab("") +
  ylab("log2FC vs. mean RPM\nof timepoints") +
  facet_wrap(~ap.cluster,scales="free_y", nrow = 2) +
  scale_color_nejm() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot(p)
#Save
for(i in c(".pdf",".png")){ggsave(file=paste0("DE-expression_patterns_probing-vs-feeding_AP-cluster",i), plot=p, width=28, height=10, units = c("cm"), dpi = 300)}

#Print each affinity propagation cluster
for(cluster in unique(dfCts$ap.cluster)){
  #Calculate average l2FC expression per gene across leaf
  dfCtsA <- dfCts %>%
    filter(ap.cluster %in% c(cluster)) %>%
    dplyr::select("V1","name","meanCts","V2","Timepoint","square","replicate","ap.cluster") %>%
    group_by(V1, Timepoint) %>%
    mutate(l2vsMean=log2(meanCts/mean(V2))) %>% #L2FC expression of each gene in each square normalised to its RPM across leaf
    ungroup() %>%
    dplyr::select("V1","name","l2vsMean","Timepoint","square","ap.cluster") %>%
    unique()

  dfCtsA$name <- factor(dfCtsA$name, levels=sapply(strsplit(str_sort(unique(paste0("Cl.",dfCtsA$ap.cluster,": ",dfCtsA$name)), numeric=T),": "), '[[', 2))
  p <- ggplot(dfCtsA, aes(x=square, y=name, fill=l2vsMean, color=ap.cluster)) + 
    geom_tile() +
    coord_equal() +
    scale_fill_viridis_c() +
    #geom_text(aes(label=paste0(nGenesGroupTimeHigher,"/-",nGenesGroupTimeLower))) +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, size=17),
          axis.text=element_text(size=20),
          axis.title=element_text(size=20),
          axis.text.x=element_text(angle=60, hjust = 1.0, size=20),
          axis.text.y=element_text(size=20),
          strip.text = element_text(size = 20),
          legend.text=element_text(size=20),
          legend.title=element_text(size=20),
          legend.position="right") +
    labs(fill = "log2FC vs. mean RPM (gene)") +
    #scale_y_discrete(labels=function(x) str_wrap(x, width=85)) +
    xlab("") +
    facet_grid(~Timepoint) +
    ylab ("") + #This is the axis label and title
    guides(color=guide_legend(title=""))
  plot(p)
  #Save
  for(i in c(".pdf",".png")){ggsave(file=paste0("DE-expression_patterns_probing-vs-feeding_ap-cluster",cluster,i), plot=p, width=35, height=20, units = c("cm"), dpi = 300)}
  
  #GOenrich
  GO_enrich(filename=paste0("GO-BP_ap-cluster",cluster),genes=unique(dfCtsA$V1),type="BP")
  GO_enrich(filename=paste0("GO-MF_ap-cluster",cluster),genes=unique(dfCtsA$V1),type="MF")
}

#Export for cytoscape
stringNw <- merge(stringNw, unique(dfCts[c("V1","ap.cluster")]), by.x="from", by.y="V1")
write.table(stringNw,"for-cytoscape-interactions_table_STRING_MCL_AP.txt",row.names=FALSE,sep="\t", quote = FALSE)

