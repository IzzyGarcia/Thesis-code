DESEQ2 downstream analysis following inital RNA-seq 

#downloading R packages
library(DESeq2) #DESeq2 DE analysis
library(gplots) #doing the heatmap and PCA plots
library(RColorBrewer) #doing the heatmap plot
library(calibrate) #doing the PCA
library(EnhancedVolcano) #obviously for the volcano plots (bioconductor)
library(tidyverse) #helps with the annotation join
library(AnnotationDbi) #for annotating the ensembl IDs to gene symbols (bioconductor)
library(org.Hs.eg.db) #for annotating the ensembl IDs to gene symbols (bioconductor)

#Read in the count data
countdata <- read.delim("counts_PE.txt", header=TRUE)
x=countdata[,-1]
rownames(x)=countdata[,1]
countdata=x
countdata <- as.matrix(countdata)
#head(countdata)

#Setting contrasts
condition <- factor(c(rep("Control",3),rep("Transformed",3),rep("ZYFS.GFP",3),rep("ZFYL.GFP",3),rep("ZYFS.HA",3),rep("ZFYL.HA",3)))
(coldata <- data.frame(row.names=colnames(countdata), condition))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds <- DESeq(dds) 

#Outputs based on contrasts
#res<-results(dds) #res is some sort of object!!
#resultsNames(dds) will get a list of contrast names to call out with

#resultsNames(dds)
#"Intercept"
#"condition_Transformed_vs_Control" 
#"condition_ZFYL.GFP_vs_Control"    
#"condition_ZFYL.HA_vs_Control"    
#"condition_ZYFS.GFP_vs_Control"
#"condition_ZYFS.HA_vs_Control" 

#res <- results(dds, name="")
res.TransformedvsCont <- results(dds, name="condition_Transformed_vs_Control")
  res.TransformedvsCont.DF<-(as.data.frame(res.TransformedvsCont[which(res.TransformedvsCont$padj < 0.05),]))
res.ZYFS_GFFvsCont <- results(dds, name="condition_ZYFS.GFP_vs_Control")
  res.ZYFS_GFFvsCont.DF<-(as.data.frame(res.ZYFS_GFFvsCont[which(res.ZYFS_GFFvsCont$padj < 0.05),]))
res.ZYFS_HAvsCont <- results(dds, name="condition_ZYFS.HA_vs_Control")
  res.ZYFS_HAvsCont.DF<-(as.data.frame(res.ZYFS_HAvsCont[which(res.ZYFS_HAvsCont$padj < 0.05),]))
res.ZFYL_GFFvsCont <- results(dds, name="condition_ZFYL.GFP_vs_Control")
  res.ZFYL_GFFvsCont.DF<-(as.data.frame(res.ZFYL_GFFvsCont[which(res.ZFYL_GFFvsCont$padj < 0.05),]))
res.ZFYL_HAvsCont <- results(dds, name="condition_ZFYL.HA_vs_Control")
  res.ZFYL_HAvsCont.DF<-(as.data.frame(res.ZFYL_HAvsCont[which(res.ZFYL_HAvsCont$padj < 0.05),]))
  
#Annotate the dataframes with gene symbols
TransformedvsCont.annotation <- AnnotationDbi::select(org.Hs.eg.db, keys = row.names(res.TransformedvsCont.DF), keytype = 'ENSEMBL', columns = c("SYMBOL","MAP"))
  res.TransformedvsCont.DF <- left_join(res.TransformedvsCont.DF %>% mutate(ENSEMBL = rownames(res.TransformedvsCont.DF)), TransformedvsCont.annotation) 
ZYFS_GFFvsCont.annotation <- AnnotationDbi::select(org.Hs.eg.db, keys = row.names(res.ZYFS_GFFvsCont.DF), keytype = 'ENSEMBL', columns = c("SYMBOL","MAP"))
  res.ZYFS_GFFvsCont.DF <- left_join(res.ZYFS_GFFvsCont.DF %>% mutate(ENSEMBL = rownames(res.ZYFS_GFFvsCont.DF)), ZYFS_GFFvsCont.annotation) 
ZYFS_HAvsCont.annotation <- AnnotationDbi::select(org.Hs.eg.db, keys = row.names(res.ZYFS_HAvsCont.DF), keytype = 'ENSEMBL', columns = c("SYMBOL","MAP"))
  res.ZYFS_HAvsCont.DF <- left_join(res.ZYFS_HAvsCont.DF %>% mutate(ENSEMBL = rownames(res.ZYFS_HAvsCont.DF)), ZYFS_HAvsCont.annotation) 
ZFYL_GFFvsCont.annotation <- AnnotationDbi::select(org.Hs.eg.db, keys = row.names(res.ZFYL_GFFvsCont.DF), keytype = 'ENSEMBL', columns = c("SYMBOL","MAP"))
  res.ZFYL_GFFvsCont.DF <- left_join(res.ZFYL_GFFvsCont.DF %>% mutate(ENSEMBL = rownames(res.ZFYL_GFFvsCont.DF)), ZFYL_GFFvsCont.annotation) 
ZFYL_HAvsCont.annotation <- AnnotationDbi::select(org.Hs.eg.db, keys = row.names(res.ZFYL_HAvsCont.DF), keytype = 'ENSEMBL', columns = c("SYMBOL","MAP"))
  res.ZFYL_HAvsCont.DF <- left_join(ZFYL_HAvsCont.annotation, res.ZFYL_HAvsCont.DF %>% mutate(ENSEMBL = rownames(res.ZFYL_HAvsCont.DF)))
  
#Annotated write out
write.csv(res.TransformedvsCont.DF,'Transformed_vs_Control_results.csv')
write.csv(res.ZYFS_GFFvsCont.DF,'ZYFS_GFF_vs_Control_results.csv')
write.csv(res.ZYFS_HAvsCont.DF,'ZYFS_HA_vs_Control_results.csv')
write.csv(res.ZFYL_GFFvsCont.DF,'ZFYL_GFF_vs_Control_results.csv')
write.csv(res.ZFYL_HAvsCont.DF,'ZFYL_HA_vs_Control_results.csv')

#Make a dispersal plot
plotDispEsts(dds, main="Dispersion plot")

#Make heatmap
rld <- rlogTransformation(dds)
sampleDists <- as.matrix(dist(t(assay(rld))))
(mycols <- brewer.pal(8, "Dark2")[1:length(unique(condition))])
## [1] "#1B9E77" "#D95F02" "#7570B3" "#E7298A" "#66A61E" "#E6AB02"
heatmap.2(as.matrix(sampleDists), key=F, trace="none",

          col=colorpanel(100, "blue", "green"),
          ColSideColors=mycols[condition], RowSideColors=mycols[condition],
          margin=c(10, 10), main="Sample Distance Matrix")
          
#PCA plot
rld_pca <- function (rld, intgroup = "condition", ntop = 5000, colors=NULL, legendpos="topright", main="PCA Biplot", textcx=0.5, ...) {
  rv = rowVars(assay(rld))
  select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(rld)[select, ]))
  fac = factor(apply(as.data.frame(colData(rld)[, intgroup, drop = FALSE]), 1, paste, collapse = " : "))
  if (is.null(colors)) {
    if (nlevels(fac) >= 3) {
      colors = brewer.pal(nlevels(fac), "Paired")
    } else {
      colors = c("black", "red")
    }
  }
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC2 (",as.character(pc2var),"%)")
  plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=21, xlab=pc1lab, ylab=pc2lab, main=main, cex=1,...)
  with(as.data.frame(pca$x), textxy(PC1, PC2, labs=rownames(as.data.frame(pca$x)), cex=textcx))
  legend(legendpos, legend=levels(fac), col=colors, pch=20)
  # rldyplot(PC2 ~ PC1, groups = fac, data = as.data.frame(pca$rld),
  # pch = 16, cerld = 2, aspect = "iso", col = colours, main = draw.key(key = list(rect = list(col = colours),
  # terldt = list(levels(fac)), rep = FALSE)))
}
rld_pca(rld, colors=mycols, intgroup="condition", xlim=c(-45, 60), ylim=c(-30,30))

#MA plots
DESeq2::plotMA(res.TransformedvsCont, main = "Transformed vs Control", ylim=c(-5,5))
DESeq2::plotMA(res.ZYFS_GFFvsCont, main = "ZYFS GFF vs Control", ylim=c(-5,5))
DESeq2::plotMA(res.ZYFS_HAvsCont, main = "ZYFS HA vs Control", ylim=c(-5,5))
DESeq2::plotMA(res.ZFYL_GFFvsCont, main = "ZFYL GFF vs Control", ylim=c(-5,5))
DESeq2::plotMA(res.ZFYL_HAvsCont, main = "ZFYL HA vs Control", ylim=c(-5,5))

#Volcano plot (plotted on padj < 0.05 filtered data)
EnhancedVolcano(res.TransformedvsCont.DF,
  lab = res.TransformedvsCont.DF$SYMBOL,
  x = 'log2FoldChange',
  y = 'padj',
  xlim = c(-4, 4),
  title = 'Transformed vs Control',
  pCutoff = 0.05,
  FCcutoff = 1.0
  )
  
 EnhancedVolcano(res.ZYFS_GFFvsCont.DF,
  lab = res.ZYFS_GFFvsCont.DF$SYMBOL,
  x = 'log2FoldChange',
  y = 'padj',
  xlim = c(-4, 4),
  title = 'ZYFS GFF vs Control',
  pCutoff = 0.05,
  FCcutoff = 1.0
  )
  
  EnhancedVolcano(res.ZYFS_HAvsCont.DF,
  lab = res.ZYFS_HAvsCont.DF$SYMBOL,
  x = 'log2FoldChange',
  y = 'padj',
  xlim = c(-4, 4),
  title = 'ZYFS HA vs Control',
  pCutoff = 0.05,
  FCcutoff = 1.0
  )
  
  EnhancedVolcano(res.ZFYL_GFFvsCont.DF,
  lab = res.ZFYL_GFFvsCont.DF$SYMBOL,
  x = 'log2FoldChange',
  y = 'padj',
  xlim = c(-4, 4),
  title = 'ZFYL GFF vs Control',
  pCutoff = 0.05,
  FCcutoff = 1.0
  )
  
  EnhancedVolcano(res.ZFYL_HAvsCont.DF,
  lab = res.ZFYL_HAvsCont.DF$SYMBOL,
  x = 'log2FoldChange',
  y = 'padj',
  xlim = c(-4, 4),
  title = 'ZFYL HA vs Control',
  pCutoff = 0.05,
  FCcutoff = 1.0
  )
  
