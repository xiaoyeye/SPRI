library(clusterProfiler)
library(org.Hs.eg.db)

prescreen_type = 'results' #prescreen_negcorr
prefix = 'BC2020_tic_pvals'

all_prefix = paste(prescreen_type,'/BC2020/reference_genelist.csv',sep = "", collapse = NULL)
read_path = paste(prescreen_type,'/BC2020/',prefix,'.csv',sep = "", collapse = NULL)
GO_save_path = paste(prescreen_type,'/clusterProfiler/BC2020/',prefix,'_GO_refall.csv',sep = "", collapse = NULL)
KEGG_save_path = paste(prescreen_type,'/clusterProfiler/BC2020/',prefix,'_KEGG_refall.csv',sep = "", collapse = NULL)

BC.df<-read.csv(read_path)
gene<- BC.df$gene[which(as.numeric(BC.df$qvalue)<0.05)] #FDR = 0.05
BC.df<-read.csv(all_prefix)
geneList<-BC.df$gene

# SYMBOL- ENSEMBL,ENTREZID
gene.df <- bitr(gene, fromType = "ENSEMBL",
                toType = "ENTREZID",
                OrgDb = org.Hs.eg.db)

geneList.df <- bitr(geneList, fromType = "ENSEMBL",
                    toType = "ENTREZID",
                    OrgDb = org.Hs.eg.db)

head(gene.df$ENTREZID)

#GO over-representation test
ego <- enrichGO(gene = gene.df$ENTREZID,
               universe = geneList.df$ENTREZID,
               OrgDb = org.Hs.eg.db,
               keyType = 'ENTREZID',
               ont = "ALL",
               pvalueCutoff  = 0.05,
               qvalueCutoff  = 0.05,
               pAdjustMethod = "BH")
head(ego)

# save results
write.csv(ego,file=GO_save_path,quote=F,row.names = F)

# KEGG analysis
kk <- enrichKEGG(gene = gene.df$ENTREZID,
                 universe = geneList.df$ENTREZID,
                 organism = 'hsa',
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05)

#save results
write.csv(kk,file=KEGG_save_path,quote=F,row.names = F)