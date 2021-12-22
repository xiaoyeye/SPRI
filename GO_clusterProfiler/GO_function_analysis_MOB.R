library(clusterProfiler)
library(org.Mm.eg.db)

prescreen_type = 'results' #prescreen_negcorr
MOB = 'MOB11' #MOB12
prefix = 'MOB11_tic_pvals'

all_prefix = paste(prescreen_type,'/',MOB,'/reference_genelist.csv',sep = "", collapse = NULL)
read_path = paste(prescreen_type,'/',MOB,'/',prefix,'.csv',sep = "", collapse = NULL)
GO_save_path = paste(prescreen_type,'/clusterProfiler/',MOB,'/',prefix,'_GO_refall.csv',sep = "", collapse = NULL)
KEGG_save_path = paste(prescreen_type,'/clusterProfiler/',MOB,'/',prefix,'_KEGG_refall.csv',sep = "", collapse = NULL)

BC.df<-read.csv(read_path)
gene<- BC.df$gene[which(as.numeric(BC.df$qvalue)<0.05)] #FDR = 0.05
BC.df<-read.csv(all_prefix)
geneList<-BC.df$gene

# SYMBOL- ENSEMBL,ENTREZID
gene.df <- bitr(gene, fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Mm.eg.db)

geneList.df <- bitr(geneList, fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = org.Mm.eg.db)

head(gene.df$ENTREZID)

# GO over-representation test
ego <- enrichGO(gene = gene.df$ENTREZID,
               universe = geneList.df$ENTREZID,
               OrgDb = org.Mm.eg.db,
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
                 organism = 'mmu',
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff = 0.05)

#save results
write.csv(kk,file=KEGG_save_path,quote=F,row.names = F)
