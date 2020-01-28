source('http://bioconductor.org/biocLite.R')

biocLite('clusterProfiler')#富集分析
biocLite('topGO')#go分析
biocLite('Rgraphviz')#可视化作图
biocLite('pathview')#kegg pathway分析
biocLite('org.Hs.eg.db')#人的数据库注释

library('clusterProfiler')
library('topGO')
library('Rgraphviz')
library('pathview')
library('org.Hs.eg.db')

#load selected_gens
load("~/Documents/lab/R/goAna/GEO0630.RData")

#pick out gene symbol
DEGs=selected_gens$Gene.symbol

#translate official gene symbol to entrezid
DEGs_entrez_id<-mapIds(x=org.Hs.eg.db,column = 'ENTREZID',
                       keys = as.character(DEGs),keytype = 'SYMBOL')

enrich_go_BP<-enrichGO(gene = DEGs_entrez_id,
                       OrgDb = 'org.Hs.eg.db',
                       ont = 'BP',pvalueCutoff = 0.05)

enrich_go_MF<-enrichGO(gene = DEGs_entrez_id,
                       OrgDb = 'org.Hs.eg.db',
                       ont = 'MF',pvalueCutoff = 0.05)


enrich_go_CC<-enrichGO(gene = DEGs_entrez_id,
                       OrgDb = 'org.Hs.eg.db',
                       ont = 'CC',pvalueCutoff = 0.05)


enrich_go_KEGG<-enrichKEGG(gene = DEGs_entrez_id,organism = 'hsa',
                           keyType = 'kegg',pvalueCutoff = 0.05)



barplot(enrich_go_BP)
barplot(enrich_go_CC)
barplot(enrich_go_MF)
barplot(enrich_go_KEGG)

dotplot(enrich_go_BP,title='Biological Process of DEGs')
dotplot(enrich_go_MF,title='Molecular function of DEGs')
dotplot(enrich_go_KEGG,title='KEGG pathway of DEGs')
dotplot(enrich_go_CC,title='Cell Component of DEGs')







###自定义作气泡图

x<-read.csv(file.choose(),stringsAsFactors = F)
#筛选p<0.05
x<-x[x$PValue<0.05,]
x_go=x[,1:5]

xbp=x_go[grep("BP",x_go$Category),]
xcc=x_go[grep("CC",x_go$Category),]
xmf=x_go[grep("MF",x_go$Category),]
xkegg=x_go[grep("KEGG",x_go$Category),]

xbp$Term=gsub(".*\\~","",xbp$Term)#Biological Process
xcc$Term=gsub(".*\\~","",xcc$Term)#Cell Component
xmf$Term=gsub(".*\\~","",xmf$Term)#Molecular Function
xkegg$Term=gsub(".*\\:","",xkegg$Term)#KEGG pathway

#加载ggplot2
library(ggplot2)
make_GO_bubble<-function(go_data,term_name){
  
  #选择top10的数据（count）
  GO_DATA=go_data[order(go_data$Count,decreasing = T),]
  GO_DATA=head(GO_DATA,10)
  
  # 四维数据的展示
  p = ggplot(GO_DATA,aes(X.,Term))
  bubble=p+ geom_point(aes(size=Count,color=-log10(PValue)))
  
  # 自定义渐变颜色
  bubble =bubble+ scale_colour_gradient(low="green",high="red")
  
  # 改变图片的样式（主题）
  pr=bubble + theme_test(base_size = 16,base_rect_size = 1)
  
  pr=pr+labs(x="Rich factor",y=term_name,title="Enrichment of DEGs")
  
  return(pr) 
}  


#BP
make_GO_bubble(xbp,term_name="Biological Process")#HEIGHT 550
make_GO_bubble(xcc,term_name = "Cell Component")
make_GO_bubble(xmf,term_name = "Molecular Function")
make_GO_bubble(xkegg,term_name = "KEGG pathway")


##DAVID+R