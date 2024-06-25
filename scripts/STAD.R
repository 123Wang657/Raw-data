
if (T) {
  dir.create("scripts")
  dir.create("results")
  dir.create("files")
  dir.create("figures")
  dir.create("origin_datas/GEO",recursive = T)
  dir.create("origin_datas/TCGA")
}

library(stringr)
library(tidydr)
library(openxlsx)
library(data.table)
library(reshape2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(clusterProfiler)
library(pheatmap)
library(ComplexHeatmap)
library(GSVA)
library(GSEABase)
library(fgsea)
library(corrplot)
library(colorspace)
library(survival)
library(survminer)
library(maftools)
library(vegan)
library(forcats)
library(ggpubr)
library(ggsci)
library(ggplot2)
library(rstatix)
library(ggstatsplot)
library(ggcor)
library(ggstance)
options(stringsAsFactors = F)
source('/home/pub252/projects/codes/mg_base.R')
bioForest=function(rt=null,col){
  #
  gene <- rownames(rt)
  hr <- sprintf("%.3f",rt[,2])
  hrLow  <- sprintf("%.3f",rt[,3])
  hrHigh <- sprintf("%.3f",rt[,4])
  Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
  pVal <- ifelse(rt[,1]<0.001, "<0.001", sprintf("%.3f", rt[,1]))
  
  #
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3,2.5))
  
  #
  xlim = c(0,3)
  par(mar=c(4,2.5,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
  text(3.1,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3.1,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1)
  
  #
  par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
  xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
  abline(v=1,col="black",lty=2,lwd=2)
  boxcolor = ifelse(as.numeric(hr) > 1, col[2], col[1])
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.5)
  axis(1)
}
my_volcano=function(dat,p_cutoff=0.05,fc_cutoff=1,col=c("red","blue","black"),
                    ylab='-log10 (adj.PVal)',xlab='log2 (FoldChange)',leg.pos='right'){
  degs_dat=dat$DEG
  degs_dat$type=factor(ifelse(degs_dat$adj.P.Val<p_cutoff & abs(degs_dat$logFC) > fc_cutoff, 
                              ifelse(degs_dat$logFC> fc_cutoff ,'Up','Down'),'No Signif'),levels=c('Up','Down','No Signif'))
  p=ggplot(degs_dat,aes(x=logFC,y=-log10(adj.P.Val),color=type))+
    geom_point()+
    scale_color_manual(values=col)+#
    # geom_text_repel(
    #   data = tcga.diff$DEG[tcga.diff$DEG$adj.P.Val<p_fit & abs(tcga.diff$DEG$logFC)>fc_fit,],
    #   #aes(label = Gene),
    #   size = 3,
    #   segment.color = "black", show.legend = FALSE )+#
    theme_bw()+#
    theme(
      legend.title = element_blank(),#
      legend.position = leg.pos,
      text = element_text(family = 'Times')
    )+
    ylab(ylab)+#
    xlab(xlab)+#
    geom_vline(xintercept=c(-fc_cutoff,fc_cutoff),lty=3,col="black",lwd=0.5) +#
    geom_hline(yintercept = -log10(p_cutoff),lty=3,col="black",lwd=0.5)#
  return(p)
}

custom_theme <- function() {
  theme_survminer() %+replace%
    theme(text = element_text(family = 'Times'),panel.grid = element_blank())
}
getGeneFC=function(gene.exp,group,ulab=ulab,dlab=dlab){
  degs_C1_C3=mg_limma_DEG(gene.exp, 
                          group,
                          ulab=ulab,
                          dlab=dlab)
  
  ## 
  degs_C1_C3_sig<-degs_C1_C3$DEG[which(degs_C1_C3$DEG$adj.P.Val <= 1),]
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(stringr)
  degs_C1_C3_sig_gene<-rownames(degs_C1_C3_sig)
  
  degs_gene_entz=bitr(degs_C1_C3_sig_gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
  degs_gene_entz <- dplyr::distinct(degs_gene_entz,SYMBOL,.keep_all=TRUE)
  
  gene_df <- data.frame(logFC=degs_C1_C3_sig$logFC,
                        SYMBOL = rownames(degs_C1_C3_sig))
  gene_df <- merge(gene_df,degs_gene_entz,by="SYMBOL")
  head(gene_df)
  
  geneList<-gene_df$logFC
  names(geneList)=gene_df$ENTREZID 
  head(geneList)
  
  geneList=sort(geneList,decreasing = T)
  head(geneList) ## 
  return(geneList)
}
my_boxplot=function(dat,group,group_cols=ggsci::pal_aaas()(10),#test_method='kruskal.test',
                    fill= "Group",label=c("p.format",'p.signif')[1],notch=F,size=10,
                    xlab='',ylab='',title='',x.size=10,y.size=10,legend.position='top'){
  
  # dat=ARGs.score[tcga.subtype$Samples,'score']
  # group=tcga.subtype$Cluste
  data=data.frame(Group=group,value=dat)
  data=crbind2DataFrame(data)
  data=melt(data)
  data=data[which(!is.na(data[,1])),]
  if(length(names(table(group)))>2){
    test_method='kruskal.test'
  }else{
    test_method='wilcox.test'
  }
  p=ggplot(data, aes(x=Group, y=value,fill=Group)) +
    geom_boxplot(notch = notch) +
    scale_fill_manual(values = group_cols)+   #
    # if(length(names(table(group)))>2){
    #   test_method=''
    # }
    ggpubr::stat_compare_means(aes(group=Group), label = label, method = test_method)+
    labs(x=xlab, y = ylab, fill = fill) +
    theme_bw()+
    theme(legend.position =legend.position,                 #
          plot.title = element_text(hjust = 0.5),
          text = element_text(family = 'Times',size = size),
          axis.text.x = element_text(size = x.size),
          axis.text.y = element_text(size = y.size)) # 
  return(p)
}
my_violin=function(dat,group,group_cols=ggsci::pal_aaas()(10),#test_method='kruskal.test',
                   fill= "Group",label=c("p.format",'p.signif')[1],
                   xlab='',ylab='',title='',x.size=10,y.size=10,legend.position='top'){
  
  # dat = tcga.b.cell$GeneSet,
  # group = tcga.subtype$Cluster
  data=data.frame(Group=group,value=dat)
  data=crbind2DataFrame(data)
  data=melt(data)
  data=data[which(!is.na(data[,1])),]
  if(length(names(table(group)))>2){
    test_method='kruskal.test'
  }else{
    test_method='wilcox.test'
  }
  p=ggplot(data,aes(x=Group, y=value,fill=Group)) +
    geom_violin()+  
    geom_boxplot(width=0.2,position=position_dodge(0.8),outlier.colour = NA,fill="white")+
    scale_fill_manual(values = group_cols)+
    theme_classic(base_size = 20)+labs(x=xlab,y=ylab,title = title)+
    ggpubr::stat_compare_means(aes(group=Group), label = label, method =test_method)+
    theme(legend.position = legend.position,axis.text = element_text(color = 'black'),
          title = element_text(size = 12),text = element_text(family = 'Times'),
          axis.text.x = element_text(size = 12),axis.title.x = element_text(size = 15),
          axis.text.y = element_text(size = 12),axis.title.y = element_text(size = 12))
  return(p)
}

my_mutiboxplot=function(dat,group,group_cols=ggsci::pal_aaas()(10),
                        #test_method=c('t.test','wilcox.test','anova','kruskal.test')[4],
                        bw=T,xlab='',ylab='score',title='',size=10,angle = 45, hjust = 1,
                        legend.position='top',fill='group',notch=F){
  # dat=tcga.est[tcga.subtype.cli$Samples,]
  # group=tcga.subtype.cli$Cluster
  dat.bind=cbind(dat,Cluster=group)
  dat.bind=crbind2DataFrame(dat.bind)
  dat.melt=melt(dat.bind)
  #data=data[which(!is.na(data[,1])),]
  colnames(dat.melt)=c('Group','type','value')
  if(length(names(table(group)))>2){
    test_method='kruskal.test'
  }else{
    test_method=c('wilcox.test','t.test')[1]
  }
  p=dat.melt %>%
    ggplot(aes(x=type, y=value,fill=Group)) +
    geom_boxplot(notch = notch) +  
    scale_fill_manual(values =group_cols)+   #
    ggpubr::stat_compare_means(aes(group=Group), label = "p.signif", method = test_method)+
    labs(x="", y = ylab, fill =fill,title =title) +
    #theme_light()+
    # theme_bw()+
    theme_classic()+
    theme(legend.position = legend.position,                 #
          plot.title = element_text(hjust = 0.5),
          text = element_text(family = 'Times',size = size),
          axis.text.x = element_text(size = size,angle = angle, hjust = hjust)) # 
  return(p)
}
my_mutiviolin=function(dat,group,group_cols=ggsci::pal_aaas()(10),
                       #test_method=c('t.test','wilcox.test','anova','kruskal.test')[4],
                       bw=T,xlab='',ylab='score',title='',size=10,angle = 45, hjust = 1,
                       legend.position='top',fill='group'){
  # dat=tcga.est[tcga.subtype.cli$Samples,]
  # group=tcga.subtype.cli$Cluster
  dat.bind=cbind(dat,Cluster=group)
  dat.bind=crbind2DataFrame(dat.bind)
  dat.melt=melt(dat.bind)
  #data=data[which(!is.na(data[,1])),]
  colnames(dat.melt)=c('Group','type','value')
  if(length(names(table(group)))>2){
    test_method='kruskal.test'
  }else{
    test_method='wilcox.test'
  }
  p=ggviolin(dat.melt, x = "type", y = "value", fill = "Group",
             add = "boxplot",palette = group_cols)+
    #scale_fill_manual(values =group_cols)+   #
    ggpubr::stat_compare_means(aes(group=Group), label = "p.signif", method = test_method)+
    labs(x="", y = ylab, fill =fill,title =title) +
    theme(text = element_text(family = 'Times',size = size),
          axis.text.x = element_text(size = size,angle = angle, hjust = hjust))
  return(p)
}
my_riskplot=function(cli_dat,cols=c("red","blue"),xlab='Samples',
                     a.ylab="Risk score",b.labs="Survival time(year)",cutoff=0,labs=c('A','B')){
  # cli_dat=tcga.risktype.cli
  cli.dat.order=cli_dat[order(cli_dat$Riskscore),c('OS.time','Status','Riskscore','Risktype')]
  fp_dat=data.frame(Samples=1:nrow(cli_dat),cli.dat.order)
  p1=ggplot(fp_dat,aes(x=Samples,y=Riskscore))+geom_point(aes(color=Risktype))+
    scale_colour_manual(values =cols)+
    theme_bw()+labs(x=xlab,y=a.ylab)+
    geom_hline(yintercept=cutoff,colour="black", linetype="dotted",size=0.8)+
    geom_vline(xintercept=sum(fp_dat$Risktype=="Low"),colour="black", linetype="dotted",size=0.8)
  p1
  p2=ggplot(fp_dat,aes(x=Samples,y=OS.time))+geom_point(aes(col=Status))+theme_bw()+
    scale_colour_manual(values =cols)+
    labs(x=xlab,y=b.labs)+
    geom_vline(xintercept=sum(fp_dat$Risktype=="Low"),colour="black", linetype="dotted",size=0.8)
  p2
  p.merge=mg_merge_plot(p1,p2,nrow=2,ncol=1,labels = labs)
  return(p.merge)
}
coxFun <- function(dat){
  library(survival)
  colnames(dat)=c('time','status','gene')
  fmla=as.formula("Surv(time,status)~gene")
  cox=coxph(fmla,data=dat)
  p=summary(cox)[[7]][5]
  result=c(p,summary(cox)[[8]][1],summary(cox)[[8]][3],summary(cox)[[8]][4])
  return(result)
}


#00.data.pre####
#
dir.create('results/00.data.pre')
genecode=read.delim('00.data.pre/GeneTag.genecode.v32.txt',sep='\t',header = T)
table(genecode$TYPE)
mrna_genecode=genecode[which(genecode$TYPE=='protein_coding'),]
head(mrna_genecode)


#######
tcga.pancancer.cli=read.xlsx('00.data.pre/TCGA/TCGA_pancancer_cli_PMID_29625055.xlsx')
head(tcga.pancancer.cli)
tcga.cli=tcga.pancancer.cli[which(tcga.pancancer.cli$type=='STAD'),]
head(tcga.cli)
tcga.cli=data.frame(Samples=paste0(tcga.cli$bcr_patient_barcode,'-01'),
                    Age=tcga.cli$age_at_initial_pathologic_diagnosis,
                    Gender=tcga.cli$gender,
                    AJCC_stage=tcga.cli$ajcc_pathologic_tumor_stage,
                    # Clinical_Stage=tcga.cli$clinical_stage,
                    Grade=tcga.cli$histological_grade,
                    tcga.cli[,c('OS','OS.time','DSS','DSS.time','DFI','DFI.time','PFI','PFI.time')])
rownames(tcga.cli)=tcga.cli$Samples
head(tcga.cli)
tcga.cli$OS.time
tcga.cli=tcga.cli %>% drop_na(OS.time)
tcga.cli=tcga.cli[tcga.cli$OS.time>0,]
dim(tcga.cli)

head(tcga.cli)
fivenum(tcga.cli$Age)
tcga.cli$Age1=ifelse(tcga.cli$Age>67,'>67','<=67')
table(tcga.cli$Gender)
table(tcga.cli$AJCC_stage)
tcga.cli$AJCC_stage[tcga.cli$AJCC_stage=='[Discrepancy]'|tcga.cli$AJCC_stage=='[Not Available]']=NA
tcga.cli$AJCC_stage=gsub('[ABC]','',tcga.cli$AJCC_stage)
tcga.cli$AJCC_stage=gsub('Stage ','',tcga.cli$AJCC_stage)

table(tcga.cli$Grade)
tcga.cli$Grade[tcga.cli$Grade=='GX']=NA


tcga.data=read.delim('00.data.pre/TCGA/Merge_RNA_seq_TPM.txt',row.names = 1,check.names = F)
tcga.data[1:4,1:4]
table(substr(colnames(tcga.data),14,15))
dim(tcga.data)

sample_T=colnames(tcga.data)[which(as.numeric(substr(colnames(tcga.data),14,15))==1)]#
sample_T=intersect(sample_T,tcga.cli$Samples)
sample_N=colnames(tcga.data)[which(as.numeric(substr(colnames(tcga.data),14,15))==11)]#
tcga_type=data.frame(Samples=c(sample_T,sample_N),Type=rep(c('Tumor','Normal'),c(length(sample_T),length(sample_N))))
rownames(tcga_type)=tcga_type$Samples
table(tcga_type$Type)


length(sample_T)
#353

range(tcga.data)
tcga.data=log2(tcga.data+1)
range(tcga.data)
tcga.data[1:5,1:5]


tcga.exp=tcga.data[,sample_T]
tcga.cli=tcga.cli[sample_T,]



##GSE66229#######
GSE66229=readRDS('results/00.data.pre/GSE66229.rds')

GSE66229.pheno=pData(GSE66229)
GSE66229.pheno=data.frame(Samples=GSE66229.pheno$geo_accession,tissue=GSE66229.pheno$`tissue:ch1`)
rownames(GSE66229.pheno)=GSE66229.pheno$Samples
head(GSE66229.pheno)

GSE66229.df=exprs(GSE66229)
GSE66229.df[1:5,1:5]
GSE66229.exp=exp_probe2symbol_v2(datExpr = GSE66229.df,GPL = 'GPL570')
rm(GSE66229.df)
range(GSE66229.exp)
dim(GSE66229.exp)

GSE66229.cli=read.xlsx('results/00.data.pre/GSE66229_cli.xlsx')
head(GSE66229.cli)
GSE66229.cli=data.frame(Samples=GSE66229.cli$GEO_ID,
                        Age=GSE66229.cli$age,Gender=GSE66229.cli$sex,
                        Stage=GSE66229.cli$Stage,
                        OS=GSE66229.cli$Death,OS.time=GSE66229.cli$OS.m/12*365)
rownames(GSE66229.cli)=GSE66229.cli$Samples
GSE66229.cli$OS.time
dim(GSE66229.cli)






#01.##########
dir.create('results/01.scRNA')
library(Seurat)
library(dplyr)
library(ggplot2)
library(magrittr)
library(gtools)
library(stringr)
library(Matrix)
library(tidyverse)
library(patchwork)
library(data.table)
library(RColorBrewer)
library(ggpubr)
library(dplyr)
library(harmony)
dir_name=list.files('origin_datas/GEO/GSE167297_RAW/')
datalist=list()
for (i in 1:length(dir_name)){
  files = paste0("origin_datas/GEO/GSE167297_RAW/",dir_name[i])
  counts=fread(file = files,sep = '\t',check.names = F)
  counts[1:5,1:5]
  counts=as.data.frame(counts)
  rownames(counts)=counts[,1]
  counts=counts[,-1]
  rownames(counts) <- gsub("_","-", rownames(counts))
  Samples=stringr::str_split_fixed(dir_name[i],'_',4)[,1]
  patient=stringr::str_split_fixed(dir_name[i],'_',4)[,2]
  tissue=stringr::str_split_fixed(dir_name[i],'_',4)[,3]
  datalist[[i]]<- CreateSeuratObject(counts=counts,project = Samples,min.cells = 3, min.features = 200)
  datalist[[i]] <- AddMetaData(datalist[[i]] , Samples,col.name = "Samples")
  datalist[[i]] <- AddMetaData(datalist[[i]] , patient,col.name = "patient")
  datalist[[i]] <- AddMetaData(datalist[[i]] , tissue,col.name = "tissue")
}
names(datalist)=stringr::str_split_fixed(dir_name,'_',4)[,1]


####
for (i in 1:length(datalist)){
  sce <- datalist[[i]]
  sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")# 
  sce[["percent.Ribo"]] <- PercentageFeatureSet(sce, pattern = "^RP[SL]")# 
  datalist[[i]] <- sce
  rm(sce)
}
sce <- merge(datalist[[1]],y=datalist[2:length(datalist)])
#
raw_meta=sce@meta.data
raw_count <- table(raw_meta$Samples)
raw_count
sum(raw_count)#  23531
pearplot_befor<-VlnPlot(sce,group.by ='Samples',
                        features = c("nFeature_RNA", "nCount_RNA","percent.mt"),
                        pt.size = 0,
                        ncol = 3)
pearplot_befor
ggsave('results/01.scRNA/pearplot_befor.pdf',pearplot_befor,height = 5,width = 10)

#
sce=subset(sce, subset=nFeature_RNA>200 & nFeature_RNA<2000 & percent.mt<15)
rm(datalist)
clean_meta=sce@meta.data
clean_count <- table(clean_meta$Samples)
clean_count
sum(clean_count)#17397
pearplot_after <- VlnPlot(sce,group.by ='Samples',
                          features = c("nFeature_RNA", "nCount_RNA","percent.mt"),
                          pt.size = 0,
                          ncol = 3)
pearplot_after
ggsave('results/01.scRNA/pearplot_after.pdf',pearplot_after,height = 6,width = 15)

# 
sce <- NormalizeData(sce)
sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
sce <- ScaleData(sce, features = rownames(sce))
#
sce <- RunPCA(sce, features = VariableFeatures(sce))
colnames(sce@meta.data)
##
library(harmony)
sce = RunHarmony(sce, group.by.vars="Samples")
#
pca.plot=ElbowPlot(sce,ndims = 50)+theme(text = element_text(family = 'Times',size = 12))
pca.plot
ggsave('results/01.scRNA/PCA.pdf',pca.plot,height = 5,width = 5)
###
sce <- RunTSNE(sce, dims=1:10, reduction="harmony")
library(clustree)
sce <- FindNeighbors(sce, dims = 1:10, reduction="harmony")
#
sce <- FindClusters(object = sce,resolution = .2)
DefaultAssay(sce) <- "RNA"
colnames(sce@meta.data)
length(table(sce@meta.data$seurat_clusters))
sce=subset(sce,seurat_clusters%in%c(0:8))
####
after_batch=DimPlot(sce,group.by='Samples',reduction="tsne",label = F,pt.size = 0.2)+
  theme_dr(xlength = 0.3, ylength = 0.3,arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(text = element_text(family = 'Times',size = 12),panel.grid = element_blank())
after_batch
ggsave('results/01.scRNA/after_batch.pdf',after_batch,height = 6,width = 6.5)

p=DimPlot(sce,group.by='seurat_clusters',reduction="tsne",label = F,pt.size = 0.2)+
  theme_dr(xlength = 0.3, ylength = 0.3,arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(text = element_text(family = 'Times',size = 12),panel.grid = element_blank())
p=LabelClusters(p,id = 'seurat_clusters',family='Times')
ggsave('results/01.scRNA/cluster_tsne.pdf',p,height = 5,width = 5.5)


# #
# Logfc = 0.5
# #
# Minpct = 0.25
# DefaultAssay(sce) <- "RNA"
# colnames(sce@meta.data)
# Idents(sce)<-'seurat_clusters'
# 
# sce.markers <- FindAllMarkers(object = sce,logfc.threshold = Logfc, min.pct = Minpct,only.pos = T)
# head(sce.markers)
# sce.markers["pct.diff"]=sce.markers$pct.1-sce.markers$pct.2
# sce.markers <- sce.markers[sce.markers$p_val_adj<0.05,]
# table(sce.markers$cluster)
# length(unique(sce.markers$gene))#762
# head(sce.markers)
# table(sce.markers$cluster)
# write.csv(sce.markers,'results/01.scRNA/cluster_diff_marker_gene.csv')
# ### 
# Top5 <- sce.markers %>% group_by(cluster) %>% slice_max(n =5, order_by = avg_logFC)
# length(Top5$gene)
# length(unique(Top5$gene))
# ###
# diff.marker.dotplot=DotPlot(object = sce, features = unique(Top5$gene),
#                             cols=c("snow", "blue"),scale = T)+
#   RotatedAxis()+ ggtitle("Marker Genes")+
#   theme(plot.title = element_text(hjust = 0.5)) +
#   xlab('')+ylab('')+coord_flip()
# diff.marker.dotplot
# ggsave('results/01.scRNA/cluster_diff_marker_gene.pdf',diff.marker.dotplot,height = 10,width = 10)
# 
# VlnPlot(sce,features = c('CD14','CD68','CD11c'),pt.size = 0,group.by = 'seurat_clusters')

##1.1 ######
marker <- data.frame(cluster = 0:8,cell = 0:8)
marker[marker$cluster %in% c(0,2),2] <- 'T cells'
marker[marker$cluster %in% c(1),2] <- 'Dendritic cells'
marker[marker$cluster %in% c(3),2] <- 'Macrophage'
marker[marker$cluster %in% c(4),2] <- 'B cells'
marker[marker$cluster %in% c(5),2] <- 'Endothelial cells'
marker[marker$cluster %in% c(6),2] <- 'Epithelial cells'
marker[marker$cluster %in% c(7),2] <- 'Fibroblast'
marker[marker$cluster %in% c(8),2] <- 'Mast cells'
marker
sce@meta.data$cell_type <- sapply(sce@meta.data$seurat_clusters,function(x){marker[x,2]})
cell_type_umap=DimPlot(sce,group.by='cell_type',reduction="tsne",label = F,pt.size = 0.1,cols =pal_d3()(8))+
  theme_dr(xlength = 0.3, ylength = 0.3,arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(panel.grid = element_blank(),text = element_text(family = 'Times',size = 12))
cell_type_umap
ggsave('results/01.scRNA/cell_type_tsne.pdf',cell_type_umap,height = 5,width = 5.5)

marker_gene=c('MZB1','IGLL5','IGJ',
              'CD2','CD3D','CD3E','NKG7','GZMA',
              'KRT8','KRT18','KRT19',
              'IL8','CD68','FCER1G',
              'HLA-DQB1','HLA-DRA',
              'CD123','CD83',
              'PLVAP','TM4SF1',
              'APOD','DCN','ACTA2',
              'TPSAB1','CPA3','CTSG'
              )
Idents(sce)='cell_type'
marker.dot=DotPlot(object = sce, features = marker_gene,
        cols=c("snow", "red"),scale = T)+
  RotatedAxis()+ ggtitle("Marker Genes")+
  theme(plot.title = element_text(hjust = 0.5),axis.text.y = element_text(size=12),text = element_text(family = 'Times',size = 12)) +
  xlab('')+ylab('')+coord_flip()
ggsave('results/01.scRNA/marker_dot.pdf',marker.dot,height = 6,width = 6)



cell_freq1=data.frame(t(prop.table(table(sce$cell_type,sce$Samples),margin=2)))
cell_freq1
colnames(cell_freq1)<-c('Samples','cell_type','Freq')
cell_prop_fig1=ggplot(cell_freq1,aes(x=Samples,y=Freq,fill=cell_type))+
  scale_fill_manual(values = pal_d3()(9))+
  # facet_grid(~Tissue,scales = 'free',space='free')+
  geom_bar(position = "fill",stat="identity")+
  xlab('')+ylab('Proportion')+
  theme(text = element_text(family = 'Times',size=12),
        axis.text.x = element_text(angle = 30,hjust = 1),
        legend.text =element_text(family = 'Times'),
        legend.title = element_text(family = 'Times'))
cell_prop_fig1

###########
#
Logfc = 0.5
#
Minpct = 0.25
DefaultAssay(sce) <- "RNA"
colnames(sce@meta.data)
Idents(sce)<-'cell_type'

sce.markers <- FindAllMarkers(object = sce,logfc.threshold = Logfc, min.pct = Minpct,only.pos = T)
head(sce.markers)
sce.markers["pct.diff"]=sce.markers$pct.1-sce.markers$pct.2
sce.markers <- sce.markers[sce.markers$p_val_adj<0.05,]
table(sce.markers$cluster)
length(unique(sce.markers$gene))#762
head(sce.markers)
table(sce.markers$cluster)
write.csv(sce.markers,'results/01.scRNA/cell_type_DEGs.csv')
write.csv(sce.markers[marker_gene,],'results/01.scRNA/cell_type_marker_genes.csv')

##1.2 ####
hp.geneset1=read.gmt('origin_datas/HP_HELICOBACTER_PYLORI_INFECTION.v2023.2.Hs.gmt')
hp.geneset2=read.gmt('origin_datas/KEGG_EPITHELIAL_CELL_SIGNALING_IN_HELICOBACTER_PYLORI_INFECTION.v2023.2.Hs.gmt')

hp.genes=unique(c(hp.geneset1$gene,hp.geneset2$gene))
length(hp.genes)
#73

library(AUCell)
library(GSEABase)
countexp<-sce@assays$RNA@counts
countexp<-data.frame(as.matrix(countexp))
cells_rankings <- AUCell_buildRankings(as.matrix(countexp)) #rank
geneSets <- list(HP=hp.genes) #signature read
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings) #calc

HP_AUCell <- as.numeric(getAUC(cells_AUC)['HP', ])
sce$HP_AUCell  <- HP_AUCell

sc_df=data.frame(sce@meta.data, sce@reductions$tsne@cell.embeddings)
head(sc_df)


hp.boxplot=ggplot(sc_df) +
  geom_boxplot(aes(x = reorder(cell_type, HP_AUCell, FUN = median), y = HP_AUCell,fill=cell_type))+
  scale_fill_manual(values = pal_d3()(10))+
  xlab('')+ylab('HP AUCell')+ theme_classic()+
  theme(text = element_text(family = 'Times',size = 12),legend.position = 'none',
        axis.text.x = element_text(color = "black", size = 12,angle = 30,hjust = 1))
hp.boxplot
  
fig1=mg_merge_plot(mg_merge_plot(p,cell_type_umap,marker.dot,cell_prop_fig1,ncol=2,nrow=2,labels = LETTERS[1:4],heights = c(1,1.5)),
                   hp.boxplot,nrow=2,labels = c('','E'),heights = c(2,1))

ggsave('figures/Fig1.pdf',fig1,height = 15,width = 10)


##1.3 ###########
dir.create('results/01.scRNA/Macrophage')
table(sce$cell_type)
immuce_cells=subset(sce,cell_type %in% 'Macrophage')

immuce_cells <- NormalizeData(immuce_cells)
immuce_cells <- FindVariableFeatures(immuce_cells, selection.method = "vst", nfeatures = 2000)
immuce_cells <- ScaleData(immuce_cells, features = rownames(immuce_cells))
immuce_cells <- RunPCA(immuce_cells, features = VariableFeatures(immuce_cells))
immuce_cells = RunHarmony(immuce_cells, group.by.vars="Samples")
ElbowPlot(immuce_cells,ndims = 50)+theme(text = element_text(family = 'Times',size = 12))
dev.off()
immuce_cells <- RunTSNE(immuce_cells, dims=1:10, reduction="harmony")
immuce_cells <- FindNeighbors(immuce_cells, dims = 1:10, reduction="harmony")
immuce_cells <- FindClusters(object = immuce_cells,resolution = .4)
DefaultAssay(immuce_cells) <- "RNA"
length(table(immuce_cells@meta.data$seurat_clusters))
immuce_cells=subset(immuce_cells,seurat_clusters%in%c(0:4,6),)

marker2 <- data.frame(cluster = 0:7,cell = 0:7)
marker2[marker2$cluster %in% c(0),2] <- 'Macrophage C1'
marker2[marker2$cluster %in% c(1),2] <- 'Macrophage C2'
marker2[marker2$cluster %in% c(2),2] <- 'Macrophage C3'
marker2[marker2$cluster %in% c(3),2] <- 'Macrophage C4'
marker2[marker2$cluster %in% c(4),2] <- 'Macrophage C5'
marker2[marker2$cluster %in% c(6),2] <- 'Macrophage C6'
marker2
immuce_cells@meta.data$cell_type <- sapply(immuce_cells@meta.data$seurat_clusters,function(x){marker2[x,2]})
Idents(immuce_cells)='cell_type'
p1=DimPlot(immuce_cells,group.by='cell_type',reduction="tsne",label = F,pt.size = 0.2,cols =pal_jama()(8))+
  theme_dr(xlength = 0.3, ylength = 0.3,arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(panel.grid = element_blank(),text = element_text(family = 'Times',size = 12))
ggsave('results/01.scRNA/Macrophage/TAM_cluster_tsne.pdf',p1,height = 5,width = 5.5)

#########
sc_df2=data.frame(immuce_cells@meta.data, immuce_cells@reductions$tsne@cell.embeddings)
head(sc_df2)
table(sc_df2$cell_type)
p2=ggplot(sc_df2) +
  geom_boxplot(aes(x = reorder(cell_type, HP_AUCell, FUN = median), y = HP_AUCell,fill=cell_type))+
  scale_fill_manual(values = pal_jama()(10))+
  xlab('')+ylab('HP AUCell')+theme_classic()+
  theme(text = element_text(family = 'Times',size = 12),legend.position = 'none',
        axis.text.x = element_text(color = "black", size = 12,angle = 30,hjust = 1))
ggsave('results/01.scRNA/Macrophage/TAM_HP_AUCell.pdf',p2,height = 5,width = 6)

saveRDS(immuce_cells,file='results/01.scRNA/Macrophage.rds')


########
#
Logfc = 0.5
#
Minpct = 0.25
DefaultAssay(immuce_cells) <- "RNA"
colnames(immuce_cells@meta.data)
Idents(immuce_cells)<-'cell_type'

sce.markers <- FindAllMarkers(object = immuce_cells,logfc.threshold = Logfc, min.pct = Minpct,only.pos = T)
head(sce.markers)
sce.markers["pct.diff"]=sce.markers$pct.1-sce.markers$pct.2
sce.markers <- sce.markers[sce.markers$p_val_adj<0.05,]
table(sce.markers$cluster)
length(unique(sce.markers$gene))#762
head(sce.markers)
table(sce.markers$cluster)
write.csv(sce.markers,'results/01.scRNA/Macrophage/TAM_cluster_diff_marker_gene.csv')

##
# Top5 <- sce.markers %>%group_by(cluster)   %>% slice_max(n =10, order_by = avg_logFC)
# Top5=Top5[Top5$cluster=='Macrophage C4','gene']
p3=DotPlot(object = immuce_cells, features = c('APOE','C1QC','C1QB','C1QA','CCL2','FCN1','CCL19','CCL22',
                                               'SPP1','CCL3','CXCL5','CCL3L3','HLA-DQA1','FCER1A'),
           cols=c("snow", "blue"),scale = T)+
  RotatedAxis()+ ggtitle("Marker Genes")+
  theme(plot.title = element_text(hjust = 0.5),text = element_text(family = 'Times')) +
  xlab('')+ylab('')+coord_flip()
ggsave('results/01.scRNA/Macrophage/TAM_marker_dot.pdf',p3,height = 5,width = 5.5)



###C#####
c4.enrichment=mg_clusterProfiler(sce.markers$gene[sce.markers$cluster=='Macrophage C4'])
p4=enrichplot::dotplot(c4.enrichment$KEGG)+ggtitle('KEGG')+theme(text = element_text(family = 'Times'))+
  scale_y_discrete(labels=function(y)str_wrap(y,width = 25))
# enrichplot::dotplot(c4.enrichment$GO_BP)+ggtitle('Biological Process')
# enrichplot::dotplot(c4.enrichment$GO_CC)+ggtitle('Cell Component')
# enrichplot::dotplot(c4.enrichment$GO_MF)+ggtitle('Molecular Function')
write.xlsx(c4.enrichment$KEGG@result,'results/01.scRNA/Macrophage/TAM_cluster_enrichment.xlsx',overwrite = T)
fig2=mg_merge_plot(p1,p2,p3,p4,ncol=2,nrow=2,labels = LETTERS[1:4],widths = c(1,1.2))
ggsave('figures/Fig2.pdf',fig2,height = 10,width = 12)

##1.4 ########
dir.create('results/01.scRNA/CellChat')
library(CellChat)
sce$cell_type2=''
sce$cell_type2=sce$cell_type
sce$cell_type2[rownames(sc_df2)]=sc_df2$cell_type
table(sce$cell_type2)
sce=subset(sce, subset=(cell_type2!='Macrophage'))
table(sce$cell_type2)

library(CellChat)
colnames(sce@meta.data)
table(sce$cell_type2)
# #
DefaultAssay(sce) <- 'RNA'
cellchat <- createCellChat(object = sce, meta = sce@meta.data, group.by = "cell_type2")
cellchat@DB <- subsetDB(CellChatDB = CellChatDB.human, search = "Secreted Signaling")
cellchat <- subsetData(cellchat)
cellchat@data.signaling
#
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
# # #
# #
# # #
cellchat <- computeCommunProb(cellchat)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat)
# # 
# # 
 library(NMF)
# selectK(cellchat, pattern = "outgoing")
# selectK(cellchat, pattern = "incoming")

cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = 3)
# save(cellchat,file='results/01.scRNA/CellChat/cellchat.RData')
# load('results/01.scRNA/CellChat/cellchat.RData')

groupSize <- as.numeric(table(cellchat@idents))

cellchat@netP$pathways
length(cellchat@netP$pathways)
#17


netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing") 
netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming") 

par(mfrow=c(1,2))
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F,
                 title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F,
                 title.name = "Interaction weights/strength")
dev.off()



netAnalysis_dot(cellchat, pattern = "outgoing") 


table(sce$cell_type2)
pdf('results/01.scRNA/CellChat/cellchat_dotplot_1.pdf',height = 5,width = 8)
netVisual_bubble(cellchat, signaling = cellchat@netP$pathways,sources.use = c('Macrophage C4'))
dev.off()

pdf('results/01.scRNA/CellChat/cellchat_dotplot_2.pdf',height = 5,width = 8)
netVisual_bubble(cellchat, signaling = cellchat@netP$pathways,targets.use = c('Macrophage C4'))
dev.off()


##1.5 #####
load('results/tcga.cellscore.RData')
tcga.cell.ssgsea[1:5,1:5]
tcga.cli.sc.merge=cbind.data.frame(tcga.cli,t(tcga.cell.ssgsea[,tcga.cli$Samples]))
tcga.cli.sc.merge=crbind2DataFrame(tcga.cli.sc.merge)
head(tcga.cli.sc.merge)
summary(coxph(formula=Surv(OS.time, OS)~`Macrophage C4` ,data=tcga.cli.sc.merge))
# tcga.cli.sc.merge$group=ifelse(tcga.cli.sc.merge$`Macrophage C4`>median(tcga.cli.sc.merge$`Macrophage C4`),'High','Low')

library(survival)
library(survminer)
cutoff<-surv_cutpoint(tcga.cli.sc.merge,time="OS.time",event="OS",variables='Macrophage C4')
summary(cutoff)
cutoff$cutpoint$cutpoint
tcga.cli.sc.merge$group=ifelse(tcga.cli.sc.merge$`Macrophage C4`> cutoff$cutpoint$cutpoint,'High','Low')
pdf('results/02.WGCNA/KM.pdf',height = 5,width = 6,onefile = F)
ggsurvplot(fit=survfit(Surv(OS.time/365, OS) ~ group,data = tcga.cli.sc.merge),
           data=tcga.cli.sc.merge,ggtheme = custom_theme(),
           conf.int = F,pval = T,fun = "pct",risk.table =F, size = 0.7,
           surv.median.line = 'hv',linetype = c("solid", "dashed","strata")[1],
           palette = pal_d3()(10),legend = c(0.8,0.75), legend.title = "")+ggtitle('Macrophage C4')
dev.off()

#02.WGCNA####
dir.create('results/02.WGCNA')
library(WGCNA)
allowWGCNAThreads(nThreads = 36)#
enableWGCNAThreads(nThreads = 36)# 


table(tcga_type$Type)
tcga.limma=mg_limma_DEG(tcga.data[intersect(mrna_genecode$SYMBOL,rownames(tcga.data)),tcga_type$Samples],
                        group = tcga_type$Type,ulab = 'Tumor',dlab = 'Normal')
tcga.limma$Summary
#"6377|969"
tcga.degs=tcga.limma$DEG[tcga.limma$DEG$adj.P.Val<0.05&abs(tcga.limma$DEG$logFC)>log2(1.5),]
dim(tcga.degs)


tcga.exp.use=tcga.exp[intersect(mrna_genecode$SYMBOL,rownames(tcga.degs)),]
tpm_T2=tcga.exp.use
dim(tpm_T2)
tpm_T2=t(tpm_T2)
dim(tpm_T2)
range(tpm_T2)

pdf('results/02.WGCNA/1.pdf',width = 10,height = 10)
tpm_T2.power=mg_wgcna_get_power(tpm_T2)
dev.off()

# minModuleSize = 30,    ##
# mergeCutHeight = 0.25, ##
tpm_T2.power$cutPower
tpm_T2.module=mg_WGCNA_getModule(tpm_T2,
                                 power = tpm_T2.power$cutPower,
                                 deepSplit=2,
                                 mergeCutHeight=0.3,
                                 minModuleSize=50)

table(tpm_T2.module$Modules[,2])
length(table(tpm_T2.module$Modules[,2]))

pdf('results/02.WGCNA/2.pdf',height = 5,width = 6)
plotDendroAndColors(tpm_T2.module$Tree, tpm_T2.module$Modules,
                    c("Dynamic Module",'Merged Module'),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
# 
#writeMatrix(tpm_T2.module$Modules,outpath = 'results/02.WGCNA/tcga.wgcna.module.genes.txt')
pdf('results/02.WGCNA/3.pdf',height = 6,width = 6)
mg_barplot_point(labels = names(table(tpm_T2.module$Modules[,2]))
                 ,values = as.numeric(table(tpm_T2.module$Modules[,2]))
                 ,point_sizes = 2
                 ,point_cols = names(table(tpm_T2.module$Modules[,2]))
                 ,xlab = 'Number of Genes',legend.pos = NULL)
dev.off()

#### 
# Calculate eigengenes
MEs = tpm_T2.module$MEs
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
pdf('results/02.WGCNA/4.pdf',height = 6,width = 12,onefile = T)
plot(METree, main = "Clustering of module eigengenes",xlab = "", sub = "")
dev.off()




##choose 
colnames(tcga.cli.sc.merge)

tcga_cli_use.part=tcga.cli.sc.merge[,c(2:5,18)]
str(tcga_cli_use.part)

tcga_cli_use.part=sapply(tcga_cli_use.part, function(x)as.numeric(as.factor(x)))
spms=tcga_cli_use.part

MEs_col<-tpm_T2.module$MEs
dim(MEs_col)

modTraitCor = cor(MEs_col[,rownames(MEDiss)[METree$order]]
                  , spms
                  ,use = 'pairwise.complete.obs')
modTraitP = corPvalueStudent(modTraitCor, dim(spms)[1])

textMatrix = paste(signif(modTraitCor, 2), " (", format(modTraitP,scientific =TRUE,digits = 3), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)
dim(textMatrix)

rownames(modTraitCor)=gsub("ME","",rownames(modTraitCor))
rownames(textMatrix)=gsub("ME","",rownames(textMatrix))
colnames(modTraitCor)

pdf('results/02.WGCNA/5.pdf',width =10,height = 5)
labeledHeatmap(Matrix = data.frame(modTraitCor),
               xLabels = colnames(modTraitCor),
               yLabels = rownames(modTraitCor),
               cex.lab = 1,
               ySymbols = colnames(t(modTraitCor)), colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = data.frame(textMatrix),
               setStdMargins = FALSE,xLabelsAngle = 0,
               cex.text = 0.8, zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

#
#
## 
## 
geneModuleMembership <- signedKME(tpm_T2
                                  , data.frame(tpm_T2.module$MEs)
                                  , outputColumnName = "")
head(geneModuleMembership)
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership)
                                           , nrow(tpm_T2.module$MEs)))
#
geneTraitSignificance <- as.data.frame(cor(tpm_T2
                                           , spms
                                           , use = 'pairwise.complete.obs'))

GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance)
                                           , nrow(spms)))

modNames<-colnames(geneModuleMembership)
modNames

module = "green"
column = match(module, modNames)
column
moduleGenes <- (tpm_T2.module$Modules[,'mergedColors']==module)
tcga.wgcna.gene=names(which(moduleGenes))
length(tcga.wgcna.gene)
#364

wgcna.enrichment=mg_clusterProfiler(genes = tcga.wgcna.gene)
p1=dotplot(wgcna.enrichment$KEGG)+ggtitle('KEGG')+scale_y_discrete(labels=function(y)str_wrap(y,width = 25))+theme(text=element_text(family = 'Times'))
p2=dotplot(wgcna.enrichment$GO_BP)+ggtitle('Biological Process')+scale_y_discrete(labels=function(y)str_wrap(y,width = 25))+theme(text=element_text(family = 'Times'))
p3=dotplot(wgcna.enrichment$GO_CC)+ggtitle('Cell Component')+scale_y_discrete(labels=function(y)str_wrap(y,width = 25))+theme(text=element_text(family = 'Times'))
p4=dotplot(wgcna.enrichment$GO_MF)+ggtitle('Molecular Function')+scale_y_discrete(labels=function(y)str_wrap(y,width = 25))+theme(text=element_text(family = 'Times'))

pdf('results/02.WGCNA/enrichment.pdf',height = 12,width = 12,onefile = F)
mg_merge_plot(p1,p2,p3,p4,ncol=2,nrow=2,labels = LETTERS[1:4])
dev.off()

write.xlsx(list(KEGG=wgcna.enrichment$KEGG@result,BP=wgcna.enrichment$GO_BP@result,
                CC=wgcna.enrichment$GO_CC@result,MF=wgcna.enrichment$GO_MF@result),
           'results/02.WGCNA/wgcna_enrichment_res.xlsx',overwrite = T)

#03.#######
dir.create('results/03.model')
com.gene=Reduce(intersect,list(tcga.wgcna.gene,rownames(tcga.exp),rownames(GSE66229.exp)))
length(com.gene)

tcga_model_data <- cbind(tcga.cli[, c("OS.time", "OS")],
                         t(tcga.exp[com.gene, tcga.cli$Samples]))
colnames(tcga_model_data) <- gsub('-', '_', colnames(tcga_model_data))

GSE66229_model_data <- cbind(GSE66229.cli[, c("OS.time", "OS")],
                             t(GSE66229.exp[com.gene, GSE66229.cli$Samples]))
colnames(GSE66229_model_data) <- gsub('-', '_', colnames(GSE66229_model_data))


num <- 1098   #OK
tra.samples <- rownames(read.delim(paste0('files/model_select/tra.dat_zscore_',num,'.txt'), 
                                   header = T, row.names = 1, stringsAsFactors = F))
test.samples <- rownames(read.delim(paste0('files/model_select/test.dat_zscore_',num,'.txt'),
                                    header = T, row.names = 1, stringsAsFactors = F))
tra.data <- tcga_model_data[tra.samples, ]
dim(tra.data)
test.data <- tcga_model_data[test.samples, ]
dim(test.data)


tra.cox=cox_batch(dat = tcga.exp[com.gene,tra.samples],
                  time = tcga.cli[tra.samples,]$OS.time,event = tcga.cli[tra.samples,]$OS)
tra.cox=na.omit(tra.cox)
head(tra.cox)
rownames(tra.cox)=gsub('-','__',rownames(tra.cox))
p_cutoff=0.05
table(tra.cox$p.value<p_cutoff)
tra.cox.fit=tra.cox[tra.cox$p.value<p_cutoff,]

#########lasso
library(glmnet)
set.seed(num)
fit1=glmnet(as.matrix(tra.data[,rownames(tra.cox.fit)])
            #,factor(samps)
            ,cbind(time=tra.data$OS.time,
                   status=tra.data$OS)
            ,family="cox"
            #,family="binomial"
            #,type.measure="deviance"
            ,nlambda=100
            , alpha=1) 

cv.fit<-cv.glmnet(as.matrix(tra.data[,rownames(tra.cox.fit)])
                  #,factor(samps)
                  ,cbind(time=tra.data$OS.time,
                         status=tra.data$OS)
                  ,family="cox"
                  ,nlambda=100
                  , alpha=1)
sig.coef <- coefficients(cv.fit,s=cv.fit$lambda.min)[which(coefficients(cv.fit,s=cv.fit$lambda.min)[,1]!=0),1]
cv.fit$lambda.min
names(sig.coef)


par(mfrow=c(1,2))
plot(fit1)
plot(cv.fit)
dev.off()

fmla <- as.formula(paste0("Surv(OS.time, OS) ~",paste0(names(sig.coef),collapse = '+')))
cox <- coxph(fmla, data =as.data.frame(tra.data))
cox=step(cox)

module.coxforest=ggforest(cox, data = tra.data, 
                          main = "Hazardratio", fontsize =1.0, 
                          noDigits = 2)
module.coxforest

lan <- coef(cox)
lan
paste0(round(lan, 3), '*', names(lan),collapse = '+')
##训练集########
risktype.col=c("orange","blue")
risk.tra=as.numeric(lan%*%as.matrix(t(tra.data[tra.samples,names(lan)])))
tra.risktype.cli=data.frame(tcga.cli[tra.samples,],Riskscore=risk.tra)

#######
tra.data.point <- surv_cutpoint(tra.risktype.cli, time = "OS.time", event = "OS",
                                variables = 'Riskscore')
tra.cutoff <- as.numeric(summary(tra.data.point)[1])
tra.cutoff
tra.risktype.cli$Risktype=ifelse(tra.risktype.cli$Riskscore>tra.cutoff,'High','Low')
tra.km=ggplotKMCox(data.frame(tra.risktype.cli$OS.time/365,
                              tra.risktype.cli$OS,
                              tra.risktype.cli$Risktype),
                   palette = risktype.col,show_confint = F,title = 'Train cohort')
tra.km

tra.roc=ggplotTimeROC(tra.risktype.cli$OS.time,
                      tra.risktype.cli$OS,
                      tra.risktype.cli$Riskscore,mks = c(1:5))
tra.roc



################
risk.test=as.numeric(lan%*%as.matrix(t(test.data[test.samples,names(lan)])))
test.risktype.cli=data.frame(tcga.cli[test.samples,],Riskscore=risk.test)

#######
test.data.point <- surv_cutpoint(test.risktype.cli, time = "OS.time", event = "OS",
                                 variables = 'Riskscore')
test.cutoff <- as.numeric(summary(test.data.point)[1])
test.cutoff
test.risktype.cli$Risktype=ifelse(test.risktype.cli$Riskscore>test.cutoff,'High','Low')
test.km=ggplotKMCox(data.frame(test.risktype.cli$OS.time/365,
                               test.risktype.cli$OS,
                               test.risktype.cli$Risktype),
                    palette = risktype.col,show_confint = F,title = 'Test cohort')
test.km


test.roc=ggplotTimeROC(test.risktype.cli$OS.time,
                       test.risktype.cli$OS,
                       test.risktype.cli$Riskscore,mks = c(1:5))
test.roc




###############
risk.tcga=as.numeric(lan%*%as.matrix(t(tcga_model_data[tcga.cli$Samples,names(lan)])))
tcga.risktype.cli=data.frame(tcga.cli,Riskscore=risk.tcga)

#######
tcga.data.point <- surv_cutpoint(tcga.risktype.cli, time = "OS.time", event = "OS",
                                 variables = 'Riskscore')
tcga.cutoff <- as.numeric(summary(tcga.data.point)[1])
tcga.cutoff
tcga.risktype.cli$Risktype=ifelse(tcga.risktype.cli$Riskscore>tcga.cutoff,'High','Low')
tcga.km=ggplotKMCox(data.frame(tcga.risktype.cli$OS.time/365,
                               tcga.risktype.cli$OS,
                               tcga.risktype.cli$Risktype),
                    palette = risktype.col,show_confint = F,title = 'TCGA cohort')
tcga.km

tcga.roc=ggplotTimeROC(tcga.risktype.cli$OS.time,
                       tcga.risktype.cli$OS,
                       tcga.risktype.cli$Riskscore,mks = c(1:5))
tcga.roc


tcga.risktype.cli$Status=ifelse(tcga.risktype.cli$OS==0,'Alive','Death')
pdf('results/03.model/part2_1.pdf',height = 5,width = 5,onefile = F)
my_riskplot(cli_dat = tcga.risktype.cli,cols = c("orange","blue"),xlab = 'sample',
            a.ylab = 'Riskscore',b.labs = 'Time',cutoff = tcga.cutoff)
dev.off()

cli.anno=tcga.risktype.cli[order(tcga.risktype.cli$Riskscore),'Risktype',drop=F]
head(cli.anno)

risktype.col.use=risktype.col
names(risktype.col.use)=c('High','Low')
pdf('results/03.model/part2_2.pdf',height = 3,width = 5,onefile = F)
pheatmap(tcga.exp[c('ABCA1','IKBIP','NPC2','TNFRSF1B','CTLA4','AKAP5'),rownames(cli.anno)],
         scale = 'row',name = 'Expr',
         color = colorRampPalette(c('navy', "white", 'yellow'))(100),
         cluster_rows = F,cluster_cols = F,
         show_rownames = T,show_colnames = F,
         annotation_col = cli.anno,annotation_colors = list(Risktype=risktype.col.use))
dev.off()

pdf('results/03.model/part2_3.pdf',height = 9,width = 5,onefile = F)
mg_merge_plot(tcga.km,tcga.roc,nrow=2)
dev.off()

##GSE66229#########
cox1 <- coxph(as.formula(paste0("Surv(OS.time, OS) ~",paste0(names(lan),collapse = '+'))), 
              data =as.data.frame(GSE66229_model_data))
lan1 <- coef(cox1)
lan1
risk.GSE66229=as.numeric(lan1%*%as.matrix(t(GSE66229_model_data[GSE66229.cli$Samples,names(lan1)])))
GSE66229.risktype.cli=data.frame(GSE66229.cli,Riskscore=risk.GSE66229)

#######
GSE66229.data.point <- surv_cutpoint(GSE66229.risktype.cli, time = "OS.time", event = "OS",
                                     variables = 'Riskscore')
GSE66229.cutoff <- as.numeric(summary(GSE66229.data.point)[1])
GSE66229.cutoff
GSE66229.risktype.cli$Risktype=ifelse(GSE66229.risktype.cli$Riskscore>GSE66229.cutoff,'High','Low')
GSE66229.roc=ggplotTimeROC(GSE66229.risktype.cli$OS.time,
                           GSE66229.risktype.cli$OS,
                           GSE66229.risktype.cli$Riskscore,mks = c(1:5))
GSE66229.roc

GSE66229.km=ggplotKMCox(data.frame(GSE66229.risktype.cli$OS.time/365,
                                   GSE66229.risktype.cli$OS,
                                   GSE66229.risktype.cli$Risktype),
                        palette = risktype.col,show_confint = F,title = 'GSE66229 cohort')
GSE66229.km

GSE66229.risktype.cli$Status=ifelse(GSE66229.risktype.cli$OS==0,'Alive','Death')
pdf('results/03.model/part3_1.pdf',height = 5,width = 5,onefile = F)
my_riskplot(cli_dat = GSE66229.risktype.cli,cols = c("orange","blue"),xlab = 'sample',
            a.ylab = 'Riskscore',b.labs = 'Time',cutoff = GSE66229.cutoff)
dev.off()

GSE66229.cli.anno=GSE66229.risktype.cli[order(GSE66229.risktype.cli$Riskscore),'Risktype',drop=F]
head(GSE66229.cli.anno)
pdf('results/03.model/part3_2.pdf',height = 3,width = 5,onefile = F)
pheatmap(GSE66229.exp[c('ABCA1','IKBIP','NPC2','TNFRSF1B','CTLA4','AKAP5'),rownames(GSE66229.cli.anno)],
         scale = 'row',name = 'Expr',
         color = colorRampPalette(c('navy', "white", 'yellow'))(100),
         cluster_rows = F,cluster_cols = F,
         show_rownames = T,show_colnames = F,
         annotation_col = GSE66229.cli.anno,
         annotation_colors = list(Risktype=risktype.col.use))
dev.off()




pdf('results/03.model/part1.pdf',height = 5.5,width = 20)
mg_merge_plot(tra.km,tra.roc,test.km,test.roc,ncol=4,labels = c('A','','B',''))
dev.off()



pdf('results/03.model/part3_3.pdf',height = 9,width = 5,onefile = F)
mg_merge_plot(GSE66229.km,GSE66229.roc,nrow=2)
dev.off()

#04.########
dir.create('results/04.nomogram')
tcga_cox_datas=tcga.risktype.cli
colnames(tcga_cox_datas)
table(tcga_cox_datas$Age1)


table(tcga_cox_datas$AJCC_stage)
tcga_cox_datas$AJCC_stage[tcga_cox_datas$AJCC_stage=='I'|tcga_cox_datas$AJCC_stage=='II']<-'I+II'
tcga_cox_datas$AJCC_stage[tcga_cox_datas$AJCC_stage=='III'|tcga_cox_datas$AJCC_stage=='IV']<-'III+IV'


table(tcga_cox_datas$Grade)
tcga_cox_datas$Grade[tcga_cox_datas$Grade=='G1'|tcga_cox_datas$Grade=='G2']<-'G1+G2'
tcga_cox_datas$Grade[tcga_cox_datas$Grade=='G3'|tcga_cox_datas$Grade=='G4']<-'G3+G4'


#######
#Age
tcga_cox_datas=crbind2DataFrame(tcga_cox_datas)
Age_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Age,
                             data=tcga_cox_datas))
Age_sig_cox_dat <- data.frame(Names=rownames(Age_sig_cox[[8]]),
                              HR = round(Age_sig_cox[[7]][,2],3),
                              lower.95 = round(Age_sig_cox[[8]][,3],3),
                              upper.95 = round(Age_sig_cox[[8]][,4],3),
                              p.value=round(Age_sig_cox[[7]][,5],3))
Age_sig_cox_dat

#Gender
Gender_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Gender,
                                data=tcga_cox_datas))
Gender_sig_cox_dat <- data.frame(Names=rownames(Gender_sig_cox[[8]]),
                                 HR = round(Gender_sig_cox[[7]][,2],3),
                                 lower.95 = round(Gender_sig_cox[[8]][,3],3),
                                 upper.95 = round(Gender_sig_cox[[8]][,4],3),
                                 p.value=round(Gender_sig_cox[[7]][,5],3))
Gender_sig_cox_dat
#AJCC_stage
Stage_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~AJCC_stage,
                               data=tcga_cox_datas))
Stage_sig_cox_dat <- data.frame(Names=rownames(Stage_sig_cox[[8]]),
                                HR = round(Stage_sig_cox[[7]][,2],3),
                                lower.95 = round(Stage_sig_cox[[8]][,3],3),
                                upper.95 = round(Stage_sig_cox[[8]][,4],3),
                                p.value=round(Stage_sig_cox[[7]][,5],3))
Stage_sig_cox_dat

#Grade
Grade_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Grade,
                               data=tcga_cox_datas))
Grade_sig_cox_dat <- data.frame(Names=rownames(Grade_sig_cox[[8]]),
                                HR = round(Grade_sig_cox[[7]][,2],3),
                                lower.95 = round(Grade_sig_cox[[8]][,3],3),
                                upper.95 = round(Grade_sig_cox[[8]][,4],3),
                                p.value=round(Grade_sig_cox[[7]][,5],3))
Grade_sig_cox_dat

#riskscore
riskscore_sig_cox<-summary(coxph(formula=Surv(OS.time, OS)~Riskscore,
                                 data=tcga_cox_datas))
riskscore_sig_cox_dat <- data.frame(Names=rownames(riskscore_sig_cox[[8]]),
                                    HR = round(riskscore_sig_cox[[7]][,2],3),
                                    lower.95 = round(riskscore_sig_cox[[8]][,3],3),
                                    upper.95 = round(riskscore_sig_cox[[8]][,4],3),
                                    p.value=round(riskscore_sig_cox[[7]][,5],3))
riskscore_sig_cox_dat

sig_cox_dat <- rbind(Age_sig_cox_dat,
                     Gender_sig_cox_dat,
                     Stage_sig_cox_dat,
                     Grade_sig_cox_dat,
                     riskscore_sig_cox_dat)
data.sig <- data.frame(Features=sig_cox_dat$Names,
                       p.value=sig_cox_dat$p.value,
                       sig_cox_dat$HR,
                       sig_cox_dat$lower.95,
                       sig_cox_dat$upper.95)
data.sig <- crbind2DataFrame(data.sig)
rownames(data.sig) <- c("Age","Gender","AJCC stage","Grade","RiskScore")
data.sig$Features=rownames(data.sig) 
data.sig
data.sig$p.value=ifelse(data.sig$p.value<0.001,'<0.001',data.sig$p.value)
pdf('results/04.nomogram/Univariate.pdf',height = 5,width = 6,onefile = F)
mg_forestplot_v2(data.sig,xlog = T,colgap =5,lineheight = 9,zero = 1,
                 boxsize = 0.3,lwd.zero=1,lwd.ci=1,lwd.xaxis=1
                 ,box_col='orange',summary_col="black",lines_col='black',zero_col='grey'
                 ,xlab='Hazard Ratio',lty.ci = 6,graph.pos =3)
dev.off()

write.csv(data.sig,'results/04.nomogram/Univariate analysis.csv')

#########
#T.stage +  N.stage +M.stage
muti_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Age+AJCC_stage+Riskscore, 
                              data=tcga_cox_datas))
muti_cox_dat <- data.frame(Names=rownames(muti_sig_cox[[8]]),
                           HR = round(muti_sig_cox[[7]][,2],3),
                           lower.95 = round(muti_sig_cox[[8]][,3],3),
                           upper.95 = round(muti_sig_cox[[8]][,4],3),
                           p.value=round(muti_sig_cox[[7]][,5],3))
data.muti <- data.frame(Features=muti_cox_dat$Names,
                        p.value=muti_cox_dat$p.value,
                        muti_cox_dat$HR,
                        muti_cox_dat$lower.95,
                        muti_cox_dat$upper.95)
data.muti <- crbind2DataFrame(data.muti)
data.muti
rownames(data.muti) <- c('Age',"AJCC stage","RiskScore")
data.muti$Features=rownames(data.muti)
data.muti
data.muti$p.value=ifelse(data.muti$p.value<0.001,'<0.001',data.muti$p.value)
pdf('results/04.nomogram/Multivariate.pdf',height = 5,width = 6,onefile = F)
mg_forestplot_v2(data.muti,xlog = T,colgap =5,lineheight = 12,zero = 1,
                 boxsize = 0.3,lwd.zero=1,lwd.ci=1,lwd.xaxis=1
                 ,box_col='orange',summary_col="black",lines_col='black',zero_col='grey'
                 ,xlab='Hazard Ratio',lty.ci = 6,graph.pos =3)
dev.off()
write.csv(data.muti,'results/04.nomogram/Multivariate analysis.csv')



pdf('results/04.nomogram/nomogram.pdf', width = 12, height = 10)
nom.plot=mg_nomogram(data.frame(RiskScore=tcga_cox_datas$Riskscore,
                                Age=tcga_cox_datas$Age,
                                AJCC_stage=tcga_cox_datas$AJCC_stage),
                     os = tcga_cox_datas$OS.time,
                     status = tcga_cox_datas$OS,
                     mks = c(1,3,5)
)
dev.off()
mg_nomogram_buti(nom.plot$Mod,cut.time = c(1*365,3*365,5*365))


#06.TME######
dir.create('results/05.TME_pathway')
df=read.csv
df[1:5,1:5]
tcga.immune=df[rownames(df)%in%tcga.cli$Samples,]
table(str_split_fixed(colnames(tcga.immune),'_',2)[,2])
p1=my_mutiboxplot( tcga.immune[tcga.risktype.cli$Samples,str_split_fixed(colnames(tcga.immune),'_',2)[,2]=='CIBERSORT'],
                group = tcga.risktype.cli$Risktype,legend.pos = 'top',group_cols = risktype.col,
                ylab = 'Fraction',fill = 'Risktype',size = 12)+ggtitle('CIBERSORT')

p2=my_mutiboxplot(tcga.immune[tcga.risktype.cli$Samples,str_split_fixed(colnames(tcga.immune),'_',2)[,2]=='TIMER'],
                group = tcga.risktype.cli$Risktype,legend.pos = 'none',group_cols = risktype.col,
                ylab = 'Score',fill = 'Risktype',size = 12,angle = 0,hjust = 0.5)+
  scale_x_discrete(labels=function(x) str_remove(x,"_TIMER"))+
  ggtitle('TIMER')+coord_flip()
  
p2

load('results/tcga.hall.ssGSEA.RData')
tcga.hall.ssGSEA[1:5,1:5]
rownames(tcga.hall.ssGSEA)=gsub('HALLMARK_','',rownames(tcga.hall.ssGSEA))
diff_pathway<-function(dat,group){
  dat=cbind.data.frame(cluster=group,t(dat))
  gr=names(table(group))
  dat1=dat[dat$cluster==gr[1],-1]
  dat2=dat[dat$cluster==gr[2],-1]
  # dat3=dat[dat$cluster==gr[3],-1]
  pathway=unique(colnames(dat)[-1])
  p_vale=data.frame(check.names = F)
  for (i in pathway){
    # x=c(dat1[,i],dat2[,i],dat3[,i])
    # g= factor(rep(names(table(group)), c(nrow(dat1), nrow(dat2), nrow(dat3))),
    #           labels = names(table(group)))
    # dd1=kruskal.test(x,g)$p.value
    dd1=wilcox.test(dat1[,i],dat2[,i])$p.value
    p_vale=rbind(p_vale,data.frame(pathway=i,p.value=dd1))
  }
  return(p_vale)
}
pathway.diff<-diff_pathway(dat=tcga.hall.ssGSEA[,tcga.risktype.cli$Samples],group=tcga.risktype.cli$Risktype)
table(pathway.diff$p.value<0.05)
pathway.diff.fit=pathway.diff[pathway.diff$p.value<0.05,]

library(ggcorrplot)
library(psych)
pathway_RS_cor <- corr.test(x =tcga.risktype.cli$Riskscore,
                            y = t(tcga.hall.ssGSEA[pathway.diff.fit$pathway,tcga.risktype.cli$Samples]),
                            method = "spearman",adjust = "BH",ci = F)


pathway_RS_cor_res=data.frame(pathway=pathway.diff.fit$pathway)
pathway_RS_cor_res$cor<-as.numeric(pathway_RS_cor$r)
pathway_RS_cor_res$p.adj<-as.numeric(pathway_RS_cor$p.adj)
head(pathway_RS_cor_res)
pathway_RS_cor_res=pathway_RS_cor_res[order(pathway_RS_cor_res$cor),]
head(pathway_RS_cor_res)
library(rcartocolor)
p3=ggplot(data=pathway_RS_cor_res,aes(x=cor,y=reorder(pathway,cor), color = -log10(p.adj))) +
  geom_point(aes(size=abs(cor)),show.legend = F) +
  scale_color_continuous(type = "gradient")+
  scale_color_gradient(low = "blue", high = "orange")+
  geom_segment(aes(yend=pathway,xend=0),size=.5) +
  labs(x='spearman Correlation',y='')+theme_classic()+
  theme(text = element_text(family = 'Times',size = 12))

pdf('results/05.TME_pathway/Fig8.pdf',height = 15,width = 14)
mg_merge_plot(p1,mg_merge_plot(p2,p3,labels = c('B','C'),widths = c(1,1.2)),
              nrow=2,labels = c('A',''),heights = c(1,1.5))
dev.off()

save.image(file = 'STAD_HP_scRNA.RData')
