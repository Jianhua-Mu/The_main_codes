setwd("C:/Users/28978/Desktop/pedi")
####GEO数据处理####
library(tidyverse)
#BiocManager::install('GEOquery')
library(GEOquery)
###下载数据，如果文件夹中有会直接读入
chooseBioCmirror()

gpl=getGEO('GPL16876', destdir=".")
gset1 = getGEO('GSE120559', destdir=".", AnnotGPL = F, getGPL = F)
gpl=getGEO('GPL16876', destdir=".")
gset7 = getGEO('GSE73517', destdir=".", AnnotGPL = F, getGPL = F)
gset4 = getGEO('GSE45547', destdir=".", AnnotGPL = F, getGPL = F)
# gset = getGEO('GSE85047', destdir=".", AnnotGPL = F, getGPL = F)
# gpl=getGEO('GPL5175', destdir=".")
colnames(Table(gpl))
head(Table(gpl))
GPL <- Table(gpl)[,c(1,18)]

grouplist=data.frame(factor(ifelse(substr(clin[,1],15,16)=="on",
                                   "Control","HSCR")))

clin1=pData(gset1[[1]])
clin4=pData(gset4[[1]])
clin7=pData(gset7[[1]])
grouplist1=data.frame(factor(ifelse(clin120559[,40]=="not amplified",
                  "MYCN-nonamplified","MYCN-amplified")))
grouplist4=data.frame(factor(ifelse(substr(clin45547[,37],53,54)==1,
                      "MYCN-nonamplified","MYCN-amplified")))
grouplist7=data.frame(factor(clin73517[,46]))
rownames(grouplist1)=colnames(exp120559)
rownames(grouplist4)=colnames(exp45547)
rownames(grouplist7)=colnames(exp73517)
colnames(grouplist1)=c("group")
colnames(grouplist7)=c("group")
colnames(grouplist4)=c("group")
save(grouplist1,file = "group_amp_1.Rda")
save(grouplist4,file = "group_amp_4.Rda")
save(grouplist7,file = "group_amp_7.Rda")
#没有symbol的处理方法
# GPL=Table(gpl)[,c(1,6)]
# GPL$GeneSymbol=str_split_fixed(GPL$gene_assignment , pattern = "//",3)[,2]
# GPL=GPL[,c(1,3)]
#
#GPL$GeneSymbol=str_split_fixed(GPL$DESCRIPTION , pattern = "'",9)[,8]

rownames(GPL)=GPL[,1]
#有时会报错  Increase it by setting `Sys.setenv("VROOM_CONNECTION_SIZE")`
Sys.setenv("VROOM_CONNECTION_SIZE"=99999999)
class(gset)
###提取子集
gset[[1]]
#读取表达谱
exp <- exprs(gset[[1]])
#把表达谱转为数据框格式
exp <- as.data.frame(exp)
##转换id
#读取GPL文件
comname <- intersect(rownames(exp),rownames(GPL))
exp <- exp[comname,]
GPL <- GPL[comname,]
exp1 <- cbind(GPL,exp)

exp1=exp1[exp1$GeneSymbol!="",]
exp1 <- exp1[!duplicated(exp1$GeneSymbol),]

rownames(exp1) <- exp1$GeneSymbol
exp1 <- exp1[,-(1:2)]


library(limma) 

boxplot(exp,outline=FALSE, notch=T, las=2,
        ylab="Expression",col = "#67C2A3")
exp=data.frame(normalizeBetweenArrays(exp))
boxplot(exp,outline=FALSE, notch=T, las=2,
        ylab="Expression",col = "#67C2A3")


library(RColorBrewer)
colors = brewer.pal(8,"Accent")
boxplot(exp,col = "white",cex=0.5, pch=20,outcol=colors[6], 
        main = "Example",ylab="Exp",names =colnames(exp))
boxplot(exp,outline=FALSE, notch=T, las=2)
range(exp)
exp <- log2(exp+1)
range(exp)
dev.off()
exp45547=data.frame(exp)

save(exp120559,file = "exp120559.Rda")
save(exp49710,file = "exp49710.Rda")
save(exp73517,file = "exp73517.Rda")
save(exp45547,file = "exp45547.Rda")

#arrayexpress
library(ArrayExpress)
# id="E-MTAB-2967"
# mexp = getAE(id, type = "full")

#去除批次效应
library(sva)
library(tidyverse)
library(limma)
load("exp120559.Rda")
load("exp49710.Rda")
load("exp73517.Rda")
load("exp45547.Rda")

#merge_eset=inner_join(exprSet2,exprSet_GSE46234,by="symbol")
merge_eset=cbind(exp120559[rownames(exp73517),],
                 exp45547[rownames(exp73517),],exp73517)
#rownames(merge_eset) <- merge_eset$symbol
#merge_eset <- merge_eset[,-1]
dim(merge_eset)
exp <- as.matrix(merge_eset)
dimnames <- list(rownames(exp),colnames(exp))
# data <- matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
# dim(data)
# class(data)
data=exp
batchType <- c(rep(1,208),rep(2,649),rep(3,105))
# modType <- c(rep("normal",6),rep("tumor",13),rep("normal",4),rep("tumor",4))
# mod  <-  model.matrix(~as.factor(modType))
# outTab <- data.frame(ComBat(data, batchType,mod, par.prior=TRUE))
#outTab <- data.frame(ComBat(data, batchType,par.prior=TRUE))
#出现负值说明有数据集有异常值(过大过小),可以比较两种方法哪个负值少
outTab2 <- data.frame(removeBatchEffect(data, batchType))
n=rep(T,ncol(outTab2))
for (i in 1:ncol(outTab2)){
  if (sum(outTab2[,i]<0)){
    n[i]=F
  }
}
exp=outTab2[,n]
save(exp,file = "exp_non_batch.Rda")
####差异分析####

library(limma)
grouplist=as.matrix(group932)
exp=exp
design=model.matrix(~grouplist)
fit=lmFit(exp,design)
fit=eBayes(fit)
deg=topTable(fit,coef=2,number = Inf)
##标记上下调基因
logFC=1
P.Value = 0.05
deg2=deg[(deg$P.Value < P.Value)&(abs(deg$logFC) > logFC),]
k1 = (deg$P.Value < P.Value)&(deg$logFC < -logFC)
k2 = (deg$P.Value < P.Value)&(deg$logFC > logFC)
deg$change = ifelse(k1,"down",ifelse(k2,"up","stable"))
table(deg$change)
gene_com=data.frame(gene_com)
colnames(gene_com)="gene"
##热图##
cg = rownames(deg)[deg$change !="stable"]
diff=exp[cg,]
library(pheatmap)
annotation_col=data.frame(group=grouplist)
rownames(annotation_col)=colnames(diff) 
pheatmap(diff,
         annotation_col=annotation_col,
         scale = "row",
         show_rownames = F,
         show_colnames =F,
         color = colorRampPalette(c("navy", "white", "red"))(50),
         fontsize = 10,
         fontsize_row=3,
         fontsize_col=3)

#RRA
library(RobustRankAggreg)
padj=0.05
logFC=1

files=c("GSE120559","GSE45547","GSE73517")
upList=list()
downList=list()
allFCList=list()
for(i in 1:3){
  inputFile=files[i]
  #rt=read.table(inputFile,header=T,sep = '\t',quote = '') 
  if (i==1){
    rt=deg11
  }
  if (i==2){
    rt=deg44
  }
  if (i==3){
    rt=deg77
  }
  #header=unlist(strsplit(inputFile,"_"))
  header=unlist(inputFile)
  downList[[header[1]]]=as.vector(rt[,1])
  upList[[header[1]]]=rev(as.vector(rt[,1]))
  fcCol=rt[,1:2]
  colnames(fcCol)=c("Gene",header[[1]])
  allFCList[[header[1]]]=fcCol
}


mergeLe=function(x,y){
  merge(x,y,by="Gene",all=T)}
newTab=Reduce(mergeLe,allFCList)
rownames(newTab)=newTab[,1]
newTab=newTab[,2:ncol(newTab)]
newTab[is.na(newTab)]=0


upMatrix = rankMatrix(upList)
upAR = aggregateRanks(rmat=upMatrix)
colnames(upAR)=c("Name","Pvalue")
upAdj=p.adjust(upAR$Pvalue,method="bonferroni")
upXls=cbind(upAR,adjPvalue=upAdj)
upFC=newTab[as.vector(upXls[,1]),]
upXls=cbind(upXls,logFC=rowMeans(upFC))

upSig=upXls[(upXls$adjPvalue<padj & upXls$logFC>logFC),]


downMatrix = rankMatrix(downList)
downAR = aggregateRanks(rmat=downMatrix)
colnames(downAR)=c("Name","Pvalue")
downAdj=p.adjust(downAR$Pvalue,method="bonferroni")
downXls=cbind(downAR,adjPvalue=downAdj)
downFC=newTab[as.vector(downXls[,1]),]
downXls=cbind(downXls,logFC=rowMeans(downFC))

downSig=downXls[(downXls$adjPvalue<padj & downXls$logFC< -logFC),]


allSig = rbind(upSig,downSig)
colnames(allSig)
allSig = allSig[,c("Name","logFC")]


hminput=newTab[c(as.vector(upSig[1:20,1]),as.vector(downSig[1:20,1])),]
library(pheatmap)
tiff(file="logFC.tiff",width = 15,height = 20,units ="cm",compression="lzw",bg="white",res=400)
pheatmap(hminput,display_numbers = TRUE,
         fontsize_row=10,
         fontsize_col=12,
         color = colorRampPalette(c("green", "white", "red"))(50),
         cluster_cols = FALSE,cluster_rows = FALSE, )
dev.off()

####WGCNA####
library("tidyverse")
library("WGCNA") 
datExpr0=
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK
###如果没有达标就需要筛选
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

# 样品聚类
# 聚类
sampleTree = hclust(dist(datExpr0), method = "average")
# 画图
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree)
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
###剪切线
abline(h = 67, col = "red")
###删除剪切线以下的样品
clust = cutreeStatic(sampleTree, cutHeight = 67, minSize = 10)
table(clust)
keepSamples = (clust==1)
datExpr0 = datExpr0[keepSamples, ]
dev.off()

# 记录基因和样本数，方便后续可视化
nGenes = ncol(datExpr0)#基因数
nSamples = nrow(datExpr0)#样本数

# 构建网络，识别模块
# power值散点图
enableWGCNAThreads()   #多线程工作
powers = c(1:20)       #幂指数范围1:20
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)

par(mfrow = c(1,2))
cex1 = 0.9
###拟合指数与power值散点图
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red") #可以修改
###平均连通性与power值散点图
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

##6.2 邻接矩阵转换
sft #查看最佳power值
softPower =sft$powerEstimate #最佳power值
softPower = as.numeric(softPower)
adjacency = adjacency(datExpr0, power = softPower)

##6.3 TOM矩阵
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM
save(TOM,file = "TOM_CTF1.Rda")

# 基因聚类
geneTree = hclust(as.dist(dissTOM), method = "average");
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)

# 动态剪切模块识别
minModuleSize = 30      #模块基因数目
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)

dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

# 相似模块聚类
MEList = moduleEigengenes(datExpr0, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average")
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

MEDissThres = 0.1 #剪切高度可修改
abline(h=MEDissThres, col = "red")

###相似模块合并
merge = mergeCloseModules(datExpr0, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs
plotDendroAndColors(geneTree, mergedColors,"Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

moduleColors = mergedColors
table(moduleColors)
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs
dev.off()
# 整理临床信息
clinical <- read.table("GBM_fpkm_mRNA_01A_estimate_score.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
clinical <- clinical[rownames(datExpr0),]
identical(rownames(clinical),rownames(datExpr0))
# 查看临床信息
head(clinical)
# 对表达矩阵进行预处理
datTraits = as.data.frame(do.call(cbind,lapply(clinical, as.numeric)))
rownames(datTraits) = rownames(clinical)

# 对样本进行聚类
sampleTree2 = hclust(dist(datExpr0), method = "average")

# 将临床信息转换为颜色，白色表示低，红色表示高，灰色表示缺失
traitColors = numbers2colors(datTraits, signed = FALSE)

# 样本聚类图与样本性状热图
plotDendroAndColors(sampleTree2, 
                    traitColors,
                    groupLabels = names(datTraits), 
                    main = "Sample dendrogram and trait heatmap")
dev.off()

#### 网络的分析
# 对模块特征矩阵进行排序
MEs=orderMEs(MEs)
#计算模型特征矩阵和样本信息矩阵的相关度。
moduleTraitCor=cor(MEs, datTraits, use="p")

moduleTraitPvalue=corPvalueStudent(moduleTraitCor, nSamples)


#使用labeledHeatmap()将上述相关矩阵和p值可视化。
textMatrix=paste(signif(moduleTraitCor,2),"\n(",signif(moduleTraitPvalue,1),")",sep="")
dim(textMatrix)=dim(moduleTraitCor)
# 基因模块与临床信息相关性图
labeledHeatmap(Matrix=moduleTraitCor,#模块和表型的相关性矩阵，这个参数最重要，其他可以不变
               xLabels=colnames(datTraits),
               yLabels=names(MEs),
               ySymbols=names(MEs),
               colorLabels=FALSE,
               colors=blueWhiteRed(50),
               textMatrix=textMatrix,
               setStdMargins=FALSE,
               cex.text=0.7,
               cex.lab=0.7,
               zlim=c(-1,1),
               main=paste("Module-trait relationships"))
dev.off()

# 不同模块与基因性状的具体分析
##矩阵一
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"))
####看一下目的基因和哪个模块相关性最高
a <- geneModuleMembership
a <- a %>% rownames_to_column()

MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

##矩阵二
traitNames=names(datTraits)
geneTraitSignificance = as.data.frame(cor(datExpr0, datTraits, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", traitNames, sep="")
names(GSPvalue) = paste("p.GS.", traitNames, sep="")

##批量输出性状和模块散点图
for (trait in traitNames){
  traitColumn=match(trait,traitNames)  
  for (module in modNames){
    column = match(module, modNames)
    moduleGenes = moduleColors==module
    if (nrow(geneModuleMembership[moduleGenes,]) > 1){
      outPdf=paste(trait, "_", module,".pdf",sep="")
      pdf(file=outPdf,width=7,height=7)
      par(mfrow = c(1,1))
      verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                         abs(geneTraitSignificance[moduleGenes, traitColumn]),
                         xlab = paste("Module Membership in", module, "module"),
                         ylab = paste("Gene significance for ",trait),
                         main = paste("Module membership vs. gene significance\n"),
                         cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
      abline(v=0.8,h=0.5,col="red")
      dev.off()
    }
  }
}

#10. 输出每个模块的基因
for (mod in 1:nrow(table(moduleColors)))
{  
  modules = names(table(moduleColors))[mod]
  probes = colnames(datExpr0)
  inModule = (moduleColors == modules)
  modGenes = probes[inModule]
  write.table(modGenes, file =paste0(modules,".txt"),sep="\t",row.names=F,col.names=F,quote=F)
}

####svm-rfe####
# 加载库和包
set.seed(100)
library(mlbench)
library(caret)

# 加载数据集
data(PimaIndiansDiabetes)

control <- rfeControl(functions = caretFuncs, method = "cv", number = 30)
#number没啥用
data=cbind(mtab_data[,c(1:8)],mtab_data[,gene_com])
#data$OS=factor(data$OS,levels = c(0,1))
# 执行SVM-RFE算法
results <- rfe(data[,9:ncol(data)], #参数
               data[,4], #目标
               sizes = c(1:50), 
               rfeControl = control,
               method = "svmRadial")

# 结果分析
print(results)
# 列出选择的变量集
svm_gene=predictors(results)
# 绘制结果
plot(results, type=c("g", "o"),main="SVM-RFE")
####randromforest####
set.seed(100)
library(randomForest)

train=data[,c(4,9:ncol(data))]
train$OS=factor(train$OS,levels = c(0,1))

err<-as.numeric()
for(i in 9:(length(names(train)))){
  mtry_test <- randomForest(OS~., data=train, mtry=i)
  err<- append( err, mean( mtry_test$err.rate ) )}
print(err)
mtry<-which.min(err)
ntree_fit<-randomForest(OS~., data=train, mtry=mtry, ntree=1000)
plot(ntree_fit,main="RandomForest")

rf<-randomForest(OS~., data=train, mtry=mtry, ntree=900, importance=T )
rf
imp=data.frame(importance(rf))
mda=rownames(imp[order(imp$MeanDecreaseAccuracy,decreasing = T),])[1:50]
mdg=rownames(imp[order(imp$MeanDecreaseGini,decreasing = T),])[1:50]
rf_gene=intersect(mda,mdg)
varImpPlot(rf)

p=predict(ntree_fit,newdata = cgga218[201:218,gene])
sum(p == rt[201:218,2])/18

#Cox分析####
library(survival)
library(forestplot)
library(tidyverse)

surv.expr=mtab_data[,c("OS","OS_time",gene)]
Coxoutput <- NULL 
for(i in 3:ncol(surv.expr)){
  g <- colnames(surv.expr)[i]
  cox <- coxph(Surv(OS_time,OS) ~ surv.expr[,i], data = surv.expr) # 单变量cox模型
  coxSummary = summary(cox)
  
  Coxoutput <- rbind.data.frame(Coxoutput,
                                data.frame(gene = g,
                                           HR = as.numeric(coxSummary$coefficients[,"exp(coef)"])[1],
                                           z = as.numeric(coxSummary$coefficients[,"z"])[1],
                                           pvalue = as.numeric(coxSummary$coefficients[,"Pr(>|z|)"])[1],
                                           lower = as.numeric(coxSummary$conf.int[,3][1]),
                                           upper = as.numeric(coxSummary$conf.int[,4][1]),
                                           stringsAsFactors = F),
                                stringsAsFactors = F)
}



###筛选top基因
pcutoff <- 0.05

topgene <- Coxoutput[Coxoutput$pvalue<pcutoff,]
  # 取出p值小于阈值的基因
topgene=topgene[order(topgene$pvalue),]
# write.table(topgene[1:100,], file = "topgene.txt",sep = "\t",
#             row.names = F,col.names = T,quote = F)
# topgene1=topgene[topgene$HR<3,]
num=as.matrix(15)
topgene=topgene[order(topgene$pvalue),][1:num,]
topgene=topgene[order(topgene$HR),]
#3. 绘制森林图
##3.1 输入表格的制作
tabletext <- cbind(c("Gene",topgene$gene),
                   c("HR",format(round(as.numeric(topgene$HR),3),nsmall = 3)),
                   c("lower 95%CI",format(round(as.numeric(topgene$lower),3),nsmall = 3)),
                   c("upper 95%CI",format(round(as.numeric(topgene$upper),3),nsmall = 3)),
                   c("pvalue",format(round(as.numeric(topgene$p),3),nsmall = 3)))
##3.2 绘制森林图

n= as.character(num+2)
forestplot(labeltext=tabletext,
           mean=c(NA,as.numeric(topgene$HR)),
           lower=c(NA,as.numeric(topgene$lower)), 
           upper=c(NA,as.numeric(topgene$upper)),
           graph.pos=5,# 图在表中的列位置
           graphwidth = unit(.25,"npc"),# 图在表中的宽度比
           fn.ci_norm="fpDrawDiamondCI",# box类型选择钻石
           col=fpColors(box="#00A896", lines="#02C39A", zero = "black"),# box颜色
           
           boxsize=0.4,# box大小固定
           lwd.ci=1,
           ci.vertices.height = 0.1,ci.vertices=T,# 显示区间
           zero=1,# zero线横坐标
           lwd.zero=1.5,# zero线宽
           xticks = c(0.5,1,1.5),# 横坐标刻度根据需要可随意设置
           lwd.xaxis=2,
           xlab="Hazard ratios",
           txt_gp=fpTxtGp(label=gpar(cex=1.2),# 各种字体大小设置
                          ticks=gpar(cex=0.85),
                          xlab=gpar(cex=1),
                          title=gpar(cex=1.5)),
           hrzl_lines=list("1" = gpar(lwd=2, col="black"), # 在第一行上面画黑色实线
                           "2" = gpar(lwd=1.5, col="black"), # 在第一行标题行下画黑色实线
                           "22" = gpar(lwd=2, col="black")), # 在最后一行上画黑色实线
           lineheight = unit(.75,"cm"),# 固定行高
           colgap = unit(0.3,"cm"),
           mar=unit(rep(1.5, times = 4), "cm"),
           new_page = F
)

#韦恩图####

#新版
library(ggVennDiagram)
library(ggplot2)
mydata<-list(Agroup=-12:294,Bgroup=197:1302)
ggVennDiagram(mydata)


#旧版
library(VennDiagram)##载入需要的程辑包:grid
##载入需要的程辑包: futile.logger
set1 <- sample( 1 : 1000,300, replace= F)
set2 <- sample( 1 : 1000,130,replace = F)
set3 <- sample( 1 : 1000,300, replace = F)
set4 <- sample( 1 : 1000,300,replace = F)
Random_forest=rf_gene
Svm_rfe=svm_gene
Cox_regression=cox_gene




#3个数据集
s2 <- list(
  set1 = Random_forest,set2 = Svm_rfe,set3 = Cox_regression)
v2 <- venn.diagram(x = s2,filename = NULL,
                   col="white",
                   fill=c(colors()[616], colors()[130],colors()[462]),
                   alpha=c(0.6,0.6,0.6),
                   lwd=c(1,1,1),
                   cex=2,
                   #cat.dist=c(0.05,8.05,-0.45),
                   cat.pos=c(0,0,0),
                   cat.cex=2,
                   scaled = F
)
cowplot::plot_grid(v2)


##4个数据集
s3 <- list(
  set1 = set1,set2 = set2,set3 = set3,set4 = set4)
v3 <- venn.diagram(x = s3, filename = NULL,
                     #height = 450,
                     #width = 450,#resolution =300,#imagetype="png",col="transparent",
                     fill=c("cornflowerblue" , "green" , "yellow" , "darkorchid1"),
                   alpha = 0.50,
                     cex=2,
                     cat.cex=2)
cowplot::plot_grid(v3)
  

#LASSO####
library(survminer)
library(survival)
library(ggpubr)
library(ROCR)
library(glmnet)
library(caret)
library(timeROC)
load("mtab表达临床.Rda")
load("C:/Users/28978/Desktop/pedi/gene_com.Rda")
rt1=data.frame(mtab_data[,c(4,6)])
rt1$OS=as.numeric(rt1$OS)

gene=c("GNAI1","CNR1","PLXNC1","NTS","GLDC","ABCC4" )#未加HD
gene=c("CNR1","GNAI1","GLDC","ABCC4")
rt=mtab_data

rt=rt[,gene]
rt=cbind(rt1,rt)

x=data.matrix(rt[,c(3:ncol(rt))])
#这里要求Surv第一个是时间，第二个参数是状态
y=data.matrix(Surv(rt1[,2],rt1[,1]))
fit=glmnet(x, y, family = "cox", maxit = 10000)
#plot(fit, xvar = "lambda", label = TRUE)

cvfit = cv.glmnet(x, y, family="cox", maxit = 10000)
#plot(cvfit)
#其中两条虚线分别指示了两个特殊的λ值

###4. 输出预测模型的相关系数与riskScore
###4.1 输出相关系数
coef=coef(fit, s = cvfit$lambda.min)
index=which(coef != 0)
actCoef=coef[index]
lassoGene=row.names(coef)[index]
geneCoef=cbind(Gene=lassoGene,Coef=actCoef)
geneCoef   #查看模型的相关系数

FinalGeneExp = rt[,lassoGene]
myFun = function(x){crossprod(as.numeric(x),actCoef)}
riskScore = apply(FinalGeneExp,1,myFun)
outCol = c("OS_time", "OS", lassoGene)
risk = as.vector(ifelse(riskScore > mean(riskScore), "high", "low"))
dat = cbind(rt[,outCol], riskScore=as.vector(riskScore), risk)

p <- ggboxplot(dat, x = "OS", y = "riskScore",
               color = "OS", palette = "jco",
               add = "jitter")
p <- p + stat_compare_means()   #  Add p-value
p   #得出预测结果

pred <- prediction(dat$riskScore, dat$OS)
perf <- performance(pred,"tpr","fpr")
AUC <- performance(pred,"auc")   #计算AUC
auc=AUC@y.values[[1]]
plot(perf,colorize=FALSE, col="red", print.auc =TRUE,main=paste("AUC=",auc)) #绘制ROC曲线
lines(c(0,1),c(0,1),col = "gray", lty = 4 )

surv <- rt
surv$OS_time=surv$OS_time/365
surv$group <- ifelse(risk =="high","High","Low")
surv$group <- factor(surv$group, levels = c("Low","High")) 
class(surv$group)
table(surv$group)

fitd <- survdiff(Surv(OS_time, OS) ~ group,
                 data      = surv,
                 na.action = na.exclude)
pValue <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)

#2.2 拟合生存曲线

fit <- survfit(Surv(OS_time, OS)~ group, data = surv)
summary(fit)

p.lab <- paste0("P", ifelse(pValue < 0.001, " < 0.001", paste0(" = ",round(pValue, 3))))
ggsurvplot(fit,
           data = surv,
           pval = p.lab,
           conf.int = TRUE, # 显示置信区间
           risk.table = TRUE, # 显示风险表
           risk.table.col = "strata",
           palette = "jco", # 配色采用jco
           legend.labs = c("Low", "High"), # 图例
           size = 1,
           xlim = c(0,15), # x轴长度，一般为0-10年
           break.time.by = 2, # x轴步长为20个月
           legend.title = "",
           surv.median.line = "hv", # 限制垂直和水平的中位生存
           ylab = "Survival probability (%)", # 修改y轴标签
           xlab = "Time (Months)", # 修改x轴标签
           ncensor.plot = TRUE, # 显示删失图块
           ncensor.plot.height = 0.25,
           risk.table.y.text = FALSE)




surv$riskScore=riskScore
ROC3 <- timeROC(T=surv$OS_time,   #结局时间
                delta=surv$OS,   #结局指标
                marker=surv$riskScore,   #预测变量
                cause=1,   #阳性结局指标数值
                weighting="marginal",   #计算方法，默认为marginal
                times=c(1,2,3,4,5),   #时间点，选取1年，3年和5年的生存率
                iid=TRUE)
ROC3

load("target表达临床.Rda")
rt_test=target_data
rt_test$OS=as.numeric(rt_test$OS)
test_exp=rt_test
FinalGeneExp = test_exp[,lassoGene]
myFun = function(x){crossprod(as.numeric(x),actCoef)}
riskScore = apply(FinalGeneExp,1,myFun)
outCol = c("OS_time", "OS", lassoGene)
risk = as.vector(ifelse(riskScore > mean(riskScore), "high", "low"))
dat2 = cbind(test_exp[,outCol], riskScore=as.vector(riskScore), risk)
surv=dat2

p <- ggboxplot(dat2, x = "OS", y = "riskScore",
               color = "OS", palette = "jco",
               add = "jitter")
p <- p + stat_compare_means()   #  Add p-value
p   #得出预测结果


surv$OS_time <- surv$OS_time/365
#median
surv$group <- ifelse(risk =="high","High","Low")
surv$group <- factor(surv$group, levels = c("Low","High")) 
table(surv$group)

fitd <- survdiff(Surv(OS_time, OS) ~ group,
                 data      = surv,
                 na.action = na.exclude)
pValue <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)

#2.2 拟合生存曲线

fit <- survfit(Surv(OS_time, OS)~ group, data = surv)
summary(fit)

p.lab <- paste0("P", ifelse(pValue < 0.001, " < 0.001", paste0(" = ",round(pValue, 3))))
ggsurvplot(fit,
           data = surv,
           pval = p.lab,
           conf.int = TRUE, # 显示置信区间
           risk.table = TRUE, # 显示风险表
           risk.table.col = "strata",
           palette = "jco", # 配色采用jco
           legend.labs = c("Low", "High"), # 图例
           size = 1,
           xlim = c(0,15), # x轴长度，一般为0-10年
           break.time.by = 2, # x轴步长为20个月
           legend.title = "",
           surv.median.line = "hv", # 限制垂直和水平的中位生存
           ylab = "Survival probability (%)", # 修改y轴标签
           xlab = "Time (Months)", # 修改x轴标签
           ncensor.plot = TRUE, # 显示删失图块
           ncensor.plot.height = 0.25,
           risk.table.y.text = FALSE)

surv$riskScore=riskScore
ROC3 <- timeROC(T=surv$OS_time,   #结局时间
                delta=surv$OS,   #结局指标
                marker=surv$riskScore,   #预测变量
                cause=1,   #阳性结局指标数值
                weighting="marginal",   #计算方法，默认为marginal
                times=c(1,2,3,4,5),   #时间点，选取1年，3年和5年的生存率0
                iid=TRUE)

ROC3


#列线图####
#Nomogram
load("dat_lasso325训练集dat数据.Rda")
clinical<- mtab_data[3:223,1:8]
  # clinical$submitter_id.samples <- gsub("-",".",clinical$submitter_id.samples)
# rownames(clinical) <- clinical$submitter_id.samples
# clinical <- clinical[,-1]
##结合dat与clinical
comgene1 <- intersect(rownames(clinical),rownames(dat))
dat <- dat[comgene1,]
clinical <- clinical[comgene1,]
dat1 <- dat[,(5:(ncol(dat)-1))]
#clinical <- clinical[,c("Gender","Age","IDH_mutation_status","MGMTp_methylation_status")]
#clinical <- clinical1
identical(rownames(dat),rownames(clinical))
rt2 <- cbind(dat,clinical)

rt2=rt2[complete.cases(rt2),]
rt2$Age=rt2$Age/365
rt2$Age=ifelse(rt2$Age>1,">1","<=1")
# write.csv(rt2, file = "rt2.csv")
#导入修改后的rt3

#加载包
library(rms)
library(foreign)
library(survival)
###3. 设置参数

rt2$Gender <- factor(rt2$Gender,labels=c("Famale", "Male"))
rt2$Age <- factor(rt2$Age,labels=c(">1","<=1"))
rt2$Tert_status <- factor(rt2$Tert_status,labels=c("wild type", "TERT-RA"))
rt2$Mycn_status <- factor(rt2$Mycn_status,
                                       labels=c("non amplified", "amplified"))

rt2$risk <- factor(rt2$risk,labels=c("low", "high"))
rt2=rt2[,5:ncol(rt2)]
rt2$OS_time=rt2$OS_time/365
ddist <- datadist(rt2)
options(datadist='ddist')   #使用函数datadist()将数据打包

###4. 构建Cox回归模型
f <- cph(Surv(OS_time, OS) ~Age +Tert_status + Mycn_status+ risk, x=T,
         y=T, surv=T, data=rt2, time.inc=1)
surv <- Survival(f)

###5. 构建Nomogram
nom2 <- nomogram(f, fun=list(function(x) surv(1, x), 
                             function(x) surv(2, x), function(x) surv(3, x)), 
                 lp=F, funlabel=c("1-year survival", "2-year survival", "3-year survival"), 
                 maxscale=100, 
                 fun.at=c(0.99, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3,0.2,0.1,0.05))
plot(nom2)

# library(regplot)#另一种画法
# regplot(f,
#         interval="confidence",
#         title="Nomogram",
#         clickable = T)

#下面绘制校准曲线
par(mfrow = c(1,3))

#不同时间节点只需要改变time.inc参数就行
f5 <- cph(Surv(OS_time, OS) ~ Age +Tert_status + Mycn_status+ risk, 
          x=T, y=T, surv=T, data=rt2, time.inc=5)
cal5 <- calibrate(f5, cmethod="KM", method="boot", u=5,m=43,B=1000)
par(mar=c(6,5,1,2),cex = 1.0)#图形参数
plot(cal5,lwd=2,lty=1,subtitles = F,#关闭副标题
     errbar.col=c(rgb(0,0,0,maxColorValue=255)),
     xlim=c(0,1),ylim=c(0,1),
     xlab="Predicted 5-year Overall Survival",ylab="Actual 5-year Overall Survival",
     col=c(rgb(255,0,0,maxColorValue=255))
)
abline(0,1,lty=3,lwd=2,col="blue")



#免疫浸润####
library(e1071)
library(parallel)
library(preprocessCore)
source("CIBERSORT.R")   
sig_matrix <- "LM22.txt"   
mixture_file = 'exp_mtab.txt'   #肿瘤患者表达谱,除去列名，第一列也要是基因名称
res_cibersort <- CIBERSORT(sig_matrix, mixture_file, perm=100, QN=TRUE)
save(res_cibersort,file = "res_cibersort_mtab.Rdata")   #保存中间文件
write.table(cbind(rownames(met_exp[comgene,]),met_exp[comgene,]),file="exp_mtab.txt",sep = '\t',row.names = T,col.names = T,quote = F)
# write.table(res_cibersort,file="CIBERSORT-Results.txt",sep = '/t',row.names = T,col.names = NA,quote = F)
res_cibersort <- res_cibersort[,1:22]   
ciber.res <- res_cibersort[,colSums(res_cibersort) > 0]   #去除丰度全为0的细胞
#可视化）
mycol <- ggplot2::alpha(rainbow(ncol(ciber.res)), 0.9) #创建彩虹色板（带70%透明度）
par(bty="o", mgp = c(2.5,0.3,0), mar = c(2.1,4.1,2.1,10.1),tcl=-.25,las = 1,xpd = F)
barplot(as.matrix(t(ciber.res)),
        border = NA, # 柱子无边框写
        names.arg = rep("",nrow(ciber.res)), # 无横坐标样本名
        yaxt = "n", # 先不绘制y轴
        ylab = "Relative percentage", # 修改y轴名称
        col = mycol) # 采用彩虹色板
axis(side = 2, at = c(0,0.2,0.4,0.6,0.8,1), # 补齐y轴添加百分号
     labels = c("0%","20%","40%","60%","80%","100%"))
legend(par("usr")[2], # 这里-20要根据实际出图的图例位置情况调整
       par("usr")[4], 
       legend = colnames(ciber.res), 
       xpd = T,
       fill = mycol,
       cex = 0.6, 
       border = NA, 
       y.intersp = 1,
       x.intersp = 0.2,
       bty = "n")
dev.off()   #关闭画板

#分组浸润程度箱线图
library(tidyverse)
a=res_cibersort
a <- read.table("CIBERSORT-Results.txt", sep = "\t",row.names = 1,check.names = F,header = T)
a=data.frame(a)
a <- a[,1:22]

load("group.Rda")
group=data.frame(as.factor(ifelse(risk=="high","high","low")))
colnames(group)="group"
rownames(group)=gsub("-",".",rownames(target_data))
identical(rownames(a),rownames(group))
#a=a[intersect(rownames(a),group$sample),]

b <- group
class(b$group)
a$group <- b$group
a <- a %>% rownames_to_column("sample")
library(ggsci)
library(tidyr)
library(ggpubr)

b <- gather(a,key=CIBERSORT,value = Proportion,-c(group,sample))

ggboxplot(b, x = "CIBERSORT", y = "Proportion",
          fill = "group", palette = "lancet")+
  stat_compare_means(aes(group = group),
                     method = "wilcox.test",
                     label = "p.signif",
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns")))+
  theme(text = element_text(size=10),
        axis.text.x = element_text(angle=45, hjust=1)) 

library(corrplot)
library(tidyverse)

gene=c("CNR1","GNAI1","GLDC","ABCC4")
exp <- mtab_data[,gene]
exp <- exp %>% as.data.frame()
load("res_cibersort_mtab.Rdata")
ciber=res_cibersort[,1:22]
ciber=data.frame(ciber)
ciber <- ciber[,1:22]
identical(rownames(ciber),rownames(exp))
class(exp$CNR1)
class(ciber$B.cells.naive)
cor<-sapply(ciber,function(x,y) cor(x,y,method="spearman"),exp)
rownames(cor)<-colnames(exp)
cor_res <- cor.mtest(cor,#计算p值
                     conf.level = 0.95)#置信区间
p=cor_res[["p"]]

corrplot(cor,
         method = "color",#相关性矩阵展示的图形
         col=colorRampPalette(c("#01468b","white","#ee0000"))(100),
         addCoef.col = "black",#为相关系数添加颜色
         tl.col="black",#设置文本标签的颜色
         number.cex = 0.5,
         tl.cex = 0.7,
         cl.align = "l",type = "lower",add = T
)

##相关性网络图####
library(igraph)
library(reshape2)
   
cutoff=0.4  #相关性阈值 
#读取输入文件
data=cbind(mtab_data[,gene],results)
cordata=cor(data)

#保留相关性矩阵的一半
mydata = cordata
upper = upper.tri(mydata)
mydata[upper] = NA

#把相关性矩阵转换为数据框
df = data.frame(gene=rownames(mydata),mydata)
dfmeltdata = melt(df,id="gene")
dfmeltdata = dfmeltdata[!is.na(dfmeltdata$value),]
dfmeltdata = dfmeltdata[dfmeltdata$gene!=dfmeltdata$variable,]
dfmeltdata = dfmeltdata[abs(dfmeltdata$value)>cutoff,]

#定义网络图的节点和边
corweight = dfmeltdata$value
weight = corweight+abs(min(corweight))+5
d = data.frame(p1=dfmeltdata$gene,p2=dfmeltdata$variable,weight=dfmeltdata$value)
g = graph.data.frame(dfmeltdata,directed = FALSE)
#设置颜色，节点大小，字体大小
E(g)$color = ifelse(corweight>0,rgb(254/255,67/255,101/255,
                                    abs(corweight)),
                    rgb(0/255,0/255,255/255,abs(corweight)))
V(g)$size = 8
V(g)$shape = "circle"
V(g)$lable.cex = 1.2
V(g)$color = "white"
E(g)$weight = weight

#可视化

layout(matrix(c(1,1,1,0,2,0),byrow=T,nc=3),height=c(6,1),width=c(3,4,3))
par(mar=c(1.5,2,2,2))
vertex.frame.color = NA
plot(g,layout=layout_nicely,
     vertex.label.cex=V(g)$lable.cex,
     edge.width = E(g)$weight,
     edge.arrow.size=0,
     vertex.label.color="black",
     vertex.frame.color=vertex.frame.color,
     edge.color=E(g)$color,
     vertex.label.cex=V(g)$lable.cex,
     vertex.label.font=2,
     vertex.size=V(g)$size,
     edge.curved=0.4)

#绘制图例
color_legend = c(rgb(254/255,67/255,101/255,seq(1,0,by=-0.01)),rgb(0/255,0/255,255/255,seq(0,1,by=0.01)))
par(mar=c(2,2,1,2),xpd = T,cex.axis=1.6,las=1)
barplot(rep(1,length(color_legend)),border = NA, space = 0,ylab="",xlab="",
        xlim=c(1,length(color_legend)),horiz=FALSE,
        axes = F, col=color_legend,main="")
axis(3,at=seq(1,length(color_legend),length=5),c(1,0.5,0,-0.5,-1),tick=FALSE)
dev.off()


#干性指数####
library(synapser)
synLogin(email = "mjh289@126.com", password = "aaa15281960846")
synRNA <- synGet( "syn2701943", downloadLocation = "~/Downloads/PCBC" )
library(tidyverse)

exp <- read.delim("C:/Users/28978/Desktop/pedi/rnaseq_norm.tsv") %>%
  # 去除 Ensembl ID 的后缀
  separate(col = "tracking_id", sep = "\\.", into = c("Ensembl_ID", "suffix")) %>%
  dplyr::select(-suffix) %>%
  column_to_rownames("Ensembl_ID") %>%
  as.matrix()

library(org.Hs.eg.db)

unimap <- mapIds(
  org.Hs.eg.db, keys = rownames(exp), keytype = "ENSEMBL", 
  column = "SYMBOL", multiVals = "filter"
)

data.exp <- exp[names(unimap),]
rownames(data.exp) <- unimap
synMeta <- synTableQuery("SELECT UID, Diffname_short FROM syn3156503")
metaInfo <- synMeta$asDataFrame() %>%
  dplyr::select(UID, Diffname_short) %>%
  column_to_rownames("UID") %>%
  filter(!is.na(Diffname_short))

X <- data.exp  
y <- metaInfo[colnames(X), ]
names(y) <- colnames(X)
library(gelnet)
gelnet(X, y, l1, l2)
# 对数据进行均值中心化
X <- data.exp
m <- apply(X, 1, mean)
X <- X - m
# 将样本分为干细胞组和非干细胞组
sc <- which(y == "SC")
X.sc <- X[, sc]
X.or <- X[, -sc]

model.RNA <- gelnet(t(X.sc), NULL, 0, 1)

save(X, y, model.RNA, file = "mRNAsi_model.rda")

load('mRNAsi_model.rda')
exp=data.frame(t(apply(target_data, 2, as.numeric)))
common <- intersect(names(model.RNA$w), rownames(exp))
X <- data.frame(exp[common, ])
w <- model.RNA$w[common]
score <- apply(X, 2, function(z) {cor(z, w, method="sp", use="complete.obs")})
score <- score - min(score)
score <- score / max(score)

predict.mRNAsi <- function(exp, modelPath='model.rda') {
  load(modelPath)
  
  common <- intersect(names(model.RNA$w), rownames(exp))
  X <- exp[common, ]
  w <- model.RNA$w[common]
  
  score <- apply(X, 2, function(z) {cor(z, w, method="sp", use="complete.obs")})
  score <- score - min(score)
  score <- score / max(score)
}
score=data.frame(score)
rownames(score)=rownames(mtab_data)
colnames(score)="mRNAsi"
save(score_tar,file = "mRNAsi_tar.Rda")

#免疫评分####
library(utils) #这个包应该不用下载，自带的
rforge <- "http://r-forge.r-project.org"
#install.packages("estimate", repos=rforge, dependencies=TRUE)
library(estimate)
library(tidyverse)
#读取肿瘤患者01A表达谱
expr <- read.table("exp_target.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)


#计算免疫评分
filterCommonGenes(input.f = "exp_target.txt",   #输入文件名
                  output.f = "exp_target.gct",   #输出文件名
                  id = "GeneSymbol")   #行名为gene symbol
estimateScore("exp_target.gct",   #刚才的输出文件名
              "exp_target_estimate_score.txt",   #新的输出文件名（即估计的结果文件）
              platform="affymetrix")   #默认平台

#3. 输出每个样品的打分
result <- read.table("exp_target_estimate_score.txt",
                     sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
result <- result[,-1]   
colnames(result) <- result[1,]   
result <- as.data.frame(t(result[-1,]))

rownames(result) <- colnames(expr)
save(result,file = "评分137.Rda")
# write.table(result, file = "exp_target_estimate_score.txt",
#             sep = "\t",row.names = T,col.names = NA,quote = F) # 保存并覆盖得分
e=data.frame(target_data[,gene])
colnames(e)=gene
e=apply(result, 2, as.numeric)
result=apply(result, 2, as.numeric)
R=cor(cbind(result,score_tar,e))


#雷达图####
library(fmsb)
data=data.frame(t(R[1:5,6:9]))
data=rbind(rep(1,ncol(data)),rep(-1,ncol(data)),data)
# data <- as.data.frame(matrix( sample( 2:20 , 10 , replace=T) , ncol=10))
# colnames(data) <- c("A" , "B" , "C" , "D" , "E", "F" , "G" ,
#                     "H", "I", "J" )
# data <- rbind(rep(20,10) , rep(0,10) , data)
# radarchart(data, axistype=1 ,
#            pcol=rgb(0.2,0.5,0.5,0.9) , pfcol=rgb(0.2,0.5,0.5,0.5) , plwd=4 ,
#            cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,20,5),
#            cglwd=0.8,
#            vlcex=0.8
# )

#多个数据，第1，2行分别是上限，下限
radarchart(data)
#定义radar图函数
create_beautiful_radarchart <-function(data, color ="#00AFBB", 
                                       vlabels = colnames(data), vlcex =0.7, 
                                       caxislabels = NULL, title = NULL, ...){
  radarchart(data, axistype = 1,
             pcol = color, pfcol = scales::alpha(color, 0.1), plwd = 2, plty = 1,
             cglcol = "grey", cglty = 1, cglwd = 0.8,
             axislabcol = "grey",
             vlcex = vlcex, vlabels = vlabels,
             caxislabels = caxislabels, title = title, ...
  )
}
#画美化后图形
op <- par(mar = c(1, 2, 2, 2))#定义绘图区域边界
create_beautiful_radarchart(
  data= data, caxislabels = c(-1, -0.5, 0, 0.5, 1),
  color= c("#00AFBB", "#E7B800", "#FC4E07"))
# 添加图例
legend(
  x ="bottom", legend = rownames(data[-c(1,2),]), horiz = TRUE,
  bty ="n", pch = 20 , col = c("#00AFBB", "#E7B800","#FC4E07"),
  text.col = "black", cex = 1, pt.cex = 1.5)
par(op)

#富集分析####
library(tidyverse)
library(org.Hs.eg.db)
library(clusterProfiler)
# DEG <- as.data.frame(res)%>% 
#   arrange(padj) %>% 
#   dplyr::filter(abs(log2FoldChange) > 1, padj < 0.05)
DEG=deg[gene_com,]

DEG=data.frame(gene_com)
rownames(DEG)=gene_com
colnames(DEG)="gene"

DEG <- DEG %>% rownames_to_column("Gene")

genelist <- bitr(DEG$Gene, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')
DEG <- inner_join(DEG,genelist,by=c("Gene"="SYMBOL"))
gene=DEG$Gene
#区分上下调基因
DEG=deg[gene_com,]
gene_up=gene[DEG$logFC>0]
gene_down=gene[DEG$logFC<0]
#修改gene为up或者down进行富集
gene=gene_up
#GO分析
ego <- enrichGO(gene = DEG$ENTREZID,
                OrgDb = org.Hs.eg.db, 
                ont = "all",
                pAdjustMethod = "BH",
                minGSSize = 1,
                pvalueCutoff =0.05, 
                qvalueCutoff =0.05,
                readable = TRUE)

ego_res <- ego@result
GO=ego_res


library(circlize)#GO的ID不用自己找，其他的需要自己总结
#ko.id=GO[,3]
ko.id=rownames(GO)
category=GO$ONTOLOGY
gene_num.min=0
#根据自己数据改max
gene_num.max=str_split_fixed(GO$GeneRatio,pattern = "/",2)[,2]
gene_num.rich=GO$Count
log10.p=-log10(GO$p.adjust)
up.regulated=rep(0,nrow(GO))
for (i in 1:nrow(GO)){
  up.regulated[i]=sum(unlist(str_split(GO$geneID[i] , pattern = "/")) %in% gene_up)
}

down.regulated=rep(0,nrow(GO))
for (i in 1:nrow(GO)){
  down.regulated[i]=sum(unlist(str_split(GO$geneID[i] , pattern = "/")) %in% gene_dpwn)
}
rich.factor=as.numeric(GO$Count)/as.numeric(
  str_split_fixed(GO$GeneRatio,pattern = "/",2)[,2])
rich.factor=(rich.factor-min(rich.factor))/(max(rich.factor)-min(rich.factor))
dat=data.frame(ko.id,category,gene_num.min,gene_num.max,gene_num.rich,
               log10.p,up.regulated,down.regulated,rich.factor)
dat[,3:9]=as.numeric(unlist(dat[,3:9]))
rownames(dat)=dat$ko.id
num=data.frame(table(dat$category))
num
n=5#展示的各个部分数量

dat=dat[c(1:n,(num[1,2]+1):(num[1,2]+5),
          (num[1,2]+num[2,2]+1):(num[1,2]+num[2,2]+5)),]
#整体布局
circos.par(gap.degree = 2, start.degree = 90)
##第一圈，绘制 ko
plot_data <- dat[c('ko.id', 'gene_num.min', 'gene_num.max')]  #选择作图数据集，定义了 ko 区块的基因总数量范围
ko_color <- c(rep('#F7CC13', 15), rep('#954572', 2), rep('#0796E0', 3))  #定义分组颜色

circos.genomicInitialize(plot_data, plotType = NULL, major.by = 1)  #一个总布局
circos.track(
  ylim = c(0, 1), track.height = 0.05, bg.border = NA, 
  bg.col = ko_color,  #圈图的高度、颜色等设置
  panel.fun = function(x, y) {
    ylim = get.cell.meta.data('ycenter')  #ylim、xlim 用于指定 ko id 文字标签添加的合适坐标
    xlim = get.cell.meta.data('xcenter')
    sector.name = get.cell.meta.data('sector.index')  #sector.name 用于提取 ko id 名称
    circos.axis(h = 'top', labels.cex = 0.6, major.tick.percentage = 0.6, labels.niceFacing = FALSE)  #绘制外周的刻度线
    circos.text(xlim, ylim, sector.name, cex = 0.6, niceFacing = FALSE)  #将 ko id 文字标签添加在图中指定位置处
  } )
##第二圈，绘制富集的基因和富集 p 值
plot_data <- dat[c('ko.id', 'gene_num.min', 'gene_num.rich', 'log10.p')]  #选择作图数据集，包括富集基因数量以及 p 值等信息
label_data <- dat['gene_num.rich']  #标签数据集，仅便于作图时添加相应的文字标识用
p_max <- round(max(dat$'log10.p')) + 1  #定义一个 p 值的极值，以方便后续作图
colorsChoice <- colorRampPalette(c('#FF906F', '#861D30'))  #这两句用于定义 p 值的渐变颜色
color_assign <- colorRamp2(breaks = 0:p_max, col = colorsChoice(p_max + 1))

circos.genomicTrackPlotRegion(
  plot_data, track.height = 0.08, bg.border = NA, stack = TRUE,  #圈图的高度、颜色等设置
  panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value, col = color_assign(value[[1]]), border = NA, ...)  #区块的长度反映了富集基因的数量，颜色与 p 值有关
    ylim = get.cell.meta.data('ycenter')  #同上文，ylim、xlim、sector.name 等用于指定文字标签（富集基因数量）添加的合适坐标
    xlim = label_data[get.cell.meta.data('sector.index'),1] / 2
    sector.name = label_data[get.cell.meta.data('sector.index'),1]
    circos.text(xlim, ylim, sector.name, cex = 0.6, niceFacing = FALSE)  #将文字标签添（富集基因数量）加在图中指定位置处
  } )

##第三圈，绘制上下调基因
#首先基于表格中上下调基因的数量，计算它们的占比
dat$all.regulated <- dat$up.regulated + dat$down.regulated
dat$up.proportion <- dat$up.regulated / dat$all.regulated
dat$down.proportion <- dat$down.regulated / dat$all.regulated

#随后，根据上下调基因的相对比例，分别计算它们在作图时的“区块坐标”和“长度”
dat$up <- dat$up.proportion * dat$gene_num.max
plot_data_up <- dat[c('ko.id', 'gene_num.min', 'up')]
names(plot_data_up) <- c('ko.id', 'start', 'end')
plot_data_up$type <- 1  #分配 1 指代上调基因

dat$down <- dat$down.proportion * dat$gene_num.max + dat$up
plot_data_down <- dat[c('ko.id', 'up', 'down')]
names(plot_data_down) <- c('ko.id', 'start', 'end')
plot_data_down$type <- 2  #分配 2 指代下调基因

#选择作图数据集（作图用）、标签数据集（添加相应的文字标识用），并分别为上下调基因赋值不同颜色
plot_data <- rbind(plot_data_up, plot_data_down)
label_data <- dat[c('up', 'down', 'up.regulated', 'down.regulated')]
color_assign <- colorRamp2(breaks = c(1, 2), col = c('#671166', '#7F8CCB'))

#继续绘制圈图
circos.genomicTrackPlotRegion(
  plot_data, track.height = 0.08, bg.border = NA, stack = TRUE,  #圈图的高度、颜色等设置
  panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value, col = color_assign(value[[1]]), border = NA, ...)  #这里紫色代表上调基因，蓝色代表下调基因，区块的长度反映了上下调基因的相对占比
    ylim = get.cell.meta.data('cell.bottom.radius') - 0.5  #同上文，ylim、xlim、sector.name 等用于指定文字标签（上调基因数量）添加的合适坐标
    xlim = label_data[get.cell.meta.data('sector.index'),1] / 2
    sector.name = label_data[get.cell.meta.data('sector.index'),3]
    circos.text(xlim, ylim, sector.name, cex = 0.5, niceFacing = FALSE)  #将文字标签（上调基因数量）添加在图中指定位置处
    xlim = (label_data[get.cell.meta.data('sector.index'),2]+label_data[get.cell.meta.data('sector.index'),1]) / 2
    sector.name = label_data[get.cell.meta.data('sector.index'),4]
    circos.text(xlim, ylim, sector.name, cex = 0.5, niceFacing = FALSE)  #类似的操作，将下调基因数量的标签也添加在图中
  } )

##第四圈，绘制富集因子

plot_data <- dat[c('ko.id', 'gene_num.min', 'gene_num.max', 'rich.factor')]  #选择作图数据集，包含富集因子列
label_data <- dat['category']  #将通路的分类信息提取出，和下一句一起，便于作图时按分组分配颜色

color_assign <- c('CC' = '#F7CC13', 'BP' = '#954572', 'MF' = '#0796E0')

circos.genomicTrack(
  plot_data, ylim = c(0, 1), track.height = 0.2, track.margin=c(0.05,0.05),
  bg.col = 'gray95', bg.border = NA,  #圈图的高度、颜色等设置
  panel.fun = function(region, value, ...) {
    sector.name = get.cell.meta.data('sector.index')  #sector.name 用于提取 ko id 名称，并添加在下一句中匹配 ko 对应的高级分类，以分配颜色
    circos.genomicRect(region, value, col = color_assign[label_data[sector.name,1]], border = NA, ytop.column = 1, ybottom = 0, ...)  #绘制矩形区块，高度代表富集因子数值，颜色代表 ko 的分类
    circos.lines(c(0, max(region)), c(0.5, 0.5), col = 'gray', lwd = 0.3)  #可选在富集因子等于 0.5 的位置处添加一个灰线
  } )

##使用 ComplexHeatmap 包绘制图例
library(ComplexHeatmap)

category_legend <- Legend(
  labels = c('CC', 'BP', 'MF'),
  type = 'points', pch = NA, background = c('#F7CC13', '#954572', '#0796E0'), 
  labels_gp = gpar(fontsize = 8), grid_height = unit(0.5, 'cm'), grid_width = unit(0.5, 'cm'))

updown_legend <- Legend(
  labels = c('Up-regulated', 'Down-regulated'), 
  type = 'points', pch = NA, background = c('#671166', '#7F8CCB'), 
  labels_gp = gpar(fontsize = 8), grid_height = unit(0.5, 'cm'), grid_width = unit(0.5, 'cm'))

pvalue_legend <- Legend(
  col_fun = colorRamp2(round(seq(0, p_max, length.out = 6), 0), 
                       colorRampPalette(c('#FF906F', '#861D30'))(6)),
  legend_height = unit(3, 'cm'), labels_gp = gpar(fontsize = 8), 
  title_gp = gpar(fontsize = 9), title_position = 'topleft', title = '-Log10(Pvalue)')

lgd_list_vertical <- packLegend(category_legend, updown_legend, pvalue_legend)
grid.draw(lgd_list_vertical)


#WGCNA####
library("tidyverse")
library("WGCNA") 
datExpr0=mtab_data[,gene_wgcna]
sampleTree = hclust(dist(datExpr0), method = "average")
#datExpr0行为样本，列为表达
# 画图
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree)
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
###剪切线
abline(h = 75, col = "red")
###删除剪切线以下的样品
clust = cutreeStatic(sampleTree, cutHeight = 75, minSize = 10)
table(clust)
keepSamples = (clust==1)
datExpr0 = datExpr0[keepSamples, ]
dev.off()
# 重新聚类
sampleTree2 = hclust(dist(datExpr0), method = "average")
plot(sampleTree2)

# 记录基因和样本数，方便后续可视化
nGenes = ncol(datExpr0)#基因数
nSamples = nrow(datExpr0)#样本数



# 构建网络，识别模块
# power值散点图
enableWGCNAThreads()   #多线程工作
powers = c(1:20)       #幂指数范围1:20
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)

par(mfrow = c(1,2))
cex1 = 0.9
###拟合指数与power值散点图
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red") #可以修改
###平均连通性与power值散点图
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

##6.2 邻接矩阵转换
sft #查看最佳power值
softPower =sft$powerEstimate #最佳power值
softPower = as.numeric(softPower)
adjacency = adjacency(datExpr0, power = softPower)

##6.3 TOM矩阵
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM


# 基因聚类
geneTree = hclust(as.dist(dissTOM), method = "average");
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)

# 动态剪切模块识别
minModuleSize = 30      #模块基因数目
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)

dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

# 相似模块聚类
MEList = moduleEigengenes(datExpr0, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average")
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

MEDissThres = 0.1 #剪切高度可修改
abline(h=MEDissThres, col = "red")

###相似模块合并
merge = mergeCloseModules(datExpr0, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs
plotDendroAndColors(geneTree, mergedColors,"Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

moduleColors = mergedColors
table(moduleColors)
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs
dev.off()
# 整理评分信息
clinical=cbind(score,result)

clinical <- clinical[rownames(datExpr0),]
identical(rownames(clinical),rownames(datExpr0))
# 查看评分信息
head(clinical)
# 对表达矩阵进行预处理
datTraits = as.data.frame(do.call(cbind,lapply(clinical, as.numeric)))
rownames(datTraits) = rownames(clinical)

# 对样本进行聚类
sampleTree2 = hclust(dist(datExpr0), method = "average")

# 将临床信息转换为颜色，白色表示低，红色表示高，灰色表示缺失
traitColors = numbers2colors(datTraits, signed = FALSE)

# 样本聚类图与样本性状热图
plotDendroAndColors(sampleTree2, 
                    traitColors,
                    groupLabels = names(datTraits), 
                    main = "Sample dendrogram and trait heatmap")
dev.off()

#### 网络的分析
# 对模块特征矩阵进行排序
MEs=orderMEs(MEs)
#计算模型特征矩阵和样本信息矩阵的相关度。
moduleTraitCor=cor(MEs, datTraits, use="p")

moduleTraitPvalue=corPvalueStudent(moduleTraitCor, nSamples)


#使用labeledHeatmap()将上述相关矩阵和p值可视化。
textMatrix=paste(signif(moduleTraitCor,2),"\n(",signif(moduleTraitPvalue,1),")",sep="")
dim(textMatrix)=dim(moduleTraitCor)
# 基因模块与临床信息相关性图
labeledHeatmap(Matrix=moduleTraitCor,#模块和表型的相关性矩阵，这个参数最重要，其他可以不变
               xLabels=colnames(datTraits),
               yLabels=names(MEs),
               ySymbols=names(MEs),
               colorLabels=FALSE,
               colors=blueWhiteRed(50),
               textMatrix=textMatrix,
               setStdMargins=FALSE,
               cex.text=0.5,
               cex.lab=0.5,
               zlim=c(-1,1),
               main=paste("Module-trait relationships"))
dev.off()

# 不同模块与基因性状的具体分析
##矩阵一
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"))
####看一下目的基因和哪个模块相关性最高
a <- geneModuleMembership
a <- a %>% rownames_to_column()

MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

##矩阵二
traitNames=names(datTraits)
geneTraitSignificance = as.data.frame(cor(datExpr0, datTraits, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", traitNames, sep="")
names(GSPvalue) = paste("p.GS.", traitNames, sep="")

##批量输出性状和模块散点图
for (trait in traitNames){
  traitColumn=match(trait,traitNames)  
  for (module in modNames){
    column = match(module, modNames)
    moduleGenes = moduleColors==module
    if (nrow(geneModuleMembership[moduleGenes,]) > 1){
      outPdf=paste(trait, "_", module,".pdf",sep="")
      pdf(file=outPdf,width=7,height=7)
      par(mfrow = c(1,1))
      verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                         abs(geneTraitSignificance[moduleGenes, traitColumn]),
                         xlab = paste("Module Membership in", module, "module"),
                         ylab = paste("Gene significance for ",trait),
                         main = paste("Module membership vs. gene significance\n"),
                         cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
      abline(v=0.8,h=0.5,col="red")
      dev.off()
    }
  }
}

#10. 输出每个模块的基因
for (mod in 1:nrow(table(moduleColors)))
{  
  modules = names(table(moduleColors))[mod]
  probes = colnames(datExpr0)
  inModule = (moduleColors == modules)
  modGenes = probes[inModule]
  write.table(modGenes, file =paste0(modules,".txt"),sep="\t",row.names=F,col.names=F,quote=F)
}

#共识聚类####
set.seed(100)
library(tidyverse)
library(ConsensusClusterPlus)

exp=data.frame(t(apply(mtab_data, 2, as.numeric)))
colnames(exp)=rownames(mtab_data)
d=as.matrix(exp)
gene=c("CCR7","CXCR6","CD8A","IL7R","CCL19","IL18","TNF","FLT3","XCL1","VCAM1")
d <- as.matrix(exp[as.matrix(intersect(t(gene),rownames(exp))),])
mads=apply(d,1,mad)
d=d[rev(order(mads))[1:as.numeric(nrow(d))],] #5是基因个数 
d = sweep(d,1, apply(d,1,median,na.rm=T))

title=("C:/Users/28978/Desktop/生信图片") ##文件夹输出图片的位置
set.seed(1) #我发现设不设置种子都一样
results = ConsensusClusterPlus(d,maxK=5,reps=50,pItem=0.8,pFeature=1,
                               title=title,clusterAlg="hc",
                               distance="pearson",seed=1,plot="pdf")

icl = calcICL(results,title=title,plot="pdf") ##画另一组图片

group<-results[[2]][["consensusClass"]]
group<-as.data.frame(group)

group$group <- factor(group$group,levels=c(1,2))


surv=exp_tar[,c(1,2)]
exp_gene <- exp[as.matrix(intersect(t(gene),rownames(exp))),]
surv=cbind(surv,group)
class(surv$group)
# 绘制ConsensusClusterPlus后的热图
library(pheatmap)
group <- group %>% 
  rownames_to_column("sample")
annotation <- group %>% arrange(group) %>% column_to_rownames("sample")
a <- group %>% arrange(group) %>% mutate(sample=substring(.$sample,1,12))
b <- t(exp_gene) %>% 
  as.data.frame() %>% 
  rownames_to_column("sample") %>% 
  mutate(sample=substring(.$sample,1,12))
c <- inner_join(a,b,"sample") %>% .[,-2] %>% column_to_rownames("sample") %>% t(.)
pheatmap(c,annotation = annotation,
         cluster_cols = F,fontsize=5,fontsize_row=5,
         scale="row",show_colnames=F,
         fontsize_col=3)
pheatmap(c,annotation = annotation,
         annotation_colors = list(group = c("1" ="#01468b","2"= "#ee0000")),
         cluster_cols = F,fontsize=5,fontsize_row=5,
         scale="row",show_colnames=F,cluster_row = F,
         fontsize_col=3)
dev.off()
#生存分析
library(survival)
surv=cbind(surv,group)
surv$futime=surv$futime/30
fitd <- survdiff(Surv(futime, fustate) ~ group,
                 data      = surv,
                 na.action = na.exclude)
pValue <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)

#2.2 拟合生存曲线
fit <- survfit(Surv(futime, fustate)~ group, data = surv)
summary(fit)
library(survminer)
p.lab <- paste0("P", ifelse(pValue < 0.001, " < 0.001", paste0(" = ",round(pValue, 3))))
ggsurvplot(fit,
           data = surv,
           pval = p.lab,
           conf.int = TRUE, # 显示置信区间
           risk.table = TRUE, # 显示风险表
           risk.table.col = "strata",
           palette = "jco", # 配色采用jco
           legend.labs = c("Cluster 1", "Cluster 2"), # 图例
           size = 1,
           xlim = c(0,120), # x轴长度，一般为0-10年
           break.time.by = 20, # x轴步长为20个月
           legend.title = "",
           surv.median.line = "hv", # 限制垂直和水平的中位生存
           ylab = "Survival probability (%)", # 修改y轴标签
           xlab = "Time (Months)", # 修改x轴标签
           ncensor.plot = TRUE, # 显示删失图块
           ncensor.plot.height = 0.25,
           risk.table.y.text = FALSE)
#分组
gene_target=read.table( "topgene.txt",sep= "\t",header=T,check.names=F)[1:15,1]
exp=cbind(group,test_exp)
exp=exp[c(gene_target[1:5],c("sample","group"))]
exp$group=ifelse(exp$group==1,"cluster1","cluster2")
exp_1=exp
colnames(exp)=c("Expression","Expression","Expression",
                "Expression","Expression","sample","group")
d_n=rbind(exp[,c(1,7)],exp[,c(2,7)],exp[,c(3,7)],exp[,c(4,7)],exp[,c(5,7)])
k=data.frame()
for (i in 1:5){
  k=rbind(k,data.frame(rep(gene_target[i],85)))
}
colnames(k)=c("Gene")
d_n=cbind(k,d_n)
library(ggpubr)  
p <- ggboxplot(d_n, x = "Gene", y = "Expression",
               color = "group", palette = "jco",
               add = "jitter")
p <- p + stat_compare_means(aes(group = group),label = "p.signif")   #  Add p-value
p #肿瘤分级 

#PCA
library(ggplot2)
library(factoextra)
library(FactoMineR)
load("325cgga.Rda")
load("group_325(85)_共识.Rda")
data <- cbind(group,test_exp)
df= k
iris.pca<- PCA(df, graph = FALSE)
fviz_pca_ind(iris.pca,
             geom.ind = "point", # show points only (nbut not "text")
             pointsize =3,
             pointshape = 21,
             fill.ind = factor(k$group,levels = c(1,2)), # color by groups
             palette = "lacent",# c("#00AFBB", "#E7B800", "#FC4E07")
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups",
             title="")+
  theme_bw() +
  theme(text=element_text(size=14,face="plain",color="black"),
        axis.title=element_text(size=16,face="plain",color="black"),
        axis.text = element_text(size=14,face="plain",color="black"),
        legend.title = element_text(size=16,face="plain",color="black"),
        legend.text = element_text(size=14,face="plain",color="black"),
        legend.background = element_blank(),
        legend.position=c(0.9,0.1)
  )


#桑葚图####
library(ggalluvial)
library(ggplot2)
library(dplyr)
#示例数据
#            ID       Tert_status   Mycn_status
# NBC-1     NBC-1        efjs     amplified
# NBC-100 NBC-100        efjs non amplified

rt=clinical
corLodes=to_lodes_form(rt, axes = 1:ncol(rt), id = "Cohort")

color <- c("#F9837B","#F5886C","#F18C5A","#EC9042","#E79419","#E09719",
           "#DA9C19","#D49F19","#CCA319","#C4A619",
           "#BBA919","#B1AC19","#A8B019","#9CB219","#8FB519","#81B819",
           "#70BA19","#59BC19","#30BE19","#19C043",
           "#19C25A","#19C36B","#19C47A","#19C587","#19C694","#19C6A0",
           "#19C7AB","#19C6B6","#19C6C1","#19C5CA",
           "#19C4D4","#19C3DC","#19C0E4","#19BEEC","#19BAF2","#19B7F9",
           "#19B3FD","#19AEFF","#58A9FF","#7AA4FF",
           "#939EFF","#A798FF","#B892FF","#FF79A3","#FF76AF","#FF73BA",
           "#FF71C4","#FF70CE","#FC70D8","#F971E0",
           "#F474E8","#EE78F0","#E77BF7","#DE81FC","#D386FF","#C68CFF",
           "#A953FF","#FF4374")
mycol=rep(color[1:nrow(rt)],nrow(rt))

mycol <- rep(c("#F9837B","#F5886C","#F18C5A","#EC9042","#E79419","#E09719","#DA9C19","#D49F19","#CCA319","#C4A619",
               "#BBA919","#B1AC19","#A8B019","#9CB219","#8FB519","#81B819","#70BA19","#59BC19","#30BE19","#19C043",
               "#19C25A","#19C36B","#19C47A","#19C587","#19C694","#19C6A0","#19C7AB","#19C6B6","#19C6C1","#19C5CA",
               "#19C4D4","#19C3DC","#19C0E4","#19BEEC","#19BAF2","#19B7F9","#19B3FD","#19AEFF","#58A9FF","#7AA4FF",
               "#939EFF","#A798FF","#B892FF","#FF79A3","#FF76AF","#FF73BA","#FF71C4","#FF70CE","#FC70D8","#F971E0",
               "#F474E8","#EE78F0","#E77BF7","#DE81FC","#D386FF","#C68CFF","#A953FF","#FF4374"),58)

ggplot(corLodes, aes(x = x, stratum = stratum, alluvium = Cohort,fill = stratum, label = stratum)) +
  scale_x_discrete(expand = c(0, 0)) +  
  #用aes.flow控制线调颜色，forward说明颜色和前面一致，backward说明与后面一致
  geom_flow(width = 2/10,aes.flow = "forward",show.legend = TRUE) +  #画流动图
  geom_stratum(alpha = .9,width = 5/10, # 格子宽度
               linetype=1,size=0.4, # 格子的边框线
               color = "#696969",inherit.aes=TRUE) + #画冲击图
  scale_fill_manual(values = mycol) +
  #size = 2代表基因名字大小
  geom_text(stat = "stratum", size = 4,color="black") +
  xlab("") + ylab("") + theme_bw() + 
  theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text.y = element_blank()) + #ȥ????????
  theme(panel.grid =element_blank()) + 
  theme(panel.border = element_blank()) + 
  ggtitle("") + guides(fill = FALSE) 

#PCA####
library(ggplot2)
library(factoextra)
library(FactoMineR)
load("325cgga.Rda")
load("group_325(85)_共识.Rda")
data <- cbind(group,test_exp)
library(factoextra)
library(FactoMineR)
iris= mtab_data

iris$Stage=ifelse(iris$Stage>=3,">=3","<3")

df <- iris[,c(gene)]
target=iris[,"Stage"]
iris.pca<- PCA(df, graph = FALSE)
fviz_pca_ind(iris.pca,
             geom.ind = "point", # show points only (nbut not "text")
             pointsize =3,
             pointshape = 21,
             fill.ind = target, # color by groups
             palette = "lacent",# c("#00AFBB", "#E7B800", "#FC4E07")
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups",
             title="")+
  theme_bw() +
  theme(text=element_text(size=14,face="plain",color="black"),
        axis.title=element_text(size=16,face="plain",color="black"),
        axis.text = element_text(size=14,face="plain",color="black"),
        legend.title = element_text(size=16,face="plain",color="black"),
        legend.text = element_text(size=14,face="plain",color="black"),
        legend.background = element_blank(),
        legend.position=c(0.9,0.1)
  )

#相关性弦图####
library(circlize)
#c("CNR1","GNAI1","GLDC","ABCC4")
ICs=c("CD48","CD40","CD27","CTLA4","LAG3","CNR1","GNAI1","GLDC","ABCC4")
data=apply(target_data[,ICs], 2,as.numeric)
r=cor(data)
gene=rep(c("CNR1","GNAI1","GLDC","ABCC4"),5)
ICS=c(rep("CD48",4),rep("CD40",4),rep("CD27",4),rep("CTLA4",4),rep("LAG3",4))
R=c(t(r[6:9,"CD48"]),t(r[6:9,"CD40"]),t(r[6:9,"CD27"]),t(r[6:9,"CTLA4"]),t(r[6:9,"LAG3"]))
df=data.frame(gene,ICS,R)
# 读取数据
df1 <- read.delim("https://www.bioladder.cn/shiny/zyp/bioladder2/model/bioladder1/Chord/demoData.txt")
#颜色
col_fun = colorRamp2(range(df$R), c('#0000FF','#FF0000'), transparency = 0.3)
# 绘图
chordDiagram(
  x = df,
  col = col_fun,
  grid.col = hcl.colors(9),                     # 颜色方案，数字向量11要和实际的数据相符
  directional = 1,                               # 箭头方向。选项有1,0,-1
  direction.type = c("arrows", "diffHeight"),    # 线条两端的形状
  diffHeight = -0.02,                            # 线条两端距离边缘的距离差
  annotationTrack = c("name", "grid", "axis"),   # 都绘制哪些内容，name标签；grid边缘形状；axis刻度
  annotationTrackHeight = c(0.05, 0.08),         # 标签距离图形的距离; 环形边缘的宽度
  link.arr.type = "big.arrow",                   # 形状"curved", "triangle", "circle", "ellipse".
  link.sort = TRUE,                              # 内部排序
  link.largest.ontop = F,                     # 控制添加链接的顺序，是否基于绝对值?
  transparency = 0.25                            # 线条透明度
)
# 更多参数?chordDiagram查看
# 附录 拓展
# 1.修改数据标签和坐标轴
#图列
grid.draw(Legend(
  col_fun = colorRamp2(round(seq(-1, 1, length.out = 6), 0), 
                       colorRampPalette(c('#0000FF','#FF0000'))(6)),
  legend_height = unit(3, 'cm'), labels_gp = gpar(fontsize = 8), 
  title_gp = gpar(fontsize = 9), title_position = 'topleft', title = 'R')
)
#
circos.trackPlotRegion(
  track.index = 1,
  bg.border = NA,
  panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    sector.index = get.cell.meta.data("sector.index")
    # 添加数据标签
    circos.text(
      x = mean(xlim),
      y = 3.2,
      labels = sector.index,
      facing = "bending",
      cex = 1
    )
    # 添加坐标轴
    circos.axis(
      h = "top",
      major.at = seq(from = 0, to = xlim[2], by = ifelse(test = xlim[2]>10, yes = 2, no = 1)),
      minor.ticks = 1,
      major.tick.percentage = 0.5,
      labels.niceFacing = FALSE)
  })



# 2.调节边距，起始角度，间隔角度等参数

circos.par(start.degree = 90,
           gap.degree = 4,
           track.margin = c(-0.1, 0.1), 
           points.overflow.warning = FALSE)
circos.clear()
