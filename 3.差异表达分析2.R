rm(list = ls()) ##一键清空
options(stringsAsFactors = F)

# 加载所需的R包
library(ggplot2)  
library(limma)    
library(pheatmap) 

# 设置过滤阈值
logFCfilter = 0.5         
adjPfilter = 0.05       

# 指定输入文件路径
expFile = "merged_after_batch_removal.txt"  

# 设置工作目录
setwd("E:/1paper/合并数据/BatchEffect  _output_20260312_141128")  

# 读取基因表达矩阵文件并进行处理
rt = read.table(expFile, sep="\t", header=T, check.names=F)  
rt = as.matrix(rt)  
rownames(rt) = rt[,1]  
exp = rt[,2:ncol(rt)]  
dimnames = list(rownames(exp), colnames(exp))  
data = matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)  
data = avereps(data)  
data = normalizeBetweenArrays(data)  

print("正在智能识别 Control(正常) 和 Treat(患病) 样本...")
# 【核心修复点】：改成识别你真实的列名特征！
all_samples = colnames(data)
# 抓取带有 "Control" 的列名（忽略大小写）
con_samples = all_samples[grep("Control", all_samples, ignore.case = TRUE)]   
# 抓取带有 "Treat" 或 "Disease" 的列名（忽略大小写，防止拼写差异）
treat_samples = all_samples[grep("Treat|Disease", all_samples, ignore.case = TRUE)] 

# 如果还是没找到，直接让程序停下来并告诉你列名到底叫什么
if(length(con_samples) == 0 | length(treat_samples) == 0) {
  print("🚨 警告：没有找到足够的 Control 或 Treat 样本！你矩阵的前5个列名是：")
  print(all_samples[1:5])
  stop("请检查列名特征！")
}

conData = data[, con_samples]     
treatData = data[, treat_samples] 
rt = cbind(conData, treatData)    
conNum = ncol(conData)            
treatNum = ncol(treatData)        
print(paste("识别成功！正常样本:", conNum, "个，患病样本:", treatNum, "个"))

print("正在进行差异表达分析...")
Type = c(rep("con", conNum), rep("treat", treatNum))  
design <- model.matrix(~0+factor(Type))  
colnames(design) <- c("con","treat")  
fit <- lmFit(rt, design)  
cont.matrix <- makeContrasts(treat-con, levels=design)  
fit2 <- contrasts.fit(fit, cont.matrix)  
fit2 <- eBayes(fit2)  

# 获取所有基因的差异分析结果并保存
allDiff = topTable(fit2, adjust='fdr', number=200000)  
write.table(allDiff, file="GEO_all.xls", sep="\t", quote=F)  

# 筛选显著差异表达基因
diffSig = allDiff[with(allDiff, (abs(logFC) > logFCfilter & adj.P.Val < adjPfilter)), ]  
diffSigOut = rbind(id=colnames(diffSig), diffSig)  
write.table(diffSigOut, file="GEO_diff.xls", sep="\t", quote=F, col.names=F)  
write.table(diffSigOut, file="GEO_diff.txt", sep="\t", quote=F, col.names=F)

print(paste("★★★ 恭喜！成功捞出", nrow(diffSig), "个显著差异基因！ ★★★"))
  

print("正在绘制热图...")
geneNum = 50  
diffSig = diffSig[order(as.numeric(as.vector(diffSig$logFC))), ]  
diffGeneName = as.vector(rownames(diffSig))  
diffLength = length(diffGeneName)  
hmGene = c()  
if(diffLength > (geneNum * 2)) {
  hmGene = diffGeneName[c(1:geneNum, (diffLength-geneNum+1):diffLength)]  
} else {
  hmGene = diffGeneName  
}
hmExp = rt[hmGene, ]  
Type = c(rep("N", conNum), rep("T", treatNum))  
names(Type) = colnames(rt)  
Type = as.data.frame(Type)  

# =========================================================
# 【新增的“救命代码”：干掉方差为0或带有 NA 的脏数据】
hmExp = na.omit(hmExp)  # 踢掉含有空值的行
hmExp = hmExp[apply(hmExp, 1, sd) > 0, ]  # 踢掉标准差为0的行
# =========================================================

pdf(file="GEO_heatmap.pdf", height=8, width=10)
pheatmap(hmExp, 
         annotation=Type, 
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),  
         cluster_cols =F,  
         show_colnames = F,  
         scale="row",  # 现在有安检了，按行标准化绝对不会报错！
         fontsize = 14,  
         fontsize_row=6,  
         fontsize_col=14)  
dev.off()

print("热图绘制成功！")

print("正在绘制热图...")
geneNum = 50  
diffSig = diffSig[order(as.numeric(as.vector(diffSig$logFC))), ]  
diffGeneName = as.vector(rownames(diffSig))  
diffLength = length(diffGeneName)  
hmGene = c()  
if(diffLength > (geneNum * 2)) {
  hmGene = diffGeneName[c(1:geneNum, (diffLength-geneNum+1):diffLength)]  
} else {
  hmGene = diffGeneName  
}
hmExp = rt[hmGene, ]  
Type = c(rep("N", conNum), rep("T", treatNum))  
names(Type) = colnames(rt)  
Type = as.data.frame(Type)  

pdf(file="GEO_heatmap.pdf", height=8, width=10)
pheatmap(hmExp, 
         annotation=Type, 
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),  
         cluster_cols =F,  
         show_colnames = F,  
         scale="row",  
         fontsize = 14,  
         fontsize_row=6,  
         fontsize_col=14)  
dev.off()


print("全部完成！请去文件夹验收热图、火山图和差异基因表格！")