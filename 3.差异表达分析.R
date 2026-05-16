rm(list = ls()) ##一键清空
options(stringsAsFactors = F)

# 加载所需的R包
library(ggplot2)  # 用于数据可视化
library(limma)    # 用于线性模型拟合和差异表达分析
library(pheatmap) # 用于绘制热图

# 设置过滤阈值
logFCfilter = 0.5         # 设置logFC过滤阈值
adjPfilter = 0.05       # 设置调整后p值的过滤阈值

# 指定输入文件路径
expFile = "merge.txt"  # 基因表达矩阵文件
conFile = "sample1.txt"     # 对照组样本文件
treatFile = "sample2.txt"   # 实验组样本文件

# 设置工作目录
setwd("D:\\项目研究--复方药效\\6.数据集合并")  # 设置当前工作目录

# 读取基因表达矩阵文件并进行处理
rt = read.table(expFile, sep="\t", header=T, check.names=F)  # 读取表达矩阵
rt = as.matrix(rt)  # 转换为矩阵格式
rownames(rt) = rt[,1]  # 将第一列设为行名
exp = rt[,2:ncol(rt)]  # 去掉第一列基因名称，保留表达数据
dimnames = list(rownames(exp), colnames(exp))  # 设置行列名称
data = matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)  # 确保数据为数值型
data = avereps(data)  # 对重复基因取平均值
data = normalizeBetweenArrays(data)  # 标准化数据，使得样本间可比

# 读取样本信息文件
sample1 = read.table(conFile, sep="\t", header=F, check.names=F)  # 读取对照组样本信息
sample2 = read.table(treatFile, sep="\t", header=F, check.names=F)  # 读取实验组样本信息
conData = data[,as.vector(sample1[,1])]  # 获取对照组的表达数据
treatData = data[,as.vector(sample2[,1])]  # 获取实验组的表达数据
rt = cbind(conData, treatData)  # 合并对照组和实验组数据
conNum = ncol(conData)  # 获取对照组样本数量
treatNum = ncol(treatData)  # 获取实验组样本数量

# 进行差异表达分析
Type = c(rep("con", conNum), rep("treat", treatNum))  # 生成样本分组标签
design <- model.matrix(~0+factor(Type))  # 创建设计矩阵
colnames(design) <- c("con","treat")  # 为设计矩阵的列命名
fit <- lmFit(rt, design)  # 线性模型拟合
cont.matrix <- makeContrasts(treat-con, levels=design)  # 设置对比矩阵
fit2 <- contrasts.fit(fit, cont.matrix)  # 应用对比
fit2 <- eBayes(fit2)  # 经验贝叶斯方法计算统计量

# 获取所有基因的差异分析结果
allDiff = topTable(fit2, adjust='fdr', number=200000)  # 获取差异分析结果
write.table(allDiff, file="GEO_all.xls", sep="\t", quote=F)  # 保存结果到文件

# 筛选显著差异表达基因
diffSig = allDiff[with(allDiff, (abs(logFC) > logFCfilter & adj.P.Val < adjPfilter)), ]  # 根据logFC和p值筛选
diffSigOut = rbind(id=colnames(diffSig), diffSig)  # 添加列名
write.table(diffSigOut, file="GEO_diff.xls", sep="\t", quote=F, col.names=F)  # 保存差异表达基因结果
write.table(diffSigOut, file="GEO_diff.txt", sep="\t", quote=F, col.names=F)

# 绘制热图
geneNum = 50  # 热图中显示的最大基因数
diffSig = diffSig[order(as.numeric(as.vector(diffSig$logFC))), ]  # 按logFC排序
diffGeneName = as.vector(rownames(diffSig))  # 提取基因名称
diffLength = length(diffGeneName)  # 获取差异基因的总数
hmGene = c()  # 用于存储热图中的基因
if(diffLength > (geneNum * 2)) {
  hmGene = diffGeneName[c(1:geneNum, (diffLength-geneNum+1):diffLength)]  # 选择上下调显著的基因
} else {
  hmGene = diffGeneName  # 如果差异基因数较少，全部展示
}
hmExp = rt[hmGene, ]  # 提取热图基因的表达数据
Type = c(rep("N", conNum), rep("T", treatNum))  # 设置样本类型标签
names(Type) = colnames(rt)  # 为样本标签命名
Type = as.data.frame(Type)  # 转换为数据框格式

# 绘制并保存热图
pdf(file="GEO_heatmap.pdf", height=8, width=10)
pheatmap(hmExp, 
         annotation=Type, 
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),  # 颜色方案更新为蓝白红渐变
         cluster_cols =F,  # 不聚类列
         show_colnames = F,  # 不显示列名
         scale="row",  # 在行间进行标准化
         fontsize = 14,  # 设置字体大小
         fontsize_row=6,  # 设置行名字体大小
         fontsize_col=14)  # 设置列名字体大小
dev.off()

# 绘制火山图
Significant = ifelse((allDiff$adj.P.Val < adjPfilter & abs(allDiff$logFC) > logFCfilter), 
                     ifelse(allDiff$logFC > logFCfilter, "Up", "Down"), "Not")  # 标记显著性

p = ggplot(allDiff, aes(logFC, -log10(adj.P.Val))) +  # 绘制火山图
  geom_point(aes(col=Significant)) +  # 设置点的颜色
  scale_color_manual(values=c("navy", "grey", "firebrick3")) +  # 颜色方案更新为红灰蓝
  labs(title = "Volcano Plot") +  # 添加图表标题
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))  # 设置标题样式
p = p + theme_bw()  # 使用白色背景主题

# 保存火山图为PDF
pdf("GEO_vol.pdf", width=5.5, height=5)
print(p)
dev.off()