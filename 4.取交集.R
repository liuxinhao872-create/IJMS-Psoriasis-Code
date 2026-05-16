rm(list = ls()) ##一键清空
options(stringsAsFactors = F)

# 加载VennDiagram包，用于绘制维恩图
library(VennDiagram)

# 设置当前工作目录，确保文件读写操作都在指定的文件夹中进行
setwd("D:\\项目研究--复方药效\\6.数据集合并")

# 获取当前目录下的所有文件名
files = dir()

# 筛选出以".txt"结尾的文件（通常是包含基因列表的文件）
files = grep(".txt$", files, value = TRUE)

# 初始化一个列表，用于存储每个文件中的基因信息
geneList = list()

# 遍历所有的txt文件，读取其中的基因信息并存储在geneList列表中
for(i in 1:length(files)) {
  inputFile = files[i]
  # 跳过名为 "interGenes.txt" 的文件（通常这个文件是用于存储交集结果的，避免干扰）
  if(inputFile == "interGenes.txt") {next}
  
  # 读取当前txt文件的内容
  rt = read.table(inputFile, header = FALSE, sep = "\t", check.names = FALSE)
  
  # 提取第一列（基因名）作为基因列表
  geneNames = as.vector(rt[, 1])
  
  # 去除基因名开头和结尾的空格
  geneNames = gsub("^ | $", "", geneNames)
  
  # 去除重复的基因名，保留唯一的基因
  uniqGene = unique(geneNames)
  
  # 从文件名中提取前缀作为列表的名字
  header = unlist(strsplit(inputFile, "\\.|\\-"))
  
  # 将当前文件的基因列表存储在geneList中，键名为文件名前缀
  geneList[[header[1]]] = uniqGene
  
  # 打印文件名前缀及基因列表的长度
  uniqLength = length(uniqGene)
  print(paste(header[1], uniqLength, sep = " "))
}

# 绘制维恩图，表示多个基因列表之间的交集情况
venn.plot = venn.diagram(geneList, filename = NULL, fill = rainbow(length(geneList)), scaled = FALSE)

# 将维恩图保存为PDF文件
pdf(file = "venn.pdf", width = 5, height = 5)
grid.draw(venn.plot)
dev.off()

# 计算多个基因列表的交集，即找到同时出现在所有文件中的基因
intersectGenes = Reduce(intersect, geneList)

# 将交集结果写入文件 "interGenes.txt"
write.table(file = "interGenes.txt", intersectGenes, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)








rm(list = ls()) ##一键清空
options(stringsAsFactors = F)

# 加载VennDiagram包，用于绘制维恩图
library(VennDiagram)

# 设置当前工作目录，确保文件读写操作都在指定的文件夹中进行
setwd("D:\\EI源刊研究\\项目研究--肝癌\\6.免疫微环境取交集")

# 获取当前目录下的所有文件名
files = dir()

# 筛选出以".txt"结尾的文件（通常是包含基因列表的文件）
files = grep("txt$", files, value = TRUE)

# 初始化一个列表，用于存储每个文件中的基因信息
geneList = list()

# 遍历所有的txt文件，读取其中的基因信息并存储在geneList列表中
for(i in 1:length(files)) {
  inputFile = files[i]
  # 跳过名为 "interGenes.txt" 的文件（通常这个文件是用于存储交集结果的，避免干扰）
  if(inputFile == "interGenes.txt") {next}
  
  # 读取当前txt文件的内容
  rt = read.table(inputFile, header = FALSE, sep = "\t", check.names = FALSE)
  
  # 提取第一列（基因名）作为基因列表
  geneNames = as.vector(rt[, 1])
  
  # 去除基因名开头和结尾的空格
  geneNames = gsub("^ | $", "", geneNames)
  
  # 去除重复的基因名，保留唯一的基因
  uniqGene = unique(geneNames)
  
  # 从文件名中提取前缀作为列表的名字
  header = unlist(strsplit(inputFile, "\\.|\\-"))
  
  # 将当前文件的基因列表存储在geneList中，键名为文件名前缀
  geneList[[header[1]]] = uniqGene
  
  # 打印文件名前缀及基因列表的长度
  uniqLength = length(uniqGene)
  print(paste(header[1], uniqLength, sep = " "))
}

# 绘制维恩图，表示多个基因列表之间的交集情况
venn.plot = venn.diagram(geneList, filename = NULL, fill = rainbow(length(geneList)), scaled = FALSE)

# 将维恩图保存为PDF文件
pdf(file = "venn.pdf", width = 5, height = 5)
grid.draw(venn.plot)
dev.off()

# 计算多个基因列表的交集，即找到同时出现在所有文件中的基因
intersectGenes = Reduce(intersect, geneList)

# 将交集结果写入文件 "interGenes.txt"
write.table(file = "interGenes.txt", intersectGenes, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)








rm(list = ls()) ##一键清空
options(stringsAsFactors = F)

library(VennDiagram)

# 设置当前工作目录，确保文件读写操作都在指定的文件夹中进行
setwd("D:\\EI源刊研究\\项目研究--肾乳头细胞癌\\6.免疫微环境取交集")

# 获取当前目录下的所有文件名
files = dir()

# 筛选出以".txt"结尾的文件（通常是包含基因列表的文件）
files = grep("txt$", files, value = TRUE)

# 初始化一个列表，用于存储每个文件中的基因信息
geneList = list()

# 遍历所有的txt文件，读取其中的基因信息并存储在geneList列表中
for(i in 1:length(files)) {
  inputFile = files[i]
  # 跳过名为 "interGenes.txt" 的文件（通常这个文件是用于存储交集结果的，避免干扰）
  if(inputFile == "interGenes.txt") {next}
  
  # 读取当前txt文件的内容
  rt = read.table(inputFile, header = FALSE, sep = "\t", check.names = FALSE)
  
  # 提取第一列（基因名）作为基因列表
  geneNames = as.vector(rt[, 1])
  
  # 去除基因名开头和结尾的空格
  geneNames = gsub("^ | $", "", geneNames)
  
  # 去除重复的基因名，保留唯一的基因
  uniqGene = unique(geneNames)
  
  # 从文件名中提取前缀作为列表的名字
  header = unlist(strsplit(inputFile, "\\.|\\-"))
  
  # 将当前文件的基因列表存储在geneList中，键名为文件名前缀
  geneList[[header[1]]] = uniqGene
  
  # 打印文件名前缀及基因列表的长度
  uniqLength = length(uniqGene)
  print(paste(header[1], uniqLength, sep = " "))
}

# 绘制维恩图，表示多个基因列表之间的交集情况
venn.plot = venn.diagram(geneList, filename = NULL, fill = rainbow(length(geneList)), scaled = FALSE)

# 将维恩图保存为PDF文件
pdf(file = "venn.pdf", width = 5, height = 5)
grid.draw(venn.plot)
dev.off()

# 计算多个基因列表的交集，即找到同时出现在所有文件中的基因
intersectGenes = Reduce(intersect, geneList)

# 将交集结果写入文件 "interGenes.txt"
write.table(file = "interGenes.txt", intersectGenes, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)






rm(list = ls()) ##一键清空
options(stringsAsFactors = F)

# 加载VennDiagram包，用于绘制维恩图
library(VennDiagram)

# 设置当前工作目录，确保文件读写操作都在指定的文件夹中进行
setwd("D:\\EI源刊研究\\项目研究--肺腺癌\\6.免疫微环境取交集")

# 获取当前目录下的所有文件名
files = dir()

# 筛选出以".txt"结尾的文件（通常是包含基因列表的文件）
files = grep("txt$", files, value = TRUE)

# 初始化一个列表，用于存储每个文件中的基因信息
geneList = list()

# 遍历所有的txt文件，读取其中的基因信息并存储在geneList列表中
for(i in 1:length(files)) {
  inputFile = files[i]
  # 跳过名为 "interGenes.txt" 的文件（通常这个文件是用于存储交集结果的，避免干扰）
  if(inputFile == "interGenes.txt") {next}
  
  # 读取当前txt文件的内容
  rt = read.table(inputFile, header = FALSE, sep = "\t", check.names = FALSE)
  
  # 提取第一列（基因名）作为基因列表
  geneNames = as.vector(rt[, 1])
  
  # 去除基因名开头和结尾的空格
  geneNames = gsub("^ | $", "", geneNames)
  
  # 去除重复的基因名，保留唯一的基因
  uniqGene = unique(geneNames)
  
  # 从文件名中提取前缀作为列表的名字
  header = unlist(strsplit(inputFile, "\\.|\\-"))
  
  # 将当前文件的基因列表存储在geneList中，键名为文件名前缀
  geneList[[header[1]]] = uniqGene
  
  # 打印文件名前缀及基因列表的长度
  uniqLength = length(uniqGene)
  print(paste(header[1], uniqLength, sep = " "))
}

# 绘制维恩图，表示多个基因列表之间的交集情况
venn.plot = venn.diagram(geneList, filename = NULL, fill = rainbow(length(geneList)), scaled = FALSE)

# 将维恩图保存为PDF文件
pdf(file = "venn.pdf", width = 5, height = 5)
grid.draw(venn.plot)
dev.off()

# 计算多个基因列表的交集，即找到同时出现在所有文件中的基因
intersectGenes = Reduce(intersect, geneList)

# 将交集结果写入文件 "interGenes.txt"
write.table(file = "interGenes.txt", intersectGenes, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)