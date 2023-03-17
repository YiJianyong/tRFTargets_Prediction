#--------------references-------------
new_ID <- as.data.frame(gsub(pattern = "^.*\\(",replacement = "",x = new_ID))
new_ID <- as.data.frame(gsub(pattern = "\\)",replacement = "",x = new_ID))
for(i in 6:24){
  ID <- as.vector(gene_id[,i])
  ID <- as.data.frame(gsub(pattern = "^.*\\(",replacement = "",x = ID))
  ID <- as.data.frame(gsub(pattern = "\\)",replacement = "",x = ID))
  new_ID <- cbind(new_ID,ID)
}

new_ID <- cbind(gene_id[,1], gene_id[ ,2],
                gene_id[,3], gene_id[ ,4],new_ID)

colnames(new_ID) <- c("tRF-5.Sequence","Rest.of.read.which.was.BLAST.searched.against.Ref-Seq.RNA",
                      "Number.of.reads","E-Value","most","V6","V7","V8","V9","V10","V11","V12","V13"
                      ,"V14","V15","V16","V17","V18","V19","V20","V21","V22","V23","V24")

save(new_ID, file = "new_trf5.RData")
write.csv(new_ID,file = "new_ID5.csv")


load("new_trf3.RData")
#写入新的csv文件中
write.csv(new_ID,file = "new_ID3.csv")

load("new_trf1.RData")
#写入新的csv文件中
write.csv(new_ID,file = "new_ID1.csv")

load("new_trf2.RData")
#写入新的csv文件中
write.csv(new_ID,file = "new_ID5.csv")


#算logFC
rm(list = ls())
library(openxlsx)
gene_info <- read.xlsx("tRF_data/GSE184690_profiles_1.xlsx",sheet = 1,startRow = 9)
gene <- read.csv(file = "tRF_data/GSE189510.csv")
gene_id <- read.xlsx("tRF_data/gene_info.xlsx")

 
change_t1 <- log(gene_info$`t-1` / gene_info$`t-NC1`)
change_t2 <- log(gene_info$`t-2` / gene_info$`t-NC2`)
change_t3 <- log(gene_info$`t-3` / gene_info$`t-NC3`)

change_data <- cbind(as.data.frame(change_t1),as.data.frame(change_t2),as.data.frame(change_t3))

#symbol ID 和其他类型ID之间的转换
rm(list = ls())
#加载数据
load("new_trf1.RData")
train_data1 <- new_ID
load("new_trf2.RData")
train_data2 <- new_ID
load("new_trf3.RData")
train_data3 <- new_ID

#处理数据：ID转换  
library("clusterProfiler")
library("org.Hs.eg.db")

EM_id <- bitr(train_data1[,6], fromType = "REFSEQ", toType = "ENTREZID",
              OrgDb = Org.Hs.eg.db)





#--------读取mrna序列------------
mrna <- read.table(file = "tRF_data/gene_info/m3.txt")
mrna <- read.table(file = "tRF_data/gene_info/m5.txt")
mrna <- read.table(file = "tRF_data/gene_info/cds.txt")
mrna_id <- as.data.frame(mrna[seq(1,nrow(mrna),2),])
mrna_seq <- as.data.frame(mrna[seq(2,nrow(mrna),2),])

utr_3 <- cbind(mrna_id, mrna_seq)
utr_5 <- cbind(mrna_id, mrna_seq)
cds <- cbind(mrna_id, mrna_seq)

save(utr_3,file = "utr_3.RData")
save(utr_5,file = "utr_5.RData")
save(cds,file = "cds.RData")

load("utr_3.RData")
load("utr_5.RData")
load("cds.RData")

#正则替换ID号
a = ">ashahshhas|dadadadad|dsdddsdsd|sdsdsd"
b <- gsub(pattern = ">",replacement = "",x = a)
b <- gsub(pattern = "\\|.*",replacement = "",x = b)

mrna_ID <- as.data.frame(utr_3[,1])
for(i in 1:nrow(mrna_ID)){
  mrna_ID[i,] <- gsub(pattern = ">",replacement = "",x = mrna_ID[i,])
  mrna_ID[i,] <- gsub(pattern = "\\|.*",replacement = "",x = mrna_ID[i,])
}


mrna_ID5 <- as.data.frame(utr_5[,1])
for(i in 1:nrow(mrna_ID5)){
  mrna_ID5[i,] <- gsub(pattern = ">",replacement = "",x = mrna_ID5[i,])
  mrna_ID5[i,] <- gsub(pattern = "\\|.*",replacement = "",x = mrna_ID5[i,])
}

mrna_ID_cds <- as.data.frame(cds[,1])
for(i in 1:nrow(mrna_ID_cds)){
  mrna_ID_cds[i,] <- gsub(pattern = ">",replacement = "",x = mrna_ID_cds[i,])
  mrna_ID_cds[i,] <- gsub(pattern = "\\|.*",replacement = "",x = mrna_ID_cds[i,])
}

save(mrna_ID,file = "mrna_ID3.RData")
save(mrna_ID5,file = "mrna_ID5.RData")
save(mrna_ID_cds,file = "mrna_ID_cds.RData")
#gene ID转换
library(stats4)
library(BiocGenerics)
library(parallel)

library("AnnotationDbi")

library("org.Hs.eg.db")

mrna3_vector <- as.vector(t(mrna_ID))
mrna_3_geneid <- mapIds(org.Hs.eg.db,keys=mrna3_vector,column="SYMBOL",keytype="ENSEMBL",multiVals="first")
mrna_3_geneid <- as.data.frame(mrna_3_geneid)

mrna5_vector <- as.vector(t(mrna_ID5))
mrna_5_geneid <- mapIds(org.Hs.eg.db,keys=mrna5_vector,column="SYMBOL",keytype="ENSEMBL",multiVals="first")
mrna_5_geneid <- as.data.frame(mrna_5_geneid)

mrna_cds_vector <- as.vector(t(mrna_ID_cds))
mrna_cds_geneid <- mapIds(org.Hs.eg.db,keys=mrna_cds_vector,column="SYMBOL",keytype="ENSEMBL",multiVals="first")
mrna_cds_geneid <- as.data.frame(mrna_cds_geneid)

utr_3_gene <- cbind(mrna_3_geneid,utr_3)
utr_3_gene <- utr_3_gene[,c(1,3)]

utr_5_gene <- cbind(mrna_5_geneid,utr_5)
utr_5_gene <- utr_5_gene[,c(1,3)]

cds_gene <- cbind(mrna_cds_geneid,cds)
cds_gene <- cds_gene[,c(1,3)]

save(utr_3_gene,file = "utr_3_gene.RData")
save(utr_5_gene,file = "utr_5_gene.RData")
save(cds_gene, file = "cds_gene.RData")

cds_gene <- cds_gene[order(cds_gene$mrna_cds_geneid),]
utr_5_gene <- utr_5_gene[order(utr_5_gene$mrna_5_geneid),]
utr_3_gene <- utr_3_gene[order(utr_3_gene$mrna_3_geneid),]

gene_seq <- cbind(utr_3_gene,utr_5_gene)
gene_seq <- cbind(gene_seq,cds_gene)
colnames(gene_seq) <- c("gene","utr3","gene","utr5","gene","cds")
gene_seq <- gene_seq[,c(1,2,4,6)]
save(gene_seq,file = "gene_seq.RData")

write.csv(gene_seq,file = "gene_seq.csv")


for(i in 1:nrow(unique_trf)){
  if(unique_trf$trf_seq[i] == "GCCCGGCTAGCTCAGTCGGTAGAGC"){
    print(unique_trf$gene[i])
  }
}



#------------正则替换ID号--------------
#处理三个表格文件
rm(list = ls())
library(openxlsx)
gene_id1 <- read.xlsx("tRF_data/gene_info.xlsx", sheet = 1)
gene_id2 <- read.xlsx("tRF_data/gene_info.xlsx", sheet = 2)
gene_id3 <- read.xlsx("tRF_data/gene_info.xlsx", sheet = 3)

#gsub函数的应用
new_ID1 <- as.data.frame(gene_id1[,5])
for(i in 1:nrow(new_ID1)){
  new_ID1[i,] <- gsub(pattern = "^.*\\(",replacement = "",x = new_ID1[i,])
  new_ID1[i,] <- gsub(pattern = "\\)",replacement = "",x = new_ID1[i,])
}

new_ID1 <- cbind(gene_id1[,1],new_ID1)
colnames(new_ID1) <- c("trf5","gene")
write.csv(new_ID1,file = "new_ID5.csv")


new_ID2 <- as.data.frame(gene_id2[,5])
for(i in 1:nrow(new_ID2)){
  new_ID2[i,] <- gsub(pattern = "^.*\\(",replacement = "",x = new_ID2[i,])
  new_ID2[i,] <- gsub(pattern = "\\)",replacement = "",x = new_ID2[i,])
}

new_ID2 <- cbind(gene_id2[,1],new_ID2)
colnames(new_ID2) <- c("trf3","gene")
write.csv(new_ID2,file = "new_ID3.csv")


new_ID3 <- as.data.frame(gene_id3[,5])
for(i in 1:nrow(new_ID3)){
  new_ID3[i,] <- gsub(pattern = "^.*\\(",replacement = "",x = new_ID3[i,])
  new_ID3[i,] <- gsub(pattern = "\\)",replacement = "",x = new_ID3[i,])
}

new_ID3 <- cbind(gene_id3[,1],new_ID3)
colnames(new_ID3) <- c("trf1","gene")
write.csv(new_ID3,file = "new_ID1.csv")

new_ID1 <- read.csv(file = "tRF_data/trf/new_ID1.csv")
new_ID3 <- read.csv(file = "tRF_data/trf/new_ID3.csv")
new_ID5 <- read.csv(file = "tRF_data/trf/new_ID5.csv")

new_ID1 <- new_ID1[,c(2,3)]
new_ID3 <- new_ID3[,c(2,3)]
new_ID5 <- new_ID5[,c(2,3)]
#------------gene-trf配对---------
#加载基因数据
load("gene_seq.RData")
load("unique_trf.RData")

gene_index <- !is.na(gene_seq$utr3)
new_gene <- gene_seq[gene_index,]

gene_index <- !is.na(new_gene$utr5)
new_gene <- new_gene[gene_index,]

gene_index <- !is.na(new_gene$cds)
new_gene <- new_gene[gene_index,]

unique_gene <- new_gene[!duplicated(new_gene[,c(1)]),]
save(unique_gene,file = "unique_gene.RData")
write.csv(unique_gene,file = "unique_gene.csv")

trf1_gene <- as.data.frame(unique_trf$gene)
trf1_test <- unique_trf
pair_gene <- as.data.frame(matrix(data = "unavaliable",nrow = nrow(unique_trf),ncol = 3))
trf1_test <- cbind(trf1_test,pair_gene)
colnames(trf1_test) <- c("trf_seq","gene","seed","comp_seed","trf_type","utr3","utr5","cds")
trf1_gene[is.na(trf1_gene)] <- 0
unique_gene$gene[is.na(unique_gene$gene)] <- 0
for(i in 1:nrow(trf1_gene)){
  for(j in 1:nrow(unique_gene)){
    if(trf1_gene[i,1] == unique_gene$gene[j]){
      trf1_test$utr3[i] <- unique_gene$utr3[j]
      trf1_test$utr5[i] <- unique_gene$utr5[j]
      trf1_test$cds[i] <- unique_gene$cds[j]
      }
    }
  }

new_trf_feature <- trf1_test

for(i in 1:nrow(new_trf_feature)){
  if(new_trf_feature$utr3[i] == "unavaliable"){
    new_trf_feature$utr3[i] = NA
  }
}

index <- !is.na(new_trf_feature$utr3)
new_trf_feature <- new_trf_feature[index,]
save(new_trf_feature,file = "new_trf_features.RData")
write.csv(new_trf_feature, file = "new_trf_features.csv")
trf_features <- trf1_test

save(trf_features,file = "trf_features.RData")
write.csv(trf_features, file = "trf_features.csv")


unique_gene <- new_gene


#----------------数据预处理-------------------

#rm(list = ls())

#测试使用rbind()函数合并文件
trf_sequence1 <- read.table("tRF_data/tRF_sequence/A5HB.txt", header = FALSE, sep = "\t")
trf_sequence2 <- read.table("tRF_data/tRF_sequence/A3E6.txt", header = FALSE, sep = "\t")
rbind_data <- rbind(trf_sequence1, trf_sequence2)

#测试使用merge()函数合并文件
all_seq <- merge(trf_sequence1,trf_sequence2, by = c("MINTbase.Unique.ID", "tRF.sequence", "tRF.type.s."),
                 all = TRUE,sort = TRUE)

#使用dir()函数获取文件夹中所有文件
list_name <- dir("./tRF_data/tRF_sequence",pattern = ".txt")

merge_data <- read.table(file = paste("tRF_data/tRF_sequence",list_name[1],sep = "/",collapse = ","),header = TRUE,sep = "\t")
#利用for循环依次读取
for(i in 2:length(list_name)){
  new_data <- read.table(file = paste("tRF_data/tRF_sequence",list_name[i],sep = "/",collapse = ","),header =TRUE,sep = "\t")
  merge_data <- merge(merge_data, new_data, by = c("MINTbase.Unique.ID", "tRF.sequence", "tRF.type.s."),
                      all = TRUE,sort = TRUE)
}
save(merge_data, file = "trf_merge_data.RData")

trf_info <- merge_data[,c(1:3)]

save(trf_info, file = "trf_info.RData")
#--------------gene数据汇总------------------
#将所有序列数据进行汇总

#1.读取基因序列数据
rm(list = ls())
cds <- read.table("tRF_data/gene_info/cds.txt", header = TRUE, sep = "\n")
utr_3 <- read.table("tRF_data/gene_info/utr_3.txt", header = TRUE, sep = "\n")
utr_5 <- read.table("tRF_data/gene_info/utr_5.txt", header = TRUE, sep = "\n")
#2.保存成RData数据格式
sava(cds, file = "cds.RData")
save(utr_3, file = "utr_3.RData")
save(utr_5, file = "utr_5.RData")

#3.

 
#--------------gene ID转换--------------------

library(stats4)
library(BiocGenerics)
library(parallel)

library("AnnotationDbi")

library("org.Hs.eg.db")

# for test

GENEID <- c(1,2,3,9,10)

ENSEMBL <- c("ENSG00000000003","ENSG00000075218","ENSG00000256069","ENSG00000171428","ENSG00000156006")

df <- data.frame(ENSEMBL,GENEID)

# ENSG转换Symbol

df$symbol <- mapIds(org.Hs.eg.db,keys=ENSEMBL,column="SYMBOL",keytype="ENSEMBL",multiVals="first")

is.vector(mrn3_vector)



#---------------计算特征-------------
#加载数据
trf1 <- new_ID1
trf3 <- new_ID3
trf5 <- new_ID5

# 计算7mer-m1 seed
library(stringr)
library(stringi)
comp_seed_list <- matrix(data = "",nrow = nrow(trf1),ncol = 2)
trf_index = 0
for(j in 1:nrow(trf1)){
  trf_index = trf_index + 1
  trf1_seed <- substr(trf1$trf1[trf_index],1,7)
  comp_seed_list[j,1] <- trf1_seed
  comp_seed <- ""
  for(i in 1:nchar(trf1_seed)) {
    if(substr(trf1_seed, i, i) == "C") {
      comp_seed <- paste(comp_seed, "G", sep = "")
    }
    else if(substr(trf1_seed, i, i) == "G"){
      comp_seed <- paste(comp_seed, "C", sep = "")
    }
    else if(substr(trf1_seed, i, i) == "T"){
      comp_seed <- paste(comp_seed, "A", sep = "")
    }
    else if(substr(trf1_seed, i, i) == "A"){
      comp_seed <- paste(comp_seed, "T", sep = "")
    }
  }
  comp_seed_list[j,2] <- comp_seed
}
trf1 <- cbind(trf1,comp_seed_list)
colnames(trf1) <- c("trf1_seq","gene","seed","comp_seed")

trf1_type <- matrix(data = 1,nrow = nrow(trf1),ncol = 1)
trf1 <- cbind(trf1,trf1_type)
colnames(trf1) <- c("trf_seq","gene","seed","comp_seed","trf_type")

#trf3
comp_seed_list <- matrix(data = "",nrow = nrow(trf3),ncol = 2)
trf_index = 0
for(j in 1:nrow(trf3)){
  trf_index = trf_index + 1
  trf_seed <- substr(trf3$trf3[trf_index],1,7)
  comp_seed_list[j,1] <- trf_seed
  comp_seed <- ""
  for(i in 1:nchar(trf_seed)) {
    if(substr(trf_seed, i, i) == "C") {
      comp_seed <- paste(comp_seed, "G", sep = "")
    }
    else if(substr(trf_seed, i, i) == "G"){
      comp_seed <- paste(comp_seed, "C", sep = "")
    }
    else if(substr(trf_seed, i, i) == "T"){
      comp_seed <- paste(comp_seed, "A", sep = "")
    }
    else if(substr(trf_seed, i, i) == "A"){
      comp_seed <- paste(comp_seed, "T", sep = "")
    }
  }
  comp_seed_list[j,2] <- comp_seed
}
trf3 <- cbind(trf3,comp_seed_list)
colnames(trf3) <- c("trf3_seq","gene","seed","comp_seed")

trf3_type <- matrix(data = 3,nrow = nrow(trf3),ncol = 1)
trf3 <- cbind(trf3,trf3_type)
colnames(trf3) <- c("trf_seq","gene","seed","comp_seed","trf_type")

#trf5
comp_seed_list <- matrix(data = "",nrow = nrow(trf5),ncol = 2)
trf_index = 0
for(j in 1:nrow(trf5)){
  trf_index = trf_index + 1
  trf_seed <- substr(trf5$trf5[trf_index],1,7)
  comp_seed_list[j,1] <- trf_seed
  comp_seed <- ""
  for(i in 1:nchar(trf_seed)) {
    if(substr(trf_seed, i, i) == "C") {
      comp_seed <- paste(comp_seed, "G", sep = "")
    }
    else if(substr(trf_seed, i, i) == "G"){
      comp_seed <- paste(comp_seed, "C", sep = "")
    }
    else if(substr(trf_seed, i, i) == "T"){
      comp_seed <- paste(comp_seed, "A", sep = "")
    }
    else if(substr(trf_seed, i, i) == "A"){
      comp_seed <- paste(comp_seed, "T", sep = "")
    }
  }
  comp_seed_list[j,2] <- comp_seed
}
trf5 <- cbind(trf5,comp_seed_list)
colnames(trf5) <- c("trf5_seq","gene","seed","comp_seed")

trf5_type <- matrix(data = 5,nrow = nrow(trf5),ncol = 1)
trf5 <- cbind(trf5,trf5_type)
colnames(trf5) <- c("trf_seq","gene","seed","comp_seed","trf_type")

save(trf1,file = "trf1.RData")
save(trf3,file = "trf3.RData")
save(trf5,file = "trf5.RData")

load("trf1.RData")
load("trf3.RData")
load("trf5.RData")


all_trf <- rbind(trf1,trf3,trf5)
save(all_trf,file = "trf_info.RData")
write.csv(all_trf,file = "trf_info.csv")

unique_trf <- all_trf[!duplicated(all_trf[,c(1,2)]),]
unique_trf <- unique_trf[order(unique_trf$gene),]
save(unique_trf,file = "unique_trf.RData")
write.csv(unique_trf,file = "unique_trf.csv")


# AU content
load("new_trf_features.RData")

AU_content <- matrix(data = 0,nrow = nrow(new_trf_feature),ncol = 1)
new_trf_feature <- cbind(new_trf_feature,AU_content)
for(i in 1:nrow(new_trf_feature)){
  upstream <- substr(new_trf_feature$utr3[i],1,30)
  downstream <- substr(new_trf_feature$utr5[i],1,30)
  upscore <- 0
  downscore <- 0
  for(j in 1:30){
    if(substr(upstream, j, j) == "A" || substr(upstream, j, j) == "T") {
      upscore <- upscore + (1/(1 + j))
    }
    if(substr(downstream, j, j) == "A" || substr(downstream, j, j) == "T") {
      downscore <- downscore + (1/(1 + j))
    }
  }
  new_trf_feature$AU_content[i] <- mean(c(upscore,downscore))
}

save(new_trf_feature,file = "new_trf_features.RData")

#AU content

AU <- matrix(data = 0,nrow = nchar(trf_info),ncol = 2)
for(i in 1:nchar(trf_info)){
  for(k in 1:length(features[,1])) {
    if(k %% 1000 == 0) {
      print(k)
    }
    site_start <- as.integer(substr(features$mrna_binding_loc[k], 1, str_locate(features$mrna_binding_loc[k], ",") - 1))
    site_end <- as.integer(substr(features$mrna_binding_loc[k], str_locate(features$mrna_binding_loc[k], ",") + 1, nchar(features$mrna_binding_loc[k])))
    if(site_start <= 30) {
      downstream <- substr(features$mrna_sequence[k], 1, site_start - 1)
    } else {
      downstream <- substr(features$mrna_sequence[k], site_start - 30, site_start - 1)
    }
    if(site_end >= (nchar(features$mrna_sequence[k]) - 30)) {
      upstream <- substr(features$mrna_sequence[k], site_end + 1, nchar(features$mrna_sequence[k]))
    } else {
      upstream <- substr(features$mrna_sequence[k], site_end + 1, site_end + 30)
    }
    downstream <- stri_reverse(downstream)
    #初始化
    upscore <- 0
    downscore <- 0
    
    
    for(i in 1:30) {
      if(substr(upstream, i, i) == "A" || substr(upstream, i, i) == "T") {
        upscore <- upscore + (1/(i + 1))
      }
      if(substr(downstream, i, i) == "A" || substr(downstream, i, i) == "T") {
        downscore <- downscore + (1/(i + 1))
      }
    }
    features$au_content[k] <- mean(c(upscore, downscore))
  }
  AU[i,1] <- upscore
  AU[i,2] <- downscore
}


# binding site position(统计结合区的匹配位点数)


#定义置换函数
isPaired <- function(nt) {
  if(nt == "A") {
    return("T")
  }
  else if(nt == "T") {
    return("A")
  }
  else if(nt == "G") {
    return("C")
  }
  else if(nt == "C") {
    return("G")
  }
  else {
    return("")
  }
}


#确定结合区的位置
start_loc <- matrix(data = 0,nrow = nrow(new_trf_feature),ncol = 1)
end_loc <- matrix(data = 0,nrow = nrow(new_trf_feature),ncol = 1)
new_trf_feature <- cbind(new_trf_feature,start_loc)
new_trf_feature <- cbind(new_trf_feature,end_loc)
for(i in 1:nrow(new_trf_feature)){
  sequence <- substr(new_trf_feature$trf_seq[i],1,1)
  for(j in 1:nchar(new_trf_feature$utr3[i])){
    if(sequence == isPaired(substr(new_trf_feature$utr3[i],j,j))){
      new_trf_feature$start_loc[i] <- j
      new_trf_feature$end_loc[i] <- j + nchar(new_trf_feature$trf_seq[i])
      break
    }
  }
}
save(new_trf_feature,file = "new_trf_features.RData")

#统计匹配位点数
count_paired <- matrix(data = 0,nrow = nrow(new_trf_feature),ncol = 1)
new_trf_feature <- cbind(new_trf_feature,count_paired)
for(i in 1:nrow(new_trf_feature)){
  for(j in 1:nchar(new_trf_feature$trf_seq[i])){
    if(isPaired(substr(new_trf_feature$trf_seq[i],j,j)) == substr(new_trf_feature$utr3[i],
                    new_trf_feature$start_loc[i] + j - 1,
                    new_trf_feature$start_loc[i] + j - 1)){
      new_trf_feature$count_paired[i] = new_trf_feature$count_paired[i] + 1
    }
  }
}

save(new_trf_feature,file = "new_trf_features.RData")



# binding site position

#定义函数
isPaired <- function(nt) {
  if(nt == "A") {
    return("T")
  }
  else if(nt == "T") {
    return("A")
  }
  else if(nt == "G") {
    return("C")
  }
  else if(nt == "C") {
    return("G")
  }
  else {
    return("")
  }
}

for(k in 1:length(features[,1])) {
  if(k %% 1000 == 0) {
    print(k)
  }
  trf_start <- as.integer(substr(features$trf_binding_loc[k], 1, str_locate(features$trf_binding_loc[k], ",") - 1)) 
  trf_end <- as.integer(substr(features$trf_binding_loc[k], str_locate(features$trf_binding_loc[k], ",") + 1, nchar(features$trf_binding_loc[k])))
  mrna_start <- as.integer(substr(features$mrna_binding_loc[k], 1, str_locate(features$mrna_binding_loc[k], ",") - 1))
  mrna_end <- as.integer(substr(features$mrna_binding_loc[k], str_locate(features$mrna_binding_loc[k], ",") + 1, nchar(features$mrna_binding_loc[k])))
  #统计匹配位点数
  count_paired <- 0
  for(i in 1:min(mrna_end-mrna_start + 1, trf_end-trf_start + 1)){
    if(isPaired(substr(features$trf_sequence[k], trf_start + i - 1, trf_start + i - 1)) == substr(features$mrna_sequence[k], mrna_start + i - 1, mrna_start + i - 1)) {
      count_paired <- count_paired + 1
    }
  }
  features$num_paired_pos[k] <- count_paired
}

#结合区的最长区域

max_length <- matrix(data = 0,nrow = nrow(new_trf_feature),ncol = 1)
pos_longest <- matrix(data = 0,nrow = nrow(new_trf_feature),ncol = 1)
new_trf_feature <- cbind(new_trf_feature,max_length)
new_trf_feature <- cbind(new_trf_feature,pos_longest)
for(i in 1:nrow(new_trf_feature)){
  paired_info <- c()
  for(j in 1:nchar(new_trf_feature$trf_seq[i])){
    if(isPaired(substr(new_trf_feature$trf_seq[i],j,j)) == substr(new_trf_feature$utr3[i],
                                                                  new_trf_feature$start_loc[i] + j - 1,
                                                                  new_trf_feature$start_loc[i] + j - 1)){
      paired_info <- c(paired_info,1)
    }
    else{
      paired_info <- c(paired_info,0)
    }
  }
  count_consec <- 0
  pos_consec <- 0
  for(k in 1:length(rle(paired_info)$values)) {
    if(rle(paired_info)$values[k] == 1 && rle(paired_info)$lengths[k] > count_consec) {
        count_consec <- rle(paired_info)$lengths[k]
        pos_consec <- sum(rle(paired_info)$lengths[1:k])
    }
  }
  new_trf_feature$max_length[i] <- count_consec
  new_trf_feature$pos_longest[i] <- pos_consec
}
save(new_trf_feature,file = "new_trf_features.RData")



# the length of binding reign(结合区的长度)

for(k in 1:length(features[,1])) {
  if(k %% 1000 == 0) {
    print(k)
  }
  mrna_start <- as.integer(substr(features$mrna_binding_loc[k], 1, str_locate(features$mrna_binding_loc[k], ",") - 1))
  mrna_end <- as.integer(substr(features$mrna_binding_loc[k], str_locate(features$mrna_binding_loc[k], ",") + 1, nchar(features$mrna_binding_loc[k])))
  features$binding_region_length[k] <- mrna_end - mrna_start + 1
}



#最长匹配长度(结合区中最长的长度)
for(k in 1:length(features[,1])) {
  if(k %% 1000 == 0) {
    print(k)
  }
  trf_start <- as.integer(substr(features$trf_binding_loc[k], 1, str_locate(features$trf_binding_loc[k], ",") - 1)) 
  trf_end <- as.integer(substr(features$trf_binding_loc[k], str_locate(features$trf_binding_loc[k], ",") + 1, nchar(features$trf_binding_loc[k])))
  mrna_start <- as.integer(substr(features$mrna_binding_loc[k], 1, str_locate(features$mrna_binding_loc[k], ",") - 1))
  mrna_end <- as.integer(substr(features$mrna_binding_loc[k], str_locate(features$mrna_binding_loc[k], ",") + 1, nchar(features$mrna_binding_loc[k])))
  paired_info <- c()
  for(i in 1:min(mrna_end-mrna_start + 1, trf_end-trf_start + 1)){
    if(isPaired(substr(features$trf_sequence[k], trf_start + i - 1, trf_start + i - 1)) == substr(features$mrna_sequence[k], mrna_start + i - 1, mrna_start + i - 1)) {
      paired_info <- c(paired_info, 1)
    }
    else {
      paired_info <- c(paired_info, 0)
    }
  }
  count_consec <- 0
  pos_consec <- 0
  for(i in 1:length(rle(paired_info)$values)) {
    if(rle(paired_info)$values[i] == 1) {
      if(rle(paired_info)$lengths[i] > count_consec) {
        count_consec <- rle(paired_info)$lengths[i]
        pos_consec <- sum(rle(paired_info)$lengths[1:i])
      }
    }
  }
  features$longest_consecutive[k] <- count_consec
  features$pos_longest_consecutive[k] <- pos_consec + trf_start - 1
}

#3'区匹配位点数
three_prime <- matrix(data = 0,nrow = nrow(new_trf_feature),ncol = 1)
new_trf_feature <- cbind(new_trf_feature,three_prime)

for(i in 1:nrow(new_trf_feature)){
  three_prime_length <- 0
  for(j in 1:7){
    if(isPaired(substr(new_trf_feature$trf_seq[i],j,j)) == substr(new_trf_feature$utr3[i],
                                                                  new_trf_feature$start_loc[i] + j - 1,
                                                                  new_trf_feature$start_loc[i] + j - 1)){
      three_prime_length <- three_prime_length + 1
    }
  }
  new_trf_feature$three_prime[i] <- three_prime_length
}
save(new_trf_feature,file = "new_trf_features.RData")

#3'区匹配位点数
for(i in 1:length(features[,1])) {
  if(i %% 1000 == 0) {
    print(i)
  }
  
  trf_start <- as.integer(substr(features$trf_binding_loc[i], 1, str_locate(features$trf_binding_loc[i], ",") - 1)) 
  trf_end <- as.integer(substr(features$trf_binding_loc[i], str_locate(features$trf_binding_loc[i], ",") + 1, nchar(features$trf_binding_loc[i])))
  three_prime_start <- trf_end - 7
  mrna_start <- as.integer(substr(features$mrna_binding_loc[i], 1, str_locate(features$mrna_binding_loc[i], ",") - 1)) + ((trf_end - trf_start) - 7)
  mrna_end <- as.integer(substr(features$mrna_binding_loc[i], str_locate(features$mrna_binding_loc[i], ",") + 1, nchar(features$mrna_binding_loc[i])))
  count_paired <- 0
  for(i in 1:min(mrna_end-mrna_start + 1, 7)){
    if(isPaired(substr(features$trf_sequence[i], trf_start + i - 1, trf_start + i - 1)) == substr(features$mrna_sequence[i], mrna_start + i - 1, mrna_start + i - 1)) {
      count_paired <- count_paired + 1
    }
  }
  features$three_prime_pairs[i] <- count_paired
}


#seed 和trf在3'区的不配位数

for(k in 1:length(features[,1])) {
  if(k %% 1000 == 0) {
    print(k)
  }
  if(features$seed[k] == 1) {
    features$seed_end_diff[k] <- abs(features$three_prime_pairs[k] - 6)
  } else {
    trf_seed <- substr(features$trf_sequence[k], 2, 7)
    comp_seed <- ""
    for(i in 1:nchar(trf_seed)) {
      comp_seed <- paste(comp_seed, isPaired(substr(trf_seed, i, i)), sep = "")
    }
    mrna_start <- as.integer(substr(features$mrna_binding_loc[k], 1, str_locate(features$mrna_binding_loc[k], ",") - 1))
    mrna_comp_seed <- substr(features$mrna_sequence[k], mrna_start + 1, mrna_start + 6)
    count_seed <- 0
    for(i in 1:nchar(trf_seed)) {
      if(substr(comp_seed, i, i) == substr(mrna_comp_seed, i, i)) {
        count_seed <- count_seed + 1
      }
    }
    features$seed_end_diff[k] <- abs(features$three_prime_pairs[k] - count_seed)
  }
}



#phylop scores

for(k in 1:length(features[,1])) {
  mrna_start <- as.integer(substr(features$mrna_binding_loc[k], 1, str_locate(features$mrna_binding_loc[k], ",") - 1))
  mrna_end <- as.integer(substr(features$mrna_binding_loc[k], str_locate(features$mrna_binding_loc[k], ",") + 1, nchar(features$mrna_binding_loc[k])))
  cat(paste("chr", as.character(features$chr[k]), sep = ""))
  cat("\t")
  if(!is.na(suppressWarnings(as.integer(features$start_loc[k])))) {
    cat(as.character(as.integer(features$start_loc[k]) + mrna_start - 1))
    cat("\t")
    cat(as.character(as.integer(features$end_loc[k]) - (features$length_utr[k] - mrna_end)))
    cat("\t")
    cat(features$tran_id[k])
    cat("\n")
  } else {
    # Get piece-wise 3' UTR chromosomal coordinates
    x = unlist(str_locate_all(features$start_loc[k], ";"))
    x = x[1:(length(x)/2)]
    y = unlist(str_locate_all(features$end_loc[k], ";"))
    y = y[1:(length(y)/2)]
    all_start_locs = c(as.integer(substr(features$start_loc[k], 1, x[1] - 1)))
    if(length(x)-1 > 0) {
      for(i in 1:(length(x)-1)) {
        all_start_locs = c(all_start_locs, as.integer(substr(features$start_loc[k], x[i] + 1, x[i + 1] - 1)))
      }
    }
    all_start_locs = c(all_start_locs, as.integer(substr(features$start_loc[k], x[length(x)] + 1, nchar(features$start_loc[k]))))
    all_end_locs = c(as.integer(substr(features$end_loc[k], 1, x[1] - 1)))
    if(length(y)-1 > 0) {
      for(i in 1:(length(x)-1)) {
        all_end_locs = c(all_end_locs, as.integer(substr(features$end_loc[k], x[i] + 1, x[i + 1] - 1)))
      }
    }
    all_end_locs = c(all_end_locs, as.integer(substr(features$end_loc[k], x[length(x)] + 1, nchar(features$end_loc[k]))))
    all_start_locs = sort(all_start_locs)
    all_end_locs = sort(all_end_locs)
    differences = all_end_locs - all_start_locs + 1
    for(i in 2:length(differences)) {
      differences[i] = sum(differences[(i-1):i])
    }
    if(sum(mrna_start > differences) == 0) {
      cat(as.character(all_start_locs[1] + (mrna_start - 1)))
      cat("\t")
    } else{ 
      which_start_loc = max(which(mrna_start > differences))
      cat(as.character(all_start_locs[which_start_loc + 1] + (mrna_start-differences[which_start_loc] - 1)))
      cat("\t")
    }
    if(sum(mrna_end > differences) == 0) {
      cat(as.character(all_start_locs[1] + (mrna_start - 1) + (mrna_end - mrna_start)))
      cat("\t")
    } else {
      which_end_loc = max(which(mrna_end > differences))
      cat(as.character(all_start_locs[which_end_loc + 1] + (mrna_end - differences[which_end_loc] - 1)))
      cat("\t")  
    }
    cat(features$tran_id[k])
    cat("\n")
  }
}





