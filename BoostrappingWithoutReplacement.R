########################
#sn: The repeating times for each strapping
#matrix1: the read counts matrix of group 1, the rows represent genes and the columns represent samples
#matrix2: the read counts matrix of group 2, the rows represent genes and the columns represent samples
#group1Inf: the vector of the sample library.
#group2Inf: the vector of the sample library.
#outDir: The output directory.
#parallel_TF: to parallelly run DESeq2.
#ncores: the number of cores used for running DESeq2.
#paired_TF: to paired the samples in comparison or not.
RandomPooling = function(matrix1, matrix2, sn, group1Inf, group2Inf, outDir, parallel_TF, ncores, paired_TF){
  
  difgeneNum.df = c(0, 0)
  difgeneNum.Group1 = c(0, 0)
  difgeneNum.Group2 = c(0, 0)
  difgene.Group1.df = rep("x", nrow(matrix1))
  difgene.Group2.df = rep("x", nrow(matrix1))
  difgene.Group1.p = rep(2, nrow(matrix1))
  difgene.Group2.p = rep(2, nrow(matrix1))
  
  poolingVec = c(3 : (min(col(matrix1), col(matrix2)) - 1))
  
  for ( sampleNum in poolingVec) {
    for(r in (1:sn)) {
      
      index.spl = sample(1:ncol(matrix1),sampleNum)
      matrix1.tmp = matrix1[,index.spl]; 
      group1.tmp = data.frame(Type = c(group1Inf)[index.spl],Time = c(group1Inf)[index.spl])
      
      if ( paired_TF == FALSE ) {index.spl = sample(1:ncol(matrix2),sampleNum)}
      matrix2.tmp = matrix2[,index.spl]; 
      group2.tmp = data.frame(Type = c(group2Inf)[index.spl],Time = c(group2Inf)[index.spl])
      
      #QC
      matrix1.tmp.qc = myGeneExpFilter(rpkm = matrix1.tmp, rpkmCutOff = 1, expRatio = 0.33)
      matrix2.tmp.qc = myGeneExpFilter(rpkm = matrix2.tmp, rpkmCutOff = 1, expRatio = 0.33)
      
      read_r1 = cbind(matrix1.tmp, matrix2.tmp)
      sampleInf_r1 = c(group1.tmp, group2.tmp)
      
      read_r1.qc = read_r1[intersect(rownames(matrix1.tmp.qc), rownames(matrix2.tmp.qc)),]
      
      #DESeq2
      metaDataForD = data.frame(type = sampleInf_r1)
      
      res = myDESeq2(reads = read_r1.qc, GroupDf = metaDataForD, parallel_TF = parallel_TF, ncores = ncores)
      
      #resInf = res@elementMetadata
      
      difgene = as.data.frame(res@listData)
      row.names(difgene) = res@rownames
      
      difgene = difgene[!is.na(difgene$padj),]
      difgene = difgene[order(difgene$padj, decreasing = F),]
      
      difgene = difgene[difgene$padj<0.05,]
      
      difgene.Group1 = difgene[difgene$log2FoldChange<0,]
      difgene.Group2 = difgene[difgene$log2FoldChange>0,]
      
      difgeneNum.df = rbind(difgeneNum.df, c(sampleNum, nrow(difgene)))
      difgeneNum.Group1 = rbind(difgeneNum.Group1, c(sampleNum, nrow(difgene.Group1)))
      difgeneNum.Group2 = rbind(difgeneNum.Group2, c(sampleNum, nrow(difgene.Group2)))
      difgene.Group1.df = cbind(difgene.Group1.df, c(rownames(difgene.Group1), rep("x", nrow(matrix1)-nrow(difgene.Group1))))
      difgene.Group2.df = cbind(difgene.Group2.df, c(rownames(difgene.Group2), rep("x", nrow(matrix1)-nrow(difgene.Group2))))
      difgene.Group1.p = cbind(difgene.Group1.p, c((difgene.Group1$padj), rep(2, nrow(matrix1)-nrow(difgene.Group1))))
      difgene.Group2.p = cbind(difgene.Group2.p, c((difgene.Group2$padj), rep(2, nrow(matrix1)-nrow(difgene.Group2))))
    }
  }
  
  #Group1 DEGs
  deg.Group1 = difgene.Group1.df
  deg.Group1 = c(as.matrix(deg.Group1))
  deg.Group1 = deg.Group1[!deg.Group1 == "x"]
  write.table(x = table(deg.Group1), file = paste(outDir, "Group1.DEGs.Frequency.txt"), sep = "\t")
  
  #Group2 DEGs
  deg.Group2 = difgene.Group2.df
  deg.Group2 = c(as.matrix(deg.Group2))
  deg.Group2 = deg.Group2[!deg.Group2 == "x"]
  write.table(x = table(deg.Group2), file = paste(outDir, "Group2.DEGs.Frequency.txt"), sep = "\t")
  
}

######################
myGeneExpFilter = function(rpkm, rpkmCutOff, expRatio) {
  ExpGeneRatio = expRatio
  gene_expression_cell_num_1 <- NULL
  for (i in 1:nrow(rpkm)) {
    gene_expression_cell_num_1[i] <- sum(rpkm[i,]>rpkmCutOff)
  }
  gene_expression_cell_num_1_f <- (gene_expression_cell_num_1 >= ncol(rpkm)*ExpGeneRatio)#second filtering row for gene
  #table(gene_expression_cell_num_1_f)
  rpkm <- rpkm[gene_expression_cell_num_1_f,]
  return(rpkm)
}

######################
#DESeq2
myDESeq2 = function(reads, GroupDf, parallel_TF, ncores) {
  
  if(parallel_TF) {
    library("DESeq2")
    cts = reads
    coldata = GroupDf
    colnames(coldata) = c("Group")
    dds <- DESeqDataSetFromMatrix(countData = cts,
                                  colData = coldata,
                                  design = ~ Group)
    #dds
    library(BiocParallel)
    register(MulticoreParam(ncores))
    dds = DESeq(dds, parallel = parallel_TF)
    res = results(dds)
  } else {
    library("DESeq2")
    cts = reads
    coldata = GroupDf
    colnames(coldata) = c("Group")
    dds <- DESeqDataSetFromMatrix(countData = cts,
                                  colData = coldata,
                                  design = ~ Group)
    #dds
    dds = DESeq(dds)
    res = results(dds)
  }
  return(res)
}
