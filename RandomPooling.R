
########################
# sn: The duplicate times for each pooling
# matrix1: the read counts matrix of group 1, the rows represent genes and the columns represent samples
# matrix2: the read counts matrix of group 2, the rows represent genes and the columns represent samples
#group1Inf: the vector of the sample library.
#group2Inf: the vector of the sample library.
#outDir: The output directory.
myRandomPooling = function(matrix1, matrix2, sn, group1Inf, group2Inf, outDir){
  
  difgeneNum.df = c(0, 0)
  difgeneNum.SNc = c(0, 0)
  difgeneNum.VTA = c(0, 0)
  difgene.SNc.df = rep("x", nrow(matrix1))
  difgene.VTA.df = rep("x", nrow(matrix1))
  difgene.SNc.p = rep(2, nrow(matrix1))
  difgene.VTA.p = rep(2, nrow(matrix1))
  
  poolingVec = c(3:min(col(matrix1), col(matrix2)))
  
  for ( sampleNum in poolingVec) {
    for(r in (1:sn)) {
      
      index.spl = sample(1:ncol(matrix1),sampleNum)
      matrix1.tmp = matrix1[,index.spl]; 
      group1.tmp = data.frame(Type = c(group1Inf)[index.spl],Time = c(group1Inf)[index.spl])

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

      res = myDESeq2(reads = read_r1.qc, GroupDf = metaDataForD)
      
      #resInf = res@elementMetadata
      
      difgene = as.data.frame(res@listData)
      row.names(difgene) = res@rownames
      
      difgene = difgene[!is.na(difgene$padj),]
      difgene = difgene[order(difgene$padj, decreasing = F),]
      
      difgene = difgene[difgene$padj<0.05,]
      
      difgene.SNc = difgene[difgene$log2FoldChange<0,]
      difgene.VTA = difgene[difgene$log2FoldChange>0,]
      
      difgeneNum.df = rbind(difgeneNum.df, c(sampleNum, nrow(difgene)))
      difgeneNum.SNc = rbind(difgeneNum.SNc, c(sampleNum, nrow(difgene.SNc)))
      difgeneNum.VTA = rbind(difgeneNum.VTA, c(sampleNum, nrow(difgene.VTA)))
      difgene.SNc.df = cbind(difgene.SNc.df, c(rownames(difgene.SNc), rep("x", nrow(matrix1)-nrow(difgene.SNc))))
      difgene.VTA.df = cbind(difgene.VTA.df, c(rownames(difgene.VTA), rep("x", nrow(matrix1)-nrow(difgene.VTA))))
      difgene.SNc.p = cbind(difgene.SNc.p, c((difgene.SNc$padj), rep(2, nrow(matrix1)-nrow(difgene.SNc))))
      difgene.VTA.p = cbind(difgene.VTA.p, c((difgene.VTA$padj), rep(2, nrow(matrix1)-nrow(difgene.VTA))))
    }
  }
  
  #SNc DEGs
  deg.SNc = difgene.SNc.df
  deg.SNc = c(as.matrix(deg.SNc))
  deg.SNc = deg.SNc[!deg.SNc == "x"]
  write.table(x = table(deg.SNc), file = paste(outDir, "Group1.DEGs.Frequency.txt"), sep = "\t")
  
  #VTA DEGs
  deg.VTA = difgene.VTA.df
  deg.VTA = c(as.matrix(deg.VTA))
  deg.VTA = deg.VTA[!deg.VTA == "x"]
  write.table(x = table(deg.VTA), file = paste(outDir, "Group2.DEGs.Frequency.txt"), sep = "\t")
  
}

######################
myGeneExpFilter = function(rpkm, rpkmCutOff, expRatio) {
  rpkm_l1_order = rpkm
  ExpGeneRatio = expRatio
  gene_expression_cell_num_1 <- NULL
  for (i in 1:nrow(rpkm_l1_order)) {
    gene_expression_cell_num_1[i] <- sum(rpkm_l1_order[i,]>rpkmCutOff)
  }
  gene_expression_cell_num_1_f <- (gene_expression_cell_num_1 >= ncol(rpkm_l1_order)*ExpGeneRatio)#second filtering row for gene
  #table(gene_expression_cell_num_1_f)
  rpkm_l1_order <- rpkm_l1_order[gene_expression_cell_num_1_f,]
  return(rpkm_l1_order)
}

######################
#DESeq2
myDESeq2 = function(reads, GroupDf) {
  
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
  return(res)
}
