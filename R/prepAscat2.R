preAscat2 = function(tumour_count_file = "tumor_ac.tsv", normal_count_file = "normal_ac.tsv", tumor_sample = "sample"){

  tumorcounts = read.table(tumour_count_file,sep="\t")
  normalcounts = read.table(normal_count_file,sep="\t")

  SNPpos = data.frame(row.names = paste0("S", 1:nrow(tumorcounts)), Chr = tumorcounts$V1, Position = tumorcounts$V2)

  Tumor_BAF = matrix(nrow = dim(normalcounts)[1],ncol = 1)
  rownames(Tumor_BAF) =
  colnames(Tumor_BAF) = tumor_sample

  acgt = tumorcounts[,c(3:6)]
  acgts = t(apply(acgt,1,sort))
  Tumor_BAF[,1] = acgts[,4]/(acgts[,3]+acgts[,4])
  Tumor_BAF[,1] = ifelse(runif(length(Tumor_BAF[,1]))<0.5,Tumor_BAF[,1],1-Tumor_BAF[,1])
  Tumor_BAF[is.nan(Tumor_BAF)]=NA

  Germline_BAF = matrix(nrow = dim(normalcounts)[1],ncol = 1)
  rownames(Germline_BAF) = paste0("S", 1:nrow(tumorcounts))
  colnames(Germline_BAF) = tumor_sample
  acgt = normalcounts[,c(3:6)]
  acgts = t(apply(acgt,1,sort))
  Germline_BAF[,1] = acgts[,4]/(acgts[,3]+acgts[,4])
  Germline_BAF[,1] = ifelse(runif(length(Germline_BAF[,1]))<0.5,Germline_BAF[,1],1-Germline_BAF[,1])
  Germline_BAF[is.nan(Germline_BAF)]=NA

  Tumor_LogR = matrix(nrow = dim(normalcounts)[1],ncol = 1)
  Germline_LogR = matrix(nrow = dim(normalcounts)[1],ncol = 1)
  rownames(Tumor_LogR) = paste0("S", 1:nrow(tumorcounts))
  colnames(Tumor_LogR) = tumor_sample
  rownames(Germline_LogR) = paste0("S", 1:nrow(tumorcounts))
  colnames(Germline_LogR) = tumor_sample
  Tumor_LogR[,1] = log(tumorcounts[,7]/normalcounts[,7],2)
  Germline_LogR[,1] = 0
  Tumor_LogR[! is.finite(Tumor_LogR)]=NA # infinite = coverage in normal only, NaN = no coverage in tumour or normal
  Tumor_LogR[,1] = Tumor_LogR[,1] - median(Tumor_LogR[,1],na.rm=T)

  write.table(cbind(SNPpos,Tumor_LogR),paste(tumor_sample,".tumour.LogR.txt",sep=""),sep="\t",row.names=T,col.names=NA,quote=F)
  write.table(cbind(SNPpos,Tumor_BAF),paste(tumor_sample,".tumour.BAF.txt",sep=""),sep="\t",row.names=T,col.names=NA,quote=F)
  write.table(cbind(SNPpos,Germline_LogR),paste(tumor_sample,".normal.LogR.txt",sep=""),sep="\t",row.names=T,col.names=NA,quote=F)
  write.table(cbind(SNPpos,Germline_BAF),paste(tumor_sample,".normal.BAF.txt",sep=""),sep="\t",row.names=T,col.names=NA,quote=F)

  ascat.bc = ascat.loadData(paste(tumor_sample,".tumour.LogR.txt",sep=""),paste(tumor_sample,".tumour.BAF.txt",sep=""),
                            paste(tumor_sample,".normal.LogR.txt",sep=""),paste(tumor_sample,".normal.BAF.txt",sep=""),
                            chrs=paste0("chr", c(1:22, "X", "Y")), sexchromosomes = paste0("chr", c("X","Y")))

  ascat.plotRawData(ascat.bc, img.prefix = "xx")
  ascat.bc = ascat.aspcf(ascat.bc, out.prefix = "xx")
  ascat.plotSegmentedData(ascat.bc, img.prefix = "xx")


}
