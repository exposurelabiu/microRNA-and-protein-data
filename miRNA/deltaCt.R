# https://support.bioconductor.org/p/p134303/#p134305
# https://support.bioconductor.org/p/105596/

pcr <- read.table("raw4Norm.tsv",sep="\t",quote="", row.names=1,header=TRUE, na.strings = "Undetermined")
assay_type <- pcr$Assay_type
pcr <- pcr[,-1]
# colMeans(pcr[grep("Ref",assay_type),])
# mean(pcr[grep("Ref",assay_type),1]
dCt <- apply(pcr,2,norm<-function(x){return (x-mean(x[grep("Ref",assay_type)]))})
dCtLimmaMatrix <- round(max(pcr,na.rm=TRUE)) - dCt
CtLimmaMatrix <- round(max(pcr,na.rm=TRUE)) - pcr
save(dCtLimmaMatrix,CtLimmaMatrix,file="forLimma.RData")
