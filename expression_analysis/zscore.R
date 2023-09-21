suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(distributions3)))

args = commandArgs(trailingOnly=TRUE)

TPM_matrix <- args[1]
FPKM_matrix <- args[2]
TMM_matrix <- args[3]
assembly <- args[4]
tissue <- args[5]
histone <- args[6]
effect <- args[7]
TE <- args[8]
gene <- args[9]
position <- args[10]
distance <- args[11]

TMM_matrix <- fread(TMM_matrix)

if ( sum(dim(TMM_matrix)) == 0 ) {

print(paste(TE,gene,histone,effect,tissue,assembly,position,distance,"NA","NA","NA"))
write(paste(TE,gene,histone,effect,tissue,assembly,position,distance,"NA","NA","NA"),file=paste0("zscoreTMM/",assembly,"/",tissue,"/",histone,"/",effect,"/zscore_",histone,"-",tissue,"-",effect,"-",assembly,".tab"),append=TRUE)

} else {


colnames(TMM_matrix) <- c("gene","type","TMM","sample","tissue","strain","TE")

zscore<-(mean(TMM_matrix[TMM_matrix$TE=="TE",]$TMM)-mean(TMM_matrix[TMM_matrix$TE=="noTE",]$TMM))/(sqrt((sd(TMM_matrix[TMM_matrix$TE=="TE",]$TMM)^2/length(TMM_matrix[TMM_matrix$TE=="TE",]$sample))+(sd(TMM_matrix[TMM_matrix$TE=="noTE",]$TMM)^2/length(TMM_matrix[TMM_matrix$TE=="noTE",]$sample))))
Z<-Normal(0,1)
pval<-1 - cdf(Z, abs(zscore)) + cdf(Z, -abs(zscore))
fc <- log2(mean(TMM_matrix[TMM_matrix$TE=="TE",]$TMM)/mean(TMM_matrix[TMM_matrix$TE=="noTE",]$TMM)+0.000000000000000000000001)
print(paste(TE,gene,histone,effect,tissue,assembly,position,distance,round(zscore,3),signif(pval,3),round(fc,3)))
write(paste(TE,gene,histone,effect,tissue,assembly,position,distance,round(zscore,3),signif(pval,3),round(fc,3)),file=paste0("zscoreTMM/",assembly,"/",tissue,"/",histone,"/",effect,"/zscore_",histone,"-",tissue,"-",effect,"-",assembly,".tab"),append=TRUE)
}
