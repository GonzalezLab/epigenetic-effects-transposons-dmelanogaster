suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(distributions3)))

args = commandArgs(trailingOnly=TRUE)


TMM_matrix <- args[1]
assembly <- args[2]
tissue <- args[3]
gene <- args[4]
position <- args[5]
distance <- args[6]


TMM_matrix <- fread(TMM_matrix)

if ( sum(dim(TMM_matrix)) == 0 ) {

print(paste(gene,tissue,assembly,"NA","NA","NA"))
write(paste(gene,tissue,assembly,"NA","NA","NA"),file=paste0("zscoreTMM_geneTE/",assembly,"/",tissue,"/geneTE2kb/zscore_",tissue,"-",assembly,".tab"),append=TRUE)

} else {


colnames(TMM_matrix) <- c("gene","type","TMM","sample","tissue","strain","TE")

zscore<-(mean(TMM_matrix[TMM_matrix$TE=="TE",]$TMM)-mean(TMM_matrix[TMM_matrix$TE=="noTE",]$TMM))/(sqrt((sd(TMM_matrix[TMM_matrix$TE=="TE",]$TMM)^2/length(TMM_matrix[TMM_matrix$TE=="TE",]$sample))+(sd(TMM_matrix[TMM_matrix$TE=="noTE",]$TMM)^2/length(TMM_matrix[TMM_matrix$TE=="noTE",]$sample))))
Z<-Normal(0,1)
pval<-1 - cdf(Z, abs(zscore)) + cdf(Z, -abs(zscore))
fc <- log2(mean(TMM_matrix[TMM_matrix$TE=="TE",]$TMM)/mean(TMM_matrix[TMM_matrix$TE=="noTE",]$TMM)+0.000000000000000000000001)
print(paste(gene,tissue,assembly,round(zscore,3),signif(pval,3),round(fc,3)))
write(paste(gene,tissue,assembly,round(zscore,3),signif(pval,3),round(fc,3)),file=paste0("zscoreTMM_geneTE/",assembly,"/",tissue,"/geneTE2kb/zscore_",tissue,"-",assembly,".tab"),append=TRUE)
}
