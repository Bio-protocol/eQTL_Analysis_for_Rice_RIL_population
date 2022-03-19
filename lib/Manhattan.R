library('qqman')
data <- read.table(" Ehd1_eQTL_result.txt ")
data <-data[,c(4,5,6,7)]
colnames(data)<-c('SNP','CHR','BP','P')
manhattan(data,main = "Manhattan Plot", suggestiveline =FALSE,genomewideline = FALSE,csi=1.5,cex.lab=1.5, cex.axis=1.5,cex.main=1.5,annotatePval = 5e-40, annotateTop = FALSE,xlim=c(32284500,155000000),ylim=c(0,50),cex = 1.5)
abline(h=-log10(6.17e-04), col="blue", lty=1, lwd=3) 
abline(h=-log10(6.12e-04), col="red", lty=2, lwd=3)

