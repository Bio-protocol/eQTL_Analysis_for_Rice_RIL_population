source('Manhattan_function.R')
data <- read.table(" Ehd1_eQTL_result.txt ")
data <-data[,c(4,5,6,7)]
colnames(data)<-c('SNP','CHR','BP','P')
manhattan(data,main = "Manhattan Plot",suggestiveline =-log10(6.17e-04),genomewideline = -log10(6.12e-04),cex.lab=1.2, cex.axis=1.5,cex.main=1.5,annotatePval = 5e-40,annotateTop = FALSE,xlim=c(32284500,155000000),ylim=c(0,50),cex = 1.5)

